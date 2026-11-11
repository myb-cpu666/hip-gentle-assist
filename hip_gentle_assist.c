#include "hip_gentle_assist.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#define PI 3.14159265358979323846
#define TAU (2.0 * PI)

static double clip(double x, double lo, double hi)
{
    return fmin(hi, fmax(lo, x));
}

static double wrap(double x)
{
    const double y = fmod(x, TAU);
    return y < 0.0 ? y + TAU : y;
}

static bool parameters_valid(const hip_gentle_parameters_t *p)
{
    if (!p) return false;
    const double positive[] = {
        p->torque_limit_nm, p->start_angle_rad, p->derivative_window_s,
        p->startup_lift_angle_rad, p->startup_lift_speed_rad_s, p->startup_lift_confirm_s,
        p->min_frequency_hz, p->max_frequency_hz, p->quiet_speed_rad_s,
        p->quiet_duration_s, p->low_speed_release_confirm_s, p->hold_confirm_s, p->hold_range_rad,
        p->reversal_speed_rad_s, p->reversal_confirm_s,
        p->peak_drop_ratio, p->progress_timeout_s,
        p->tracking_error_rad, p->tracking_error_duration_s, p->phase_error_limit_rad,
        p->phase_correction_rad_s, p->slew_nm_s, p->reset_gap_s, p->width,
        p->motion_rise_s, p->motion_fall_s,
        p->inner_rise_s, p->outer_rise_s
    };
    const double nonnegative[] = {p->amplitude_nm, p->amplitude_gain, p->phase_gain,
                                  p->frequency_gain, p->bias_gain, p->extension_torque_ratio};
    for (size_t i = 0; i < sizeof(positive)/sizeof(positive[0]); ++i)
        if (!isfinite(positive[i]) || positive[i] <= 0.0) return false;
    for (size_t i = 0; i < sizeof(nonnegative)/sizeof(nonnegative[0]); ++i)
        if (!isfinite(nonnegative[i]) || nonnegative[i] < 0.0) return false;
    if (p->adaptive_enabled) {
        const double adaptive_positive[] = {p->adaptive_max_nm, p->swing_up_threshold_rad, p->swing_down_threshold_rad,
            p->amplitude_increase_nm_s, p->amplitude_decrease_nm_s};
        for (size_t i = 0; i < sizeof(adaptive_positive)/sizeof(adaptive_positive[0]); ++i)
            if (!isfinite(adaptive_positive[i]) || adaptive_positive[i] <= 0) return false;
        if (p->adaptive_max_nm < fmin(p->amplitude_nm, p->torque_limit_nm) ||
            p->swing_down_threshold_rad >= p->swing_up_threshold_rad) return false;
    }
    return isfinite(p->asymmetry) && p->asymmetry > -1.0 && p->asymmetry < 1.0 &&
           isfinite(p->response) && p->response >= 0.0 && p->response <= 100.0 &&
           p->peak_drop_ratio < 1.0 && p->weak_half_cycles >= 2 && p->extension_torque_ratio <= 1.0 &&
           p->reversal_speed_rad_s > p->quiet_speed_rad_s &&
           p->low_speed_release_confirm_s < p->quiet_duration_s &&
           p->reversal_confirm_s < p->quiet_duration_s &&
           p->startup_lift_confirm_s < p->progress_timeout_s &&
           p->hold_confirm_s > p->low_speed_release_confirm_s &&
           p->hold_confirm_s < p->quiet_duration_s && p->hold_confirm_s < p->progress_timeout_s &&
           p->hold_range_rad < p->start_angle_rad &&
           p->motion_fall_s <= .2 &&
           p->width >= 0.2 && p->width <= 2.0 && p->min_frequency_hz < p->max_frequency_hz &&
           p->reset_gap_s >= 0.01 && p->reset_gap_s <= 0.02 &&
           /* 毫秒时钟最小非零间隔 1ms，固定缓冲必须容纳整个窗口及边界样本。 */
           p->derivative_window_s <= (HIP_GENTLE_HISTORY_CAPACITY-3U)/1000.0;
}

void hip_gentle_reset(hip_gentle_assist_t *s)
{
    if (!s || !s->initialized) return;
    const hip_gentle_parameters_t p = s->parameters;
    memset(s, 0, sizeof(*s));
    s->parameters = p;
    s->initialized = true;
    s->state = HIP_GENTLE_WAIT;
    s->motion_ramp_pending = true;
    s->omega = PI;
    s->reason = HIP_GENTLE_WAITING_CROSSING;
}

bool hip_gentle_init(hip_gentle_assist_t *s, const hip_gentle_parameters_t *p)
{
    if (!s) return false;
    if (!parameters_valid(p)) { memset(s, 0, sizeof(*s)); return false; }
    s->parameters = *p;
    s->initialized = true;
    hip_gentle_reset(s);
    return true;
}

double hip_gentle_waveform(double phase, double c, double d)
{
    if (!isfinite(phase) || !isfinite(c) || !isfinite(d) ||
        c < 0.2 || c > 2.0 || d <= -1.0 || d >= 1.0) return NAN;
    const double cc = cos(PI/c);
    const double arg = c < 1.0
        ? (d >= 0.0 ? 2*d/(1+d)-1 : -2*d/(1-d)-1)
        : (d >= 0.0 ? (cc+d)/(1+d) : (cc-d)/(1-d));
    const double arc = asin(clip(arg, -1, 1));
    const double delta = d >= 0.0 ? c/2*(PI/2-arc)-PI/2 : c/2*(arc-PI/2)+PI/2;
    const double x = wrap(phase-PI/2+delta);
    const double l1 = (1+cos(2*(x-PI/2)/c))/2;
    const double l2 = (1+cos(2*(x-3*PI/2)/c))/2;
    if (c <= 1.0) {
        const double h = (1-c)*PI/2;
        if (x < h || x >= TAU-h) return -d;
        if (x >= PI-h && x < PI+h) return d;
    }
    double value;
    if (x < PI/2) value = -d-(1-d)*l1;
    else if (x < PI) value = d-(1+d)*l1;
    else if (x < 3*PI/2) value = d+(1-d)*l2;
    else value = -d+(1+d)*l2;
    if (c > 1.0) {
        const double k = (1+cc)/2;
        value = (value+(x < PI ? k : -k))/(1-k);
    }
    return value;
}

static double reference_peak(const hip_gentle_assist_t *s, int direction)
{
    const unsigned side = direction < 0 ? 1U : 0U;
    double reference = 0;
    for (unsigned i = 0; i < s->recent_peak_count[side]; ++i)
        reference = fmax(reference, s->recent_peaks[side][i]);
    return reference;
}

static void remember_peak(hip_gentle_assist_t *s, double peak, int direction)
{
    const unsigned side = direction < 0 ? 1U : 0U;
    s->last_peak = peak;
    s->recent_peaks[side][s->recent_peak_head[side]] = peak;
    s->recent_peak_head[side] = (s->recent_peak_head[side]+1U)%HIP_GENTLE_PEAK_HISTORY_CAPACITY;
    if (s->recent_peak_count[side] < HIP_GENTLE_PEAK_HISTORY_CAPACITY) ++s->recent_peak_count[side];
}

static void clear_peak_history(hip_gentle_assist_t *s)
{
    memset(s->recent_peak_head, 0, sizeof(s->recent_peak_head));
    memset(s->recent_peak_count, 0, sizeof(s->recent_peak_count));
    s->last_peak = 0;
}

static void seed(hip_gentle_assist_t *s, double amplitude, double phase, double omega)
{
    s->bias = 0;
    s->omega = omega;
    memset(s->amplitudes, 0, sizeof(s->amplitudes));
    memset(s->phases, 0, sizeof(s->phases));
    memset(s->rates, 0, sizeof(s->rates));
    s->amplitudes[0] = amplitude;
    s->phases[0] = wrap(phase);
    s->bootstrap_halves_remaining = 0;
    s->adaptive_high = false;
    s->adaptive_amplitude_nm = fmin(s->parameters.amplitude_nm, s->parameters.torque_limit_nm);
    clear_peak_history(s);
    remember_peak(s, amplitude, s->direction);
}

static bool update_adaptive_amplitude(hip_gentle_assist_t *s, double dt)
{
    const hip_gentle_parameters_t *p = &s->parameters;
    const double base = fmin(p->amplitude_nm, p->torque_limit_nm);
    if (!p->adaptive_enabled || base == 0) {
        s->adaptive_high = false;
        s->adaptive_amplitude_nm = p->amplitude_nm;
        return true;
    }
    const double swing = fabs(s->amplitudes[0]);
    if (!isfinite(swing)) return false;
    const bool ready = s->state == HIP_GENTLE_TRACK && !s->bootstrap_pending &&
                       !s->weak_pending && s->gate >= 1;
    if (!ready) s->adaptive_high = false;
    else if (s->adaptive_high) {
        if (swing <= p->swing_down_threshold_rad) s->adaptive_high = false;
    } else if (swing >= p->swing_up_threshold_rad) s->adaptive_high = true;
    const double goal = s->adaptive_high ? fmin(p->adaptive_max_nm, p->torque_limit_nm) : base;
    s->adaptive_amplitude_nm += clip(goal-s->adaptive_amplitude_nm,
        -p->amplitude_decrease_nm_s*dt, p->amplitude_increase_nm_s*dt);
    return isfinite(s->adaptive_amplitude_nm);
}

/* 每快照一次：积分旧变化率，再计算本拍误差/变化率。没有 500Hz 伪采样。 */
static bool track(hip_gentle_assist_t *s, double x, double dt)
{
    const hip_gentle_parameters_t *p = &s->parameters;
    s->bias += dt*s->rates[0];
    const double next_omega = s->omega+dt*s->rates[1];
    if (!isfinite(next_omega) || !isfinite(s->bias)) return false;
    s->omega = clip(next_omega, TAU*p->min_frequency_hz, TAU*p->max_frequency_hz);
    const double correction = clip(s->phase_correction, -p->phase_correction_rad_s*dt, p->phase_correction_rad_s*dt);
    s->phase_correction -= correction;
    for (unsigned k = 0; k < HIP_GENTLE_HARMONICS; ++k) {
        s->amplitudes[k] += dt*s->rates[2+k];
        s->phases[k] = wrap(s->phases[k]+dt*s->rates[7+k]+(k+1)*correction);
        if (!isfinite(s->amplitudes[k]) || !isfinite(s->phases[k])) return false;
    }
    const double a = sin(s->phases[0]), b = cos(s->phases[0]);
    const double basis = 1.5*a-0.5*a*a*a, derivative = 1.5*b*b*b;
    double harmonics = 0;
    for (unsigned k = 1; k < HIP_GENTLE_HARMONICS; ++k)
        harmonics += s->amplitudes[k]*sin(s->phases[k]);
    s->estimate = s->bias+s->amplitudes[0]*basis+harmonics;
    const double error = x-s->estimate;
    s->rates[0] = p->bias_gain*error;
    s->rates[1] = p->frequency_gain*error*derivative;
    for (unsigned k = 0; k < HIP_GENTLE_HARMONICS; ++k) {
        s->rates[2+k] = p->amplitude_gain*error*(k == 0 ? basis : 0.5*sin(s->phases[k]));
        s->rates[7+k] = (k+1)*s->omega+p->phase_gain*error*(k == 0 ? derivative : cos(s->phases[k]));
    }
    if (!isfinite(s->estimate)) return false;
    for (unsigned k = 0; k < 12; ++k) if (!isfinite(s->rates[k])) return false;
    return true;
}

static void apply_target(hip_gentle_assist_t *s, double target, double dt)
{
    const hip_gentle_parameters_t *p = &s->parameters;
    s->target = clip(target, -p->torque_limit_nm, p->torque_limit_nm);
    const double previous = s->torque;
    s->torque += clip(s->target-s->torque, -p->slew_nm_s*dt, p->slew_nm_s*dt);
    s->torque_rate_nm_s = dt > 0 ? (s->torque-previous)/dt : 0;
}

static double smooth_fraction(double u)
{
    u = clip(u, 0, 1);
    return u*u*(3-2*u);
}

static void release_motion_torque(hip_gentle_assist_t *s, double dt)
{
    if (!s->motion_releasing) {
        s->motion_release_torque = s->torque;
        s->motion_release_elapsed = 0;
        s->motion_release_duration = fmax(s->parameters.motion_fall_s,
            1.5*fabs(s->torque)/s->parameters.slew_nm_s);
        const double sign = s->torque > 0 ? 1 : -1;
        /* Hermite初始斜率限制在回零方向，归一化幅度不超过3，保证单调且不越零。 */
        const double rate_limit = fmin(s->parameters.slew_nm_s,
            3*fabs(s->torque)/s->motion_release_duration);
        s->motion_release_rate = -sign*clip(-sign*s->torque_rate_nm_s, 0, rate_limit);
        s->motion_releasing = true;
    }
    s->motion_release_elapsed += dt;
    const double u = clip(s->motion_release_elapsed/s->motion_release_duration, 0, 1);
    const double target = s->motion_release_torque*(1-smooth_fraction(u))+
        s->motion_release_duration*s->motion_release_rate*u*(1-u)*(1-u);
    apply_target(s, clip(target, fmin(0, s->motion_release_torque), fmax(0, s->motion_release_torque)), dt);
    s->target = 0;
    if (s->torque == 0) s->motion_releasing = false;
}

static void update_motion_permission(hip_gentle_assist_t *s, double speed, double dt)
{
    const hip_gentle_parameters_t *p = &s->parameters;
    const int direction = speed > p->reversal_speed_rad_s ? 1 :
                          speed < -p->reversal_speed_rad_s ? -1 : 0;
    s->motion_confirm_time = direction && direction == s->motion_direction
        ? s->motion_confirm_time+dt : direction ? dt : 0;
    s->motion_direction = direction;
    s->quiet_time = fabs(speed) < p->quiet_speed_rad_s ? s->quiet_time+dt : 0;
    if (fabs(speed) < p->quiet_speed_rad_s &&
        s->quiet_time+1e-12 >= p->low_speed_release_confirm_s) {
        /* 短暂低速不撤销许可，避免误回零后被完成方向锁挡住整个半波。
         * 卸力确认与正式停止共用连续低速计时，各自使用独立门限。 */
        s->motion_waiting = true;
        s->torque_direction = 0;
        s->motion_ramp_pending = true;
    }
}

static void update_startup_permission(hip_gentle_assist_t *s,
                                      double left, double right, double vl, double vr, double dt)
{
    if (s->startup_released) return;
    const hip_gentle_parameters_t *p = &s->parameters;
    const double angles[2] = {left, right}, velocities[2] = {vl, vr};
    bool lifted[2];
    for (unsigned i = 0; i < 2; ++i) {
        if (!s->startup_reference_valid || velocities[i] <= p->startup_lift_speed_rad_s) {
            s->startup_min_angle[i] = angles[i];
            s->startup_lift_time[i] = 0;
        } else {
            s->startup_min_angle[i] = fmin(s->startup_min_angle[i], angles[i]);
            s->startup_lift_time[i] += dt;
        }
        const double rise = angles[i]-s->startup_min_angle[i];
        lifted[i] = isfinite(rise) && rise >= p->startup_lift_angle_rad &&
            s->startup_lift_time[i]+1e-12 >= p->startup_lift_confirm_s;
    }
    s->startup_reference_valid = true;
    if (!s->startup_first_leg) {
        if (s->state == HIP_GENTLE_COLLECT || s->state == HIP_GENTLE_START ||
            s->state == HIP_GENTLE_TRACK) {
            const unsigned first = s->direction > 0 ? 0U : 1U;
            const unsigned other = 1U-first;
            /* 角度差只确定候选腿，必须有该腿实际屈曲；同时抬腿不计作交替。 */
            if (lifted[first] && velocities[other] <= 0) {
                s->startup_first_leg = (uint8_t)(first+1U);
                s->startup_min_angle[other] = angles[other];
                s->startup_lift_time[other] = 0;
            }
        }
        return;
    }
    const unsigned second = s->startup_first_leg == 1U ? 1U : 0U;
    if (lifted[second] && velocities[1U-second] <= 0) s->startup_second_confirmed = true;
    if ((s->state == HIP_GENTLE_START || s->state == HIP_GENTLE_TRACK) && s->startup_second_confirmed) {
        s->startup_released = true;
        s->enabled_time = s->time_s;
        s->restart_angle_rad = 0;
    }
}

static void apply_motion_target(hip_gentle_assist_t *s, double target, double dt)
{
    const int direction = target > 0 ? 1 : target < 0 ? -1 : 0;
    const bool confirmed = s->motion_direction == direction &&
        s->motion_confirm_time+1e-12 >= s->parameters.reversal_confirm_s;
    if (direction != s->requested_direction) {
        s->requested_direction = direction; s->direction_wait_s = 0;
    }
    if (direction && direction != s->torque_direction && s->torque == 0 && !confirmed)
        s->direction_wait_s += dt;
    const bool interrupted = s->motion_ramp_pending || s->motion_releasing;
    /* 获准但尚未出力就换相，也要结束该许可，不能依赖微小脉冲解锁下一方向。 */
    if (s->torque == 0 && s->torque_direction && direction && direction != s->torque_direction) {
        s->completed_direction = s->torque_direction;
        s->torque_direction = 0; s->motion_waiting = true;
    }
    if (s->paused_direction && direction && direction != s->paused_direction) {
        s->completed_direction = s->paused_direction;
        s->paused_direction = 0;
    }
    /* 正常相位段结束后，不能因波形反弹或重新测速确认而再次同向启动。 */
    if (direction && direction == s->completed_direction) target = 0;
    else if (direction && (direction != s->torque_direction || s->motion_waiting)) {
        const bool progressed = s->paused_direction != direction ||
            direction*(s->previous_x-s->paused_x) >= s->parameters.startup_lift_angle_rad;
        if (confirmed && progressed && s->torque == 0) {
            s->torque_direction = direction;
            s->paused_direction = 0;
            s->motion_waiting = false;
            const bool late_entry = fabs(target) > s->parameters.slew_nm_s*dt ||
                s->direction_wait_s+1e-12 >= s->parameters.reversal_confirm_s;
            s->motion_entry_elapsed = s->motion_ramp_pending || late_entry ? 0 : s->parameters.motion_rise_s;
            s->motion_ramp_pending = false; s->direction_wait_s = 0;
        } else {
            target = 0;
        }
    }
    if (s->motion_waiting) target = 0;
    /* 一次共同限速；禁止跨过零直接换号，随后分配左右并缩放伸展侧。 */
    const double previous_torque = s->torque;
    if (target == 0 || s->motion_releasing) {
        s->motion_entry_elapsed = 0;
        if (s->torque != 0 && (s->motion_releasing || s->motion_ramp_pending)) release_motion_torque(s, dt);
        else if (s->torque != 0) apply_target(s, 0, dt);
        else { s->target = 0; s->motion_releasing = false; }
    } else {
        const double gain = smooth_fraction(s->motion_entry_elapsed/s->parameters.motion_rise_s);
        s->motion_entry_elapsed = fmin(s->parameters.motion_rise_s, s->motion_entry_elapsed+dt);
        apply_target(s, target*gain, dt);
    }
    if (previous_torque != 0 && s->torque == 0) {
        const int previous_direction = previous_torque > 0 ? 1 : -1;
        if (interrupted && (direction == 0 || direction == previous_direction)) {
            s->paused_direction = previous_direction;
            s->paused_x = s->previous_x;
        } else {
            s->completed_direction = previous_direction;
            s->paused_direction = 0;
        }
        s->torque_direction = 0;
        s->motion_waiting = true;
    }
}

static void begin_phase_pause(hip_gentle_assist_t *s, double x)
{
    s->pause_from_state = s->state;
    s->phase_paused = true; s->pause_used_this_half = true;
    s->pause_anchor_x = x; s->pause_zeroed = s->torque == 0;
    s->requested_direction = 0; s->direction_wait_s = 0;
    s->motion_waiting = true; s->torque_direction = 0; s->motion_ramp_pending = true;
}

static bool pause_can_resume(const hip_gentle_assist_t *s, double x)
{
    const hip_gentle_parameters_t *p = &s->parameters;
    if (!s->pause_zeroed || s->torque != 0) return false;
    const bool moving = s->motion_direction == -s->direction &&
        s->motion_confirm_time+1e-12 >= p->reversal_confirm_s;
    if (s->state == HIP_GENTLE_HOLD)
        return moving && s->direction*(s->pause_anchor_x-x) > fmax(p->startup_lift_angle_rad, p->hold_range_rad);
    return moving || x*s->direction <= 0;
}

static void resume_phase_pause(hip_gentle_assist_t *s, double x)
{
    const hip_gentle_parameters_t *p = &s->parameters;
    const bool held = s->state == HIP_GENTLE_HOLD;
    s->state = s->pause_from_state;
    if (s->state == HIP_GENTLE_START || s->state == HIP_GENTLE_TRACK) {
        /* shortcut: 用本半波峰值近似回落相位；若偏置/多峰导致恢复失配，再改进相位观测。 */
        const double ratio = clip(x/fmax(s->peak, p->start_angle_rad), -1, 1);
        const double phase = wrap(s->direction < 0 ? asin(ratio) : PI-asin(ratio));
        const double correction = wrap(phase-s->phases[0]+PI)-PI;
        for (unsigned k = 0; k < HIP_GENTLE_HARMONICS; ++k)
            s->phases[k] = wrap(s->phases[k]+(k+1)*correction);
        memset(s->rates, 0, sizeof(s->rates));
        s->phase = phase; s->phase_correction = 0;
        if (held) {
            s->bad_time = 0;
            if (s->state == HIP_GENTLE_START)
                s->peak_time = s->time_s-fmin(PI/2, acos(clip(ratio*s->direction, -1, 1)))/s->omega;
            else s->bootstrap_pending = true;
        }
    } else if (held && x*s->direction <= 0) {
        s->collect_wait_cycle = true; s->collect_has_zero = false;
    }
    if (held) {
        s->reason = HIP_GENTLE_RESUMED;
        s->relearn_peaks_after_hold = s->startup_released;
        s->resume_cross_seen = false;
        s->resume_sync_remaining = s->startup_released ? HIP_GENTLE_RESUME_SYNC_HALVES : 0;
        s->bootstrap_halves_remaining = 0;
        s->enabled_time = s->time_s;
        s->adaptive_high = false;
        s->adaptive_amplitude_nm = fmin(p->amplitude_nm, p->torque_limit_nm);
        s->completed_direction = s->paused_direction = 0;
        s->progress_x = x; s->progress_time = s->time_s; s->progress_returning = true;
        s->hold_quiet_time = 0;
    }
    s->phase_paused = false;
    s->requested_direction = 0; s->direction_wait_s = 0;
    s->motion_waiting = true; s->torque_direction = 0;
    s->motion_releasing = false; s->motion_entry_elapsed = 0; s->motion_ramp_pending = true;
}

static void stop(hip_gentle_assist_t *s, hip_gentle_reason_t reason, double dt)
{
    if (reason == HIP_GENTLE_WEAK_PEAK || reason == HIP_GENTLE_AMPLITUDE_DROP) {
        /* 不让刚触发退出的残余摆幅再次按低门槛启动。 */
        s->restart_angle_rad = fmax(s->restart_angle_rad,
            fmax(s->parameters.start_angle_rad, reference_peak(s, s->direction)*s->parameters.peak_drop_ratio));
    }
    s->state = HIP_GENTLE_EXIT;
    s->gate = s->phase_correction = 0;
    s->bootstrap_pending = false;
    s->startup_released = s->startup_reference_valid = false;
    s->startup_second_confirmed = false;
    s->collect_wait_cycle = s->collect_has_zero = false;
    s->hold_quiet_time = 0;
    s->pause_zeroed = false;
    s->relearn_peaks_after_hold = false;
    s->resume_cross_seen = false;
    s->resume_sync_remaining = 0; s->bootstrap_halves_remaining = 0;
    s->weak_pending = false;
    s->phase_paused = s->pause_used_this_half = false;
    s->requested_direction = 0; s->direction_wait_s = 0;
    s->startup_first_leg = 0;
    s->adaptive_high = false;
    s->adaptive_amplitude_nm = 0;
    s->torque_direction = s->motion_direction = 0;
    s->completed_direction = 0;
    s->paused_direction = 0;
    s->motion_confirm_time = 0;
    s->motion_waiting = false;
    s->motion_releasing = false;
    s->motion_entry_elapsed = 0;
    s->motion_ramp_pending = true;
    apply_target(s, 0, dt);
    s->reason = reason;
}

static double scale_joint_torque(double torque, double extension_ratio)
{
    /* 逻辑关节方向：正为屈曲，负为伸展；电机安装符号由驱动层处理。 */
    return torque < 0 ? torque*extension_ratio : torque;
}

static void output(const hip_gentle_assist_t *s, double x, hip_gentle_output_t *o)
{
    const double ratio = s->parameters.extension_torque_ratio;
    *o = (hip_gentle_output_t){
        .left_torque_nm=scale_joint_torque(s->torque, ratio),
        .right_torque_nm=scale_joint_torque(-s->torque, ratio),
        .state=s->state, .reason=s->reason, .phase_rad=s->phase,
        .frequency_hz=s->omega/TAU, .signal_rad=x, .reconstructed_rad=s->estimate,
        .swing_amplitude_rad=fabs(s->amplitudes[0]),
        .gate=s->gate, .target_left_torque_nm=scale_joint_torque(s->target, ratio),
        .target_right_torque_nm=scale_joint_torque(-s->target, ratio),
        .adaptive_amplitude_nm=s->adaptive_amplitude_nm
    };
}

static bool fail(hip_gentle_assist_t *s, hip_gentle_reason_t reason, hip_gentle_output_t *o)
{
    hip_gentle_reset(s);
    s->reason = reason;
    output(s, 0, o);
    return false;
}

static void shift_gait_times(hip_gentle_assist_t *s, double delta)
{
    s->cross_time += delta; s->peak_time += delta; s->event_time += delta; s->enabled_time += delta;
    if (s->collect_has_zero) s->collect_zero_time += delta;
}

bool hip_gentle_update(hip_gentle_assist_t *s, double left, double right,
                      uint32_t timestamp_ms, bool valid, hip_gentle_output_t *o)
{
    if (!o) return false;
    memset(o, 0, sizeof(*o));
    if (!s || !s->initialized) return false;
    const hip_gentle_parameters_t *p = &s->parameters;
    const double x = left-right;
    if (!valid || !isfinite(left) || !isfinite(right) || !isfinite(x))
        return fail(s, HIP_GENTLE_INVALID_ANGLE, o);
    double dt = s->has_timestamp ? (uint32_t)(timestamp_ms-s->timestamp_ms)/1000.0 : 0.0;
    if (s->has_timestamp && (dt <= 0 || dt > p->reset_gap_s)) {
        const double restart_angle = s->restart_angle_rad;
        hip_gentle_reset(s);
        s->restart_angle_rad = restart_angle;
        s->reason = HIP_GENTLE_TIMESTAMP_RESET;
        dt = 0;
    }
    s->timestamp_ms = timestamp_ms;
    s->has_timestamp = true;
    s->time_s += dt;
    if (s->history_count >= HIP_GENTLE_HISTORY_CAPACITY) return fail(s, HIP_GENTLE_HISTORY_OVERFLOW, o);
    const unsigned tail = (s->history_head+s->history_count)%HIP_GENTLE_HISTORY_CAPACITY;
    s->history[tail] = (hip_gentle_history_t){s->time_s, left, right};
    ++s->history_count;
    while (s->history_count > 2 &&
           s->history[(s->history_head+1U)%HIP_GENTLE_HISTORY_CAPACITY].time_s <= s->time_s-p->derivative_window_s+1e-12) {
        s->history_head = (s->history_head+1U)%HIP_GENTLE_HISTORY_CAPACITY;
        --s->history_count;
    }
    const hip_gentle_history_t *first = &s->history[s->history_head];
    const double span = s->time_s-first->time_s;
    const double vl = span > 0 ? (left-first->left_rad)/span : 0;
    const double vr = span > 0 ? (right-first->right_rad)/span : 0;
    const double speed = vl-vr;
    double old_x = s->previous_x;
    s->previous_x = x;
    if (!isfinite(vl) || !isfinite(vr) || !isfinite(speed)) return fail(s, HIP_GENTLE_INVALID_ANGLE, o);
    if (dt == 0) {
        s->progress_x = x; s->progress_time = s->time_s;
        goto done;
    }
    if (s->state == HIP_GENTLE_COLLECT || s->state == HIP_GENTLE_START || s->state == HIP_GENTLE_TRACK) {
        const double excursion = s->direction*x;
        const double peak = fmax(s->peak, excursion);
        /* 只认可本半波外摆/回落的新进展，微小往复抖动不能反复刷新看门狗。 */
        if (!s->progress_returning && peak-excursion >= p->startup_lift_angle_rad) {
            s->progress_returning = true;
            s->progress_x = x; s->progress_time = s->time_s;
        } else if (s->direction*(x-s->progress_x)*(s->progress_returning ? -1 : 1) >= p->startup_lift_angle_rad) {
            s->progress_x = x; s->progress_time = s->time_s;
        }
    }
    update_motion_permission(s, speed, dt);
    if (s->state == HIP_GENTLE_EXIT) {
        apply_target(s, 0, dt);
        if (s->torque == 0) {
            s->state = HIP_GENTLE_WAIT; s->reason = HIP_GENTLE_WAITING_NEW_CROSSING;
        }
        goto done;
    }
    const bool raised_pose = s->direction*x >= p->start_angle_rad && s->startup_first_leg != 0;
    const bool active_gait = s->state == HIP_GENTLE_COLLECT || s->state == HIP_GENTLE_START || s->state == HIP_GENTLE_TRACK;
    double suspended_dt = dt;
    if (active_gait) {
        if (!raised_pose) s->hold_quiet_time = 0;
        else {
            const double values[3] = {left, right, x};
            bool restart_window = s->hold_quiet_time == 0;
            for (unsigned i = 0; i < 3; ++i) {
                s->hold_low[i] = fmin(s->hold_low[i], values[i]);
                s->hold_high[i] = fmax(s->hold_high[i], values[i]);
                if (s->hold_high[i]-s->hold_low[i] > p->hold_range_rad) restart_window = true;
            }
            if (restart_window) {
                for (unsigned i = 0; i < 3; ++i) s->hold_low[i] = s->hold_high[i] = values[i];
                s->hold_quiet_time = 0;
            }
            s->hold_quiet_time += dt;
        }
        if (s->hold_quiet_time+1e-12 >= p->hold_confirm_s) {
            begin_phase_pause(s, x);
            s->state = HIP_GENTLE_HOLD;
            s->reason = HIP_GENTLE_HOLDING;
            s->adaptive_high = false;
            suspended_dt = s->hold_quiet_time;
            s->peak_time = fmin(s->peak_time, s->time_s-suspended_dt);
            memset(s->startup_lift_time, 0, sizeof(s->startup_lift_time));
            s->startup_min_angle[0] = left; s->startup_min_angle[1] = right;
        }
    }
    if (s->phase_paused) {
        const bool held = s->state == HIP_GENTLE_HOLD;
        if (held) {
            shift_gait_times(s, suspended_dt);
            if (s->torque != 0) release_motion_torque(s, dt);
            else s->target = 0;
            if (!s->pause_zeroed && s->torque == 0) { s->pause_zeroed = true; s->pause_anchor_x = x; }
            if (!raised_pose && s->quiet_time+1e-12 >= p->quiet_duration_s) {
                stop(s, HIP_GENTLE_QUIET, dt); goto done;
            }
            if (s->pause_zeroed && x*s->direction > s->pause_anchor_x*s->direction) {
                s->pause_anchor_x = x; s->peak = fmax(s->peak, fabs(x));
            }
        }
        if (pause_can_resume(s, x)) {
            if (held) {
                shift_gait_times(s, -s->motion_confirm_time);
                old_x = s->pause_anchor_x;
            }
            resume_phase_pause(s, x);
        } else if (held) goto done;
        /* 短停仍执行后面的无进展/静止保护；长HOLD沿原规则不限时等待。 */
    }
    const bool hold_candidate = active_gait && raised_pose && s->hold_quiet_time >= p->hold_confirm_s*.5;
    const double idle_progress = s->time_s-s->progress_time;
    const bool quiet_stop = s->quiet_time+1e-12 >= p->quiet_duration_s && idle_progress >= p->quiet_duration_s;
    if (s->state != HIP_GENTLE_WAIT && (quiet_stop || idle_progress > p->progress_timeout_s)) {
        /* 无进展时只宽限当前稳定候选；窗口被打破且仍无进展就退出。 */
        if (hold_candidate && s->quiet_time < p->quiet_duration_s) {
            s->motion_ramp_pending = true;
            apply_motion_target(s, 0, dt);
            goto done;
        }
        stop(s, s->weak_pending ? HIP_GENTLE_WEAK_PEAK :
             quiet_stop ? HIP_GENTLE_QUIET : HIP_GENTLE_NO_PROGRESS, dt);
        goto done;
    }
    if (s->state == HIP_GENTLE_WAIT) {
        if (s->quiet_time+1e-12 >= p->quiet_duration_s) s->restart_angle_rad = 0;
        const double start_angle = fmax(p->start_angle_rad, s->restart_angle_rad);
        const int direction = old_x <= start_angle && x > start_angle ? 1 :
            old_x >= -start_angle && x < -start_angle ? -1 : 0;
        if (direction) {
            s->state = HIP_GENTLE_COLLECT;
            s->direction = direction;
            s->cross_time = s->event_time = s->time_s;
            s->progress_x = x; s->progress_time = s->time_s;
            s->progress_returning = false;
            s->cross_angle_rad = start_angle;
            s->peak = fabs(x); s->peak_time = s->time_s;
            s->quiet_time = s->bad_time = 0;
            s->weak_count = 0;
            s->collect_wait_cycle = s->collect_has_zero = false;
        }
        update_startup_permission(s, left, right, vl, vr, dt);
        goto done;
    }
    update_startup_permission(s, left, right, vl, vr, dt);
    if (s->state == HIP_GENTLE_COLLECT) {
        if (x*s->direction > s->peak) { s->peak = x*s->direction; s->peak_time = s->time_s; }
        if (s->collect_wait_cycle) {
            /* 峰顶停顿使短区间定频不可信时，等两个有足够摆幅的真实过零事件。 */
            if (x*s->direction <= 0 && old_x*s->direction > 0) {
                const bool qualified = s->peak > s->cross_angle_rad;
                const double half_period = s->time_s-s->collect_zero_time;
                const double omega = half_period > 0 ? PI/half_period : 0;
                if (qualified && s->collect_has_zero &&
                    omega >= TAU*p->min_frequency_hz && omega <= TAU*p->max_frequency_hz) {
                    s->phase = s->direction > 0 ? PI : 0;
                    seed(s, s->peak, s->phase, omega);
                    s->state = HIP_GENTLE_TRACK;
                    s->collect_wait_cycle = false;
                    s->bootstrap_pending = false;
                    s->reason = HIP_GENTLE_EVENT_SEED;
                }
                s->collect_has_zero = qualified;
                if (qualified) s->collect_zero_time = s->event_time = s->time_s;
                s->direction *= -1;
                s->pause_used_this_half = false;
                s->progress_returning = false; s->progress_x = x;
                if (qualified) s->progress_time = s->time_s;
                s->peak = fabs(x); s->peak_time = s->time_s;
            }
            goto done;
        }
        if (x*s->direction <= 0) stop(s, HIP_GENTLE_INCOMPLETE_HALF_WAVE, dt);
        else if (speed*s->direction < -p->quiet_speed_rad_s) {
            const double interval = s->peak_time-s->cross_time;
            if (s->peak <= s->cross_angle_rad) stop(s, HIP_GENTLE_INVALID_FIRST_PEAK, dt);
            else {
                const double omega = interval > 0 ? (PI/2-asin(clip(s->cross_angle_rad/s->peak, 0, 1)))/interval : 0;
                if (omega < TAU*p->min_frequency_hz || omega > TAU*p->max_frequency_hz) {
                    s->collect_wait_cycle = true;
                    s->collect_has_zero = false;
                    s->reason = HIP_GENTLE_SEED_FREQUENCY_RANGE;
                }
                else {
                    s->phase = s->direction > 0 ? PI/2 : 3*PI/2;
                    seed(s, s->peak, s->phase, omega);
                    s->event_time = s->time_s;
                    s->state = HIP_GENTLE_START;
                    s->reason = HIP_GENTLE_EVENT_SEED;
                }
            }
        }
        goto done;
    }
    if (!s->phase_paused && !track(s, x, dt)) return fail(s, HIP_GENTLE_INVALID_TRACKER, o);
    if (s->state == HIP_GENTLE_START) {
        const double start_phase = s->direction > 0 ? PI/2 : 3*PI/2;
        s->phase = wrap(start_phase+fmin(PI/2, s->omega*(s->time_s-s->peak_time)));
        if (x*s->direction <= 0) {
            const double phase = s->direction > 0 ? PI : 0;
            const double correction = wrap(phase-s->phases[0]+PI)-PI;
            for (unsigned k = 0; k < HIP_GENTLE_HARMONICS; ++k)
                s->phases[k] = wrap(s->phases[k]+(k+1)*correction);
            memset(s->rates, 0, sizeof(s->rates));
            s->phase = phase; s->direction *= -1; s->peak = fabs(x);
            s->pause_used_this_half = false;
            s->progress_returning = false; s->progress_x = x; s->progress_time = s->time_s;
            s->event_time = s->time_s; s->state = HIP_GENTLE_TRACK; s->bootstrap_pending = true;
            s->bootstrap_halves_remaining = s->resume_sync_remaining ? 0 : HIP_GENTLE_STARTUP_SYNC_HALVES;
            if (s->resume_sync_remaining) {
                s->resume_cross_time = s->time_s;
                s->resume_cross_seen = true;
            }
            if (s->relearn_peaks_after_hold) {
                clear_peak_history(s);
                s->relearn_peaks_after_hold = false;
            }
        }
    } else {
        s->phase = s->phases[0];
        s->peak = fmax(s->peak, x*s->direction);
        /* 保留多半步参考，避免连续减小的摆幅把退出门槛一起带低。 */
        const double weak_limit = reference_peak(s, s->direction)*p->peak_drop_ratio;
        const bool weak_peak = s->peak < p->start_angle_rad ||
                               s->peak < weak_limit;
        if (s->weak_pending && !weak_peak) {
            s->weak_pending = false;
            if (s->startup_released) s->restart_angle_rad = 0;
        }
        if (p->weak_peak_stop_enabled && weak_peak && speed*s->direction < -p->quiet_speed_rad_s) {
            s->weak_pending = true;
            /* 候选期的周期异常也不能让残余摆动重新按低门槛启动。 */
            s->restart_angle_rad = fmax(s->restart_angle_rad, fmax(p->start_angle_rad, weak_limit));
        }
        if (x*s->direction <= 0 && old_x*s->direction > 0) {
            if (s->weak_pending) { stop(s, HIP_GENTLE_WEAK_PEAK, dt); goto done; }
            const bool weak = s->peak < p->start_angle_rad || s->peak < weak_limit;
            s->weak_count = weak ? s->weak_count+1U : 0U;
            if (s->weak_count >= p->weak_half_cycles) { stop(s, HIP_GENTLE_AMPLITUDE_DROP, dt); goto done; }
            if (!weak) {
                const double event_phase = s->direction > 0 ? PI : 0;
                double error = wrap(event_phase-s->phase+PI)-PI;
                const bool resynchronizing = s->resume_sync_remaining != 0;
                if (resynchronizing) {
                    /* 不把悬停前的节奏当作恢复后的节奏；仅有限次使用有效整半周期重同步。 */
                    if (s->resume_cross_seen) {
                        const double period = s->time_s-s->resume_cross_time;
                        const double omega = period > 0 ? PI/period : 0;
                        if (omega < TAU*p->min_frequency_hz || omega > TAU*p->max_frequency_hz) {
                            stop(s, HIP_GENTLE_SEED_FREQUENCY_RANGE, dt); goto done;
                        }
                        s->omega = omega;
                        --s->resume_sync_remaining;
                    }
                    s->rates[1] = 0;
                    for (unsigned k = 0; k < HIP_GENTLE_HARMONICS; ++k)
                        s->rates[7+k] = (k+1)*s->omega;
                    s->resume_cross_time = s->time_s;
                    s->resume_cross_seen = true;
                }
                /* 首次完整左右周期内允许纠正定频，正常跟踪不反复重置掩盖失锁。 */
                if (!resynchronizing && s->bootstrap_pending && fabs(error) > p->phase_error_limit_rad) {
                    const double half_period = s->time_s-s->event_time;
                    const double omega = half_period > 0 ? PI/half_period : 0;
                    if (omega >= TAU*p->min_frequency_hz && omega <= TAU*p->max_frequency_hz) {
                        s->omega = omega;
                        for (unsigned k = 0; k < HIP_GENTLE_HARMONICS; ++k)
                            s->phases[k] = wrap(s->phases[k]+(k+1)*error);
                        memset(s->rates, 0, sizeof(s->rates));
                        s->phase = event_phase; s->enabled_time = s->time_s;
                        s->reason = HIP_GENTLE_STARTUP_RESEED; error = 0;
                    }
                }
                if (!resynchronizing && s->bootstrap_halves_remaining) --s->bootstrap_halves_remaining;
                s->bootstrap_pending = s->resume_sync_remaining != 0 || s->bootstrap_halves_remaining != 0;
                if (!resynchronizing && fabs(error) > p->phase_error_limit_rad) { stop(s, HIP_GENTLE_EVENT_MISMATCH, dt); goto done; }
                s->phase_correction = error;
                /* 悬停落腿本身仍受原弱峰保护，完成后不把高抬腿当作新步幅基准。 */
                if (s->relearn_peaks_after_hold) {
                    clear_peak_history(s);
                    s->relearn_peaks_after_hold = false;
                } else remember_peak(s, s->peak, s->direction);
                s->event_time = s->time_s;
            }
            s->peak = fabs(x); s->direction *= -1;
            s->pause_used_this_half = false;
            s->progress_returning = false; s->progress_x = x; s->progress_time = s->time_s;
        }
        s->bad_time = fabs(x-s->estimate) > p->tracking_error_rad ? s->bad_time+dt : 0;
        if (p->tracking_error_stop_enabled && s->bad_time >= p->tracking_error_duration_s) {
            stop(s, HIP_GENTLE_TRACKING_ERROR, dt); goto done;
        }
    }
    /* 首腿回落仅用于学习相位，第二腿确认前不生成力矩，也不消耗渐入时间。 */
    if (!s->startup_released) {
        s->gate = 0;
        apply_target(s, 0, dt);
        goto done;
    }
    const double elapsed = s->time_s-s->enabled_time;
    s->gate = fmin(1, elapsed/p->inner_rise_s)*fmin(1, elapsed/p->outer_rise_s);
    if (!update_adaptive_amplitude(s, dt)) return fail(s, HIP_GENTLE_INVALID_TRACKER, o);
    double target = s->gate*s->adaptive_amplitude_nm*hip_gentle_waveform(s->phase+p->response/100-0.92, p->width, p->asymmetry);
    if (!isfinite(target)) return fail(s, HIP_GENTLE_INVALID_TRACKER, o);
    if (s->state == HIP_GENTLE_START && target*s->direction > 0) target = 0;
    const int request = target > 0 ? 1 : target < 0 ? -1 : 0;
    const bool direction_ready = s->motion_direction == request &&
        s->motion_confirm_time+1e-12 >= p->reversal_confirm_s;
    /* 交接等待也可能仍有前一方向的运动，不能只靠连续低速识别。 */
    const bool direction_waited = s->requested_direction == request &&
        s->direction_wait_s+1e-12 >= p->low_speed_release_confirm_s;
    const bool relative_pause = s->quiet_time+1e-12 >= p->reversal_confirm_s;
    if (s->state == HIP_GENTLE_TRACK && s->startup_released && !s->bootstrap_pending && s->gate >= 1 &&
        !s->weak_pending && !hold_candidate && !s->pause_used_this_half && raised_pose && s->torque == 0 &&
        (relative_pause || direction_waited) && request == -s->direction &&
        request != s->completed_direction && request != s->torque_direction && !direction_ready) {
        begin_phase_pause(s, x);
    }
    if (s->weak_pending || s->phase_paused) {
        s->motion_waiting = true; s->torque_direction = 0;
        s->motion_ramp_pending = true; target = 0;
    }
    if (hold_candidate) { s->motion_ramp_pending = true; target = 0; }
    apply_motion_target(s, target, dt);
done:
    if (!isfinite(s->torque) || !isfinite(s->target) || fabs(s->torque) > p->torque_limit_nm+1e-9)
        return fail(s, HIP_GENTLE_INVALID_TRACKER, o);
    output(s, x, o);
    return true;
}

/* 主机可执行检查，不进入 ESP32 构建：
 * cc -std=c11 -DHIP_GENTLE_ASSIST_SELF_TEST hip_gentle_assist.c -lm -o D:/AgentWork/hip_gentle_check.exe
 */
#ifdef HIP_GENTLE_ASSIST_SELF_TEST
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdio.h>

static void test_alternating_startup(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.amplitude_nm = 8; p.extension_torque_ratio = 2.0/3.0;
    p.start_angle_rad = .35; p.quiet_duration_s = 2;
    p.hold_confirm_s = ASSIST_GENTLE_HOLD_CONFIRM_S;
    p.inner_rise_s = .6; p.outer_rise_s = .45;
    p.weak_peak_stop_enabled = true;
    hip_gentle_assist_t s;
    hip_gentle_output_t o;
    for (unsigned first = 0; first < 2; ++first) {
        /* 单腿反复抬放、后摆跨过角度差零点，都不能代替另一腿抬起。 */
        assert(hip_gentle_init(&s, &p));
        for (unsigned i = 0; i <= 300; ++i) {
            double q[2] = {0, 0};
            q[first] = .7*sin(TAU*i*.01);
            assert(hip_gentle_update(&s, q[0], q[1], i*10, true, &o));
            assert(!s.startup_released && o.left_torque_nm == 0 && o.right_torque_nm == 0);
        }
        /* 第一腿完整抬放后，第二腿才从静止抬起；左右起步对称。 */
        assert(hip_gentle_init(&s, &p));
        bool released = false, assisted = false;
        double held[2] = {0, 0};
        for (unsigned i = 0; i <= 110; ++i) {
            const double t = i*.01;
            double q[2] = {0, 0};
            if (i < 80) q[first] = .65*sin(PI*t/.8);
            else q[1U-first] = .65*sin(PI*(t-.8)/.8);
            assert(hip_gentle_update(&s, q[0], q[1], i*10, true, &o));
            if (i < 85) assert(!s.startup_released && o.left_torque_nm == 0 && o.right_torque_nm == 0);
            if (s.startup_released && !released) {
                assert(s.startup_first_leg == first+1U);
                assert(o.gate == 0 && o.left_torque_nm == 0 && o.right_torque_nm == 0);
                released = true;
            }
            if (fabs(o.left_torque_nm)+fabs(o.right_torque_nm) > .01) {
                const double first_torque = first == 0 ? o.left_torque_nm : o.right_torque_nm;
                const double second_torque = first == 0 ? o.right_torque_nm : o.left_torque_nm;
                assert(second_torque > 0 && first_torque < 0);
                assert(fabs(first_torque+second_torque*2.0/3.0) < 1e-9);
                assisted = true;
            }
            assert(fabs(o.left_torque_nm) <= 8 && fabs(o.right_torque_nm) <= 8);
            held[0] = q[0]; held[1] = q[1];
        }
        assert(released && assisted);
        /* 短停保留许可；稳定抬腿长停进入HOLD并保持零输出。 */
        for (unsigned i = 111; i <= 120; ++i) {
            assert(hip_gentle_update(&s, held[0], held[1], i*10, true, &o));
            assert(s.startup_released);
        }
        for (unsigned i = 121; i <= 400; ++i)
            assert(hip_gentle_update(&s, held[0], held[1], i*10, true, &o));
        assert(s.state == HIP_GENTLE_HOLD && s.startup_released && o.left_torque_nm == 0 && o.right_torque_nm == 0);
        /* 另一腿的小幅抖动即使速度较大，也不足以放行。 */
        assert(hip_gentle_init(&s, &p));
        for (unsigned i = 0; i <= 140; ++i) {
            double q[2] = {0, 0};
            if (i < 80) q[first] = .65*sin(PI*i*.01/.8);
            else q[1U-first] = .01*sin(TAU*8*(i-80)*.01);
            assert(hip_gentle_update(&s, q[0], q[1], i*10, true, &o));
            assert(!s.startup_released && o.left_torque_nm == 0 && o.right_torque_nm == 0);
        }
        /* 首次动作无法分清先后时，后续真正交替仍能完成确认，不必重启。 */
        assert(hip_gentle_init(&s, &p));
        released = false;
        for (unsigned i = 0; i <= 205; ++i) {
            double q[2] = {0, 0};
            if (i < 80) {
                q[first] = .9*sin(PI*i*.01/.8);
                q[1U-first] = .45*sin(PI*i*.01/.8);
            } else if (i < 160) {
                q[1U-first] = .65*sin(PI*(i-80)*.01/.8);
            } else {
                q[first] = .65*sin(PI*(i-160)*.01/.8);
            }
            assert(hip_gentle_update(&s, q[0], q[1], i*10, true, &o));
            if (i < 165) assert(!s.startup_released && o.left_torque_nm == 0 && o.right_torque_nm == 0);
            if (i == 110) assert(s.startup_first_leg == 2U-first);
            released |= s.startup_released;
        }
        assert(released);
    }
    /* 同时屈曲不应被记录为先后两次抬腿。 */
    assert(hip_gentle_init(&s, &p));
    for (unsigned i = 0; i <= 160; ++i) {
        const double q = sin(PI*i*.01/.8);
        assert(hip_gentle_update(&s, .9*q, .45*q, i*10, true, &o));
        assert(!s.startup_released && o.left_torque_nm == 0 && o.right_torque_nm == 0);
    }
    hip_gentle_parameters_t invalid = p;
    invalid.startup_lift_angle_rad = 0;
    assert(!hip_gentle_init(&s, &invalid));
    invalid = p; invalid.startup_lift_speed_rad_s = NAN;
    assert(!hip_gentle_init(&s, &invalid));
    invalid = p; invalid.startup_lift_confirm_s = p.progress_timeout_s;
    assert(!hip_gentle_init(&s, &invalid));
}

static void test_small_slow_steps(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.amplitude_nm = ASSIST_GENTLE_AMPLITUDE_NM;
    p.torque_limit_nm = ASSIST_GENTLE_TORQUE_LIMIT_NM;
    p.extension_torque_ratio = ASSIST_GENTLE_EXTENSION_TORQUE_RATIO;
    p.adaptive_enabled = ASSIST_GENTLE_ADAPTIVE_ENABLED;
    p.adaptive_max_nm = ASSIST_GENTLE_ADAPTIVE_MAX_NM;
    p.swing_up_threshold_rad = ASSIST_GENTLE_SWING_UP_THRESHOLD_RAD;
    p.swing_down_threshold_rad = ASSIST_GENTLE_SWING_DOWN_THRESHOLD_RAD;
    p.amplitude_increase_nm_s = ASSIST_GENTLE_AMPLITUDE_INCREASE_NM_S;
    p.amplitude_decrease_nm_s = ASSIST_GENTLE_AMPLITUDE_DECREASE_NM_S;
    p.start_angle_rad = ASSIST_GENTLE_START_ANGLE_RAD;
    p.startup_lift_angle_rad = ASSIST_GENTLE_STARTUP_LIFT_ANGLE_RAD;
    p.startup_lift_speed_rad_s = ASSIST_GENTLE_STARTUP_LIFT_SPEED_RAD_S;
    p.startup_lift_confirm_s = ASSIST_GENTLE_STARTUP_LIFT_CONFIRM_S;
    p.quiet_speed_rad_s = ASSIST_GENTLE_QUIET_SPEED_RAD_S;
    p.reversal_speed_rad_s = ASSIST_GENTLE_REVERSAL_SPEED_RAD_S;
    p.quiet_duration_s = ASSIST_GENTLE_QUIET_DURATION_S;
    p.hold_confirm_s = ASSIST_GENTLE_HOLD_CONFIRM_S;
    p.low_speed_release_confirm_s = ASSIST_GENTLE_LOW_SPEED_RELEASE_CONFIRM_S;
    p.reversal_confirm_s = ASSIST_GENTLE_REVERSAL_CONFIRM_S;
    p.slew_nm_s = ASSIST_GENTLE_SLEW_NM_S;
    p.inner_rise_s = ASSIST_GENTLE_INNER_RISE_S;
    p.outer_rise_s = ASSIST_GENTLE_OUTER_RISE_S;
    p.weak_peak_stop_enabled = true;
    p.peak_drop_ratio = ASSIST_GENTLE_PEAK_DROP_RATIO;
    for (unsigned period = 200; period <= 250; period += 50) {
        for (int sign = -1; sign <= 1; sign += 2) {
            hip_gentle_assist_t s;
            hip_gentle_output_t o;
            assert(hip_gentle_init(&s, &p));
            unsigned active_cycles = 0;
            double previous = 0, held = 0;
            const unsigned end = 6*period+period/8;
            /* 角度差峰值0.24rad低于旧起步门槛，完整周期2/2.5秒。 */
            for (unsigned i = 0; i <= end; ++i) {
                held = sign*.12*sin(TAU*i/period);
                assert(hip_gentle_update(&s, held, -held, i*10, true, &o));
                assert(fabs(o.left_torque_nm) <= 8 && o.right_torque_nm == -o.left_torque_nm);
                assert(fabs(o.left_torque_nm-previous) <= p.slew_nm_s*.01+1e-9);
                if (fabs(o.left_torque_nm) > .01) active_cycles |= 1U << (i/period);
                if (i > period) assert(s.state == HIP_GENTLE_TRACK && s.startup_released);
                previous = o.left_torque_nm;
            }
            assert((active_cycles & 0x3EU) == 0x3EU);
            /* 真正停住后仍卸力、退出，不能为了慢走连续性保留静止出力。 */
            for (unsigned j = 1; j <= 260; ++j) {
                assert(hip_gentle_update(&s, held, -held, (end+j)*10, true, &o));
                if (j >= 30) assert(o.left_torque_nm == 0 && o.right_torque_nm == 0);
            }
            assert(s.state == HIP_GENTLE_WAIT && !s.startup_released);
            hip_gentle_reset(&s);
            for (unsigned i = 0; i < 200; ++i) {
                const double jitter = .01*sin(TAU*8*i*.01);
                assert(hip_gentle_update(&s, jitter, -jitter, i*10, true, &o));
                assert(o.left_torque_nm == 0 && o.right_torque_nm == 0);
            }
        }
    }
}

static void test_weak_peak_rearm(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.amplitude_nm = 4; p.torque_limit_nm = 8; p.extension_torque_ratio = 1;
    p.start_angle_rad = .2; p.quiet_speed_rad_s = .1; p.quiet_duration_s = 2;
    p.hold_confirm_s = ASSIST_GENTLE_HOLD_CONFIRM_S;
    p.reversal_speed_rad_s = .2; p.peak_drop_ratio = .65;
    p.weak_peak_stop_enabled = true; p.slew_nm_s = 50;
    hip_gentle_assist_t s;
    hip_gentle_output_t o;
    assert(hip_gentle_init(&s, &p));
    seed(&s, .8, PI/2, TAU);
    remember_peak(&s, .76, -1); remember_peak(&s, .72, 1); remember_peak(&s, .68, -1);
    assert(.49 > s.last_peak*p.peak_drop_ratio && .49 < reference_peak(&s, 1)*p.peak_drop_ratio);
    s.state = HIP_GENTLE_TRACK; s.direction = 1; s.peak = .49;
    s.torque = 8; s.startup_released = true;
    assert(hip_gentle_update(&s, .245, -.245, 0, true, &o));
    assert(hip_gentle_update(&s, .235, -.235, 10, true, &o));
    assert(o.state == HIP_GENTLE_TRACK && s.weak_pending && s.startup_released);
    assert(o.target_left_torque_nm == 0 && o.left_torque_nm < 8 && o.left_torque_nm > 7.5);
    assert(fabs(s.restart_angle_rad-.52) < 1e-12);
    const hip_gentle_assist_t pending = s;
    assert(hip_gentle_update(&s, .2, -.2, 40, true, &o));
    assert(o.reason == HIP_GENTLE_TIMESTAMP_RESET && o.left_torque_nm == 0);
    assert(!s.weak_pending && fabs(s.restart_angle_rad-.52) < 1e-12);
    s = pending;
    uint32_t now = 10;
    double previous = s.torque;
    for (unsigned i = 1; i <= 25; ++i) {
        now += 10;
        const double q = .235*(1-i/25.0);
        assert(hip_gentle_update(&s, q, -q, now, true, &o));
        assert(o.left_torque_nm >= 0 && o.left_torque_nm <= previous);
        assert(previous-o.left_torque_nm <= .5+1e-12);
        if (i < 25) assert(s.state == HIP_GENTLE_TRACK && s.weak_pending && s.startup_released);
        previous = o.left_torque_nm;
    }
    assert(s.state == HIP_GENTLE_EXIT && s.reason == HIP_GENTLE_WEAK_PEAK && !s.weak_pending);
    /* 比正常起步门槛大的残余摆动不能重启；周期异常不能清掉重启门槛。 */
    for (unsigned i = 1; i <= 400; ++i) {
        now += i == 100 ? 40U : 10U;
        const double q = .225*sin(TAU*2*i*.01);
        assert(hip_gentle_update(&s, q, -q, now, true, &o));
        assert(fabs(o.left_torque_nm) <= fabs(previous));
        assert(fabs(s.restart_angle_rad-.52) < 1e-12 && !s.startup_released);
        if (i >= 20) assert(o.left_torque_nm == 0 && o.state == HIP_GENTLE_WAIT);
        previous = o.left_torque_nm;
    }
    /* 站稳后恢复小步起步能力，不能永久锁住正常舒缓助力。 */
    for (unsigned i = 0; i < 240; ++i) {
        now += 10;
        assert(hip_gentle_update(&s, 0, 0, now, true, &o));
        assert(o.left_torque_nm == 0);
    }
    assert(s.restart_angle_rad == 0 && s.state == HIP_GENTLE_WAIT);
    bool assisted = false;
    for (unsigned i = 0; i < 500; ++i) {
        now += 10;
        const double q = .12*sin(TAU*.5*i*.01);
        assert(hip_gentle_update(&s, q, -q, now, true, &o));
        assisted |= fabs(o.left_torque_nm) > .01;
    }
    assert(assisted);
    /* 更明确的新步幅也能重新起步，定频必须使用实际提高后的越阈值角度。 */
    assert(hip_gentle_init(&s, &p));
    s.restart_angle_rad = .52;
    bool seeded = false;
    assisted = false;
    for (unsigned i = 0; i < 200; ++i) {
        const double q = .4*sin(TAU*i*.01);
        assert(hip_gentle_update(&s, q, -q, i*10, true, &o));
        if (s.state == HIP_GENTLE_START && !seeded) {
            assert(s.cross_angle_rad == .52 && o.frequency_hz > .8 && o.frequency_hz < 1.2);
            seeded = true;
        }
        assisted |= fabs(o.left_torque_nm) > .01;
    }
    assert(seeded && assisted && s.restart_angle_rad == 0);
    hip_gentle_reset(&s);
    assert(s.restart_angle_rad == 0 && s.recent_peak_count[0] == 0 && s.recent_peak_count[1] == 0);
}

static void test_double_peak(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.amplitude_nm = 4; p.extension_torque_ratio = 1; p.start_angle_rad = .2;
    p.peak_drop_ratio = .65; p.weak_peak_stop_enabled = true;
    p.reversal_speed_rad_s = .2; p.slew_nm_s = 50;
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_assist_t s;
        hip_gentle_output_t o;
        assert(hip_gentle_init(&s, &p));
        seed(&s, .5, sign>0 ? PI/2 : 3*PI/2, TAU);
        remember_peak(&s, .5, sign);
        s.state = HIP_GENTLE_TRACK; s.direction = sign; s.peak = .25;
        s.startup_released = true; s.startup_first_leg = sign>0 ? 1 : 2;
        s.enabled_time = -10; s.torque = sign*2;
        assert(hip_gentle_update(&s, sign*.125, -sign*.125, 0, true, &o));
        bool pending_seen = false, recovered = false;
        double previous = s.torque;
        for (unsigned i = 1; i <= 40; ++i) {
            const double x = sign*(i<=10 ? .25-.01*i : .15+.01*(i-10));
            assert(hip_gentle_update(&s, x/2, -x/2, i*10, true, &o));
            assert(s.state == HIP_GENTLE_TRACK && s.startup_released && s.enabled_time == -10);
            assert(s.gate == 1 && fabs(s.torque-previous) <= .5+1e-12);
            if (s.weak_pending) {
                pending_seen = true;
                assert(s.target == 0 && fabs(s.torque) <= fabs(previous));
                assert(s.restart_angle_rad >= .325-1e-12);
            } else if (pending_seen) {
                recovered = true;
                assert(s.peak >= .325-1e-12 && s.restart_angle_rad == 0);
            }
            previous = s.torque;
        }
        assert(pending_seen && recovered);
        /* 未过零也未继续增大的普通停步，仍按无进展保护退出。 */
        assert(hip_gentle_init(&s, &p));
        seed(&s, .5, sign>0 ? PI/2 : 3*PI/2, TAU);
        remember_peak(&s, .5, sign);
        s.state = HIP_GENTLE_TRACK; s.direction = sign; s.peak = .25;
        assert(hip_gentle_update(&s, sign*.125, -sign*.125, 0, true, &o));
        bool stopped = false;
        for (unsigned i = 1; i <= 240; ++i) {
            assert(hip_gentle_update(&s, sign*.12, -sign*.12, i*10, true, &o));
            if (s.state == HIP_GENTLE_EXIT) {
                assert(s.reason == HIP_GENTLE_WEAK_PEAK);
                stopped = true;
            }
        }
        assert(stopped && s.state == HIP_GENTLE_WAIT && s.torque == 0);
        s.weak_pending = true;
        assert(!hip_gentle_update(&s, NAN, 0, 410, true, &o));
        assert(s.state == HIP_GENTLE_WAIT && !s.weak_pending && o.left_torque_nm == 0);
    }
}

static void test_pacing_recovery(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.start_angle_rad = .2; p.peak_drop_ratio = .65; p.quiet_duration_s = 2;
    p.hold_confirm_s = ASSIST_GENTLE_HOLD_CONFIRM_S;
    p.weak_peak_stop_enabled = true;
    hip_gentle_assist_t s;
    hip_gentle_output_t o;
    assert(hip_gentle_init(&s, &p));
    seed(&s, .374, PI/2, TAU);
    remember_peak(&s, .489, -1); remember_peak(&s, .358, 1); remember_peak(&s, .359, -1);
    assert(reference_peak(&s, 1) == .374 && reference_peak(&s, -1) == .489);
    s.state = HIP_GENTLE_TRACK; s.direction = 1; s.peak = .315;
    assert(hip_gentle_update(&s, .1575, -.1575, 0, true, &o));
    assert(hip_gentle_update(&s, .155, -.155, 10, true, &o));
    assert(s.state == HIP_GENTLE_TRACK); // 另一侧较大峰值不能让本侧正常小步退出。

    for (int sign = -1; sign <= 1; sign += 2) {
        assert(hip_gentle_init(&s, &p));
        s.restart_angle_rad = .32;
        bool waited = false, seeded = false;
        for (unsigned i = 0; i <= 210; ++i) {
            const double t = i*.01;
            double x;
            if (t <= .30) x = .321*t/.30;
            else if (t <= .60) x = .321+.039*(t-.30)/.30;
            else if (t <= .90) x = .36-1.2*(t-.60);
            else if (t <= 1.40) x = -.72*(t-.90);
            else x = -.36+.72*(t-1.40);
            assert(hip_gentle_update(&s, sign*x/2, -sign*x/2, i*10, true, &o));
            assert(s.state != HIP_GENTLE_EXIT);
            if (s.collect_wait_cycle) {
                waited = true;
                assert(s.state == HIP_GENTLE_COLLECT && o.left_torque_nm == 0);
            }
            if (s.state == HIP_GENTLE_TRACK) {
                seeded = true;
                assert(o.frequency_hz > .45 && o.frequency_hz < .55);
            }
        }
        assert(waited && seeded && s.startup_released);
        hip_gentle_reset(&s);
        assert(!s.collect_wait_cycle && !s.collect_has_zero && !s.startup_second_confirmed);
    }
}

static void test_swing_tiers(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.amplitude_nm = ASSIST_GENTLE_AMPLITUDE_NM;
    p.torque_limit_nm = ASSIST_GENTLE_TORQUE_LIMIT_NM;
    p.adaptive_enabled = true; p.adaptive_max_nm = ASSIST_GENTLE_ADAPTIVE_MAX_NM;
    p.swing_up_threshold_rad = ASSIST_GENTLE_SWING_UP_THRESHOLD_RAD;
    p.swing_down_threshold_rad = ASSIST_GENTLE_SWING_DOWN_THRESHOLD_RAD;
    p.amplitude_increase_nm_s = ASSIST_GENTLE_AMPLITUDE_INCREASE_NM_S;
    p.amplitude_decrease_nm_s = ASSIST_GENTLE_AMPLITUDE_DECREASE_NM_S;
    const double middle = (p.swing_up_threshold_rad+p.swing_down_threshold_rad)/2;
    hip_gentle_assist_t s;
    for (int sign = -1; sign <= 1; sign += 2) {
        assert(hip_gentle_init(&s, &p));
        seed(&s, .4, 0, TAU*2); s.state = HIP_GENTLE_TRACK; s.gate = 1;
        s.previous_x = s.peak = 2; // 瞬时角度差或历史峰值不替代AFO幅值。
        s.amplitudes[0] = sign*.4;
        for (unsigned i = 0; i < 200; ++i) assert(update_adaptive_amplitude(&s,.01));
        assert(!s.adaptive_high && s.adaptive_amplitude_nm == 4); // 小摆幅快节奏仍低档。
        s.amplitudes[0] = sign*middle;
        assert(update_adaptive_amplitude(&s,.01) && !s.adaptive_high);
        s.omega = TAU*.3; s.amplitudes[0] = sign*p.swing_up_threshold_rad;
        assert(update_adaptive_amplitude(&s,.01) && s.adaptive_high);
        assert(fabs(s.adaptive_amplitude_nm-4.03)<1e-9); // 大摆幅慢节奏升档，幅值不跳变。
        for (unsigned i = 0; i < 200; ++i) {
            s.amplitudes[0] = sign*(middle+.04*sin(TAU*i*.1));
            assert(update_adaptive_amplitude(&s,.01) && s.adaptive_high);
        }
        assert(s.adaptive_amplitude_nm == 8);
        s.amplitudes[0] = sign*p.swing_down_threshold_rad;
        assert(update_adaptive_amplitude(&s,.01) && !s.adaptive_high);
        assert(fabs(s.adaptive_amplitude_nm-7.94)<1e-9);
        for (unsigned i = 0; i < 200; ++i) {
            s.amplitudes[0] = sign*(middle+.04*sin(TAU*i*.1));
            assert(update_adaptive_amplitude(&s,.01) && !s.adaptive_high);
        }
        assert(s.adaptive_amplitude_nm == 4);
    }
    for (unsigned guard = 0; guard < 4; ++guard) {
        assert(hip_gentle_init(&s, &p));
        seed(&s, 1, 0, TAU); s.state = HIP_GENTLE_TRACK; s.gate = 1;
        s.adaptive_high = true; s.adaptive_amplitude_nm = 8;
        if (guard == 0) s.state = HIP_GENTLE_START;
        if (guard == 1) s.bootstrap_pending = true;
        if (guard == 2) s.gate = .5;
        if (guard == 3) s.weak_pending = true;
        assert(update_adaptive_amplitude(&s,.01) && !s.adaptive_high);
        assert(fabs(s.adaptive_amplitude_nm-7.94)<1e-9);
        s.state = HIP_GENTLE_TRACK; s.bootstrap_pending = s.weak_pending = false; s.gate = 1;
        s.amplitudes[0] = middle;
        assert(update_adaptive_amplitude(&s,.01) && !s.adaptive_high);
    }
    s.adaptive_high = true;
    seed(&s, 1, 0, TAU);
    assert(!s.adaptive_high && s.adaptive_amplitude_nm == 4);
    s.adaptive_high = true;
    stop(&s, HIP_GENTLE_QUIET, .01);
    assert(!s.adaptive_high && s.target == 0 && s.adaptive_amplitude_nm == 0);
    s.adaptive_high = true; s.parameters.adaptive_enabled = false;
    assert(update_adaptive_amplitude(&s,.01) && !s.adaptive_high && s.adaptive_amplitude_nm == 4);
    s.parameters.adaptive_enabled = true; s.amplitudes[0] = NAN;
    assert(!update_adaptive_amplitude(&s,.01));
    hip_gentle_reset(&s);
    assert(!s.adaptive_high);
    hip_gentle_parameters_t invalid = p;
    invalid.swing_down_threshold_rad = p.swing_up_threshold_rad;
    assert(!hip_gentle_init(&s, &invalid));
    invalid.swing_down_threshold_rad = p.swing_up_threshold_rad+.1;
    assert(!hip_gentle_init(&s, &invalid));
    invalid = p; invalid.swing_down_threshold_rad = 0;
    assert(!hip_gentle_init(&s, &invalid));
    invalid = p; invalid.swing_up_threshold_rad = NAN;
    assert(!hip_gentle_init(&s, &invalid));
}

static void test_torque_handoff(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.slew_nm_s = 50; p.motion_rise_s = .16; p.motion_fall_s = .08;
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_assist_t s;
        assert(hip_gentle_init(&s, &p));
        s.motion_ramp_pending = false;
        for (unsigned i = 0; i < 3; ++i) {
            update_motion_permission(&s, sign, .01);
            apply_motion_target(&s, sign*.2, .01);
        }
        assert(s.torque == sign*.2 && s.motion_entry_elapsed == p.motion_rise_s);
        hip_gentle_reset(&s); s.motion_ramp_pending = false;
        for (unsigned i = 0; i < 3; ++i) {
            update_motion_permission(&s, sign, .01);
            apply_motion_target(&s, sign*3, .01);
        }
        assert(s.torque == 0);
        apply_motion_target(&s, sign*3, .01);
        assert(s.torque*sign > 0 && fabs(s.torque)<.05);
        for (unsigned i = 0; i < 30; ++i) apply_motion_target(&s, sign*3, .01);
        assert(s.torque == sign*3);
        /* 许可获准却还未出力就换相，不能锁住下一方向。 */
        hip_gentle_reset(&s); s.motion_ramp_pending = false; s.completed_direction = sign;
        for (unsigned i = 0; i < 3; ++i) {
            update_motion_permission(&s, -sign, .01);
            apply_motion_target(&s, -sign*3, .01);
        }
        assert(s.torque == 0 && s.torque_direction == -sign);
        for (unsigned i = 0; i < 3; ++i) {
            update_motion_permission(&s, sign, .01);
            apply_motion_target(&s, sign*3, .01);
        }
        assert(s.torque == 0 && s.completed_direction == -sign && s.torque_direction == sign);
        apply_motion_target(&s, sign*3, .01);
        assert(s.torque*sign > 0 && fabs(s.torque)<.05);
        const double torques[] = {.05, .5, 2, 8}, rates[] = {-50, -10, 0, 10, 50};
        for (unsigned q = 0; q < 4; ++q) for (unsigned v = 0; v < 5; ++v) {
            hip_gentle_reset(&s);
            s.torque = sign*torques[q]; s.torque_rate_nm_s = sign*rates[v];
            bool zeroed = false;
            double previous = s.torque;
            for (unsigned i = 0; i < 250; ++i) {
                release_motion_torque(&s, .001);
                assert(s.torque*sign >= 0 && fabs(s.torque)<=fabs(previous));
                assert(fabs(s.torque-previous)<=.05+1e-12 && s.target == 0);
                if (i == 0 && rates[v]<0 && -rates[v]<=3*torques[q]/s.motion_release_duration)
                    assert(fabs(s.torque_rate_nm_s-sign*rates[v])<2);
                previous = s.torque;
                if (s.torque == 0) { zeroed = true; break; }
            }
            assert(zeroed && !s.motion_releasing);
        }
    }
}

static void test_phase_paused(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.start_angle_rad = .2; p.reversal_speed_rad_s = .2; p.reversal_confirm_s = .03;
    p.quiet_duration_s = 2; p.hold_confirm_s = 1; p.hold_range_rad = .08;
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_assist_t s;
        hip_gentle_output_t o;
        assert(hip_gentle_init(&s, &p));
        seed(&s, .6, (sign>0 ? PI/2 : 3*PI/2)+.35, TAU);
        s.state = HIP_GENTLE_TRACK; s.direction = sign; s.peak = .6;
        s.startup_first_leg = sign>0 ? 1 : 2; s.startup_released = true;
        s.enabled_time = -10; s.gate = 1; s.completed_direction = sign;
        uint32_t now = 0;
        assert(hip_gentle_update(&s, sign*.3, -sign*.3, now, true, &o));
        for (unsigned i = 0; i < 3; ++i) {
            now += 10;
            assert(hip_gentle_update(&s, sign*.3, -sign*.3, now, true, &o));
        }
        assert(s.phase_paused && s.pause_used_this_half && s.torque == 0);
        const double paused_phase = s.phase, paused_time = s.time_s, event_time = s.event_time;
        for (unsigned i = 0; i < 20; ++i) {
            now += 10;
            assert(hip_gentle_update(&s, sign*.3, -sign*.3, now, true, &o));
            assert(s.phase_paused && s.phase == paused_phase && s.torque == 0);
        }
        assert(s.time_s > paused_time+.19 && s.event_time == event_time);
        const hip_gentle_assist_t held = s;
        const uint32_t held_ms = now;
        // 短停升级为长HOLD仍保持同一暂停状态，不重新出力。
        for (unsigned i = 0; i < 120; ++i) {
            now += 10;
            assert(hip_gentle_update(&s, sign*.3, -sign*.3, now, true, &o));
            assert(s.phase_paused && s.phase == paused_phase && s.torque == 0);
        }
        assert(s.state == HIP_GENTLE_HOLD && s.pause_from_state == HIP_GENTLE_TRACK);
        bool long_resumed = false;
        for (unsigned i = 1; i <= 20; ++i) {
            now += 10;
            const double x = sign*(.6-.015*i);
            assert(hip_gentle_update(&s, x/2, -x/2, now, true, &o));
            if (s.state != HIP_GENTLE_HOLD && !long_resumed) {
                assert(!s.phase_paused && s.state == HIP_GENTLE_TRACK && s.startup_released);
                assert(s.torque == 0 && s.completed_direction == 0 && s.gate == 0);
                assert(s.bootstrap_pending && s.resume_sync_remaining == HIP_GENTLE_RESUME_SYNC_HALVES);
                assert(!s.adaptive_high && s.adaptive_amplitude_nm == p.amplitude_nm);
                long_resumed = true;
            }
        }
        assert(long_resumed);
        s = held; now = held_ms;
        bool resumed = false;
        for (unsigned i = 1; i <= 12; ++i) {
            now += 10;
            const double x = sign*(.6-.01*i);
            assert(hip_gentle_update(&s, x/2, -x/2, now, true, &o));
            assert(s.state == HIP_GENTLE_TRACK && s.startup_released && s.enabled_time == -10);
            if (!s.phase_paused && !resumed) {
                const double phase = wrap(sign>0 ? PI-asin(x/.6) : asin(x/.6));
                assert(fabs(wrap(s.phase-phase+PI)-PI)<1e-10 && s.torque == 0);
                assert(s.completed_direction == sign && !s.bootstrap_pending && s.gate == 1);
                resumed = true;
            }
        }
        assert(resumed && s.pause_used_this_half);
        for (unsigned i = 1; i <= 24; ++i) {
            now += 10;
            const double x = i==24 ? 0 : sign*(.48-.02*i);
            assert(hip_gentle_update(&s, x/2, -x/2, now, true, &o));
        }
        assert(s.state == HIP_GENTLE_TRACK && s.direction == -sign && !s.pause_used_this_half);
        s = held;
        assert(!hip_gentle_update(&s, NAN, 0, held_ms+10, true, &o));
        assert(!s.phase_paused && s.torque == 0);
        s = held;
        assert(hip_gentle_update(&s, sign*.3, -sign*.3, held_ms+30, true, &o));
        assert(s.reason == HIP_GENTLE_TIMESTAMP_RESET && !s.phase_paused && s.torque == 0);
        s = held; now = held_ms;
        bool stopped = false;
        for (unsigned i = 1; i <= 240; ++i) {
            now += 10;
            // Common motion breaks the long HOLD window; short phase wait cannot stop safety timers.
            const double common = .0015*i;
            assert(hip_gentle_update(&s, sign*.3+common, -sign*.3+common, now, true, &o));
            stopped |= s.state == HIP_GENTLE_EXIT;
            assert(s.torque == 0);
        }
        assert(stopped && s.state == HIP_GENTLE_WAIT && !s.phase_paused);
        // A sub-0.5 Nm target delayed by actual direction still gets a soft handoff.
        hip_gentle_reset(&s); s.motion_ramp_pending = false;
        for (unsigned i = 0; i < 5; ++i) {
            update_motion_permission(&s, -sign, .01);
            apply_motion_target(&s, sign*.49, .01);
        }
        for (unsigned i = 0; i < 3; ++i) {
            update_motion_permission(&s, sign, .01);
            apply_motion_target(&s, sign*.49, .01);
        }
        assert(s.torque == 0 && s.motion_entry_elapsed == .01);
        apply_motion_target(&s, sign*.49, .01);
        assert(s.torque*sign > 0 && fabs(s.torque)<.01);
    }
}

static void test_direction_wait_pause(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.start_angle_rad = .2; p.reversal_speed_rad_s = .2; p.reversal_confirm_s = .03;
    p.low_speed_release_confirm_s = .05;
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_assist_t s;
        hip_gentle_output_t o;
        assert(hip_gentle_init(&s, &p));
        seed(&s, .6, (sign>0 ? PI/2 : 3*PI/2)+.35, TAU);
        s.state = HIP_GENTLE_TRACK; s.direction = sign; s.peak = .4;
        s.startup_first_leg = sign>0 ? 1 : 2; s.startup_released = true;
        s.enabled_time = -10; s.gate = 1; s.completed_direction = sign;
        assert(hip_gentle_update(&s, sign*.2, -sign*.2, 0, true, &o));
        bool paused = false;
        for (unsigned i = 1; i <= 12; ++i) {
            const double x = sign*(.4+.004*i);
            assert(hip_gentle_update(&s, x/2, -x/2, i*10, true, &o));
            assert(s.quiet_time == 0 && s.torque == 0 && s.state == HIP_GENTLE_TRACK);
            if (i<=5) assert(!s.phase_paused);
            paused |= s.phase_paused;
        }
        assert(paused && s.phase_paused && s.pause_used_this_half);
        bool resumed = false;
        for (unsigned i = 1; i <= 15; ++i) {
            const double x = sign*(.448-.01*i);
            assert(hip_gentle_update(&s, x/2, -x/2, (12+i)*10, true, &o));
            if (!s.phase_paused && !resumed) {
                assert(s.torque == 0 && s.startup_released && s.gate == 1);
                resumed = true;
            }
        }
        assert(resumed && s.state == HIP_GENTLE_TRACK);
    }
}

static void test_startup_full_cycle_window(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.start_angle_rad = .2;
    hip_gentle_assist_t s;
    hip_gentle_output_t o;
    assert(hip_gentle_init(&s, &p));
    assert(hip_gentle_update(&s, .05, -.05, 0, true, &o));
    seed(&s, .6, PI/2, TAU); s.state = HIP_GENTLE_START; s.direction = 1; s.peak = .6;
    assert(hip_gentle_update(&s, -.01, .01, 10, true, &o));
    assert(s.state == HIP_GENTLE_TRACK && s.bootstrap_pending && s.bootstrap_halves_remaining == 2);
    s.peak = .6; s.phases[0] = 0; memset(s.rates, 0, sizeof(s.rates));
    assert(hip_gentle_update(&s, .01, -.01, 20, true, &o));
    assert(s.bootstrap_pending && s.bootstrap_halves_remaining == 1);
    for (unsigned i = 1; i <= 100; ++i)
        assert(hip_gentle_update(&s, .15, -.15, 20+i*10, true, &o));
    s.peak = .6; s.phases[0] = 0; memset(s.rates, 0, sizeof(s.rates));
    assert(hip_gentle_update(&s, -.01, .01, 1030, true, &o));
    assert(s.state == HIP_GENTLE_TRACK && s.reason == HIP_GENTLE_STARTUP_RESEED);
    assert(!s.bootstrap_pending && s.bootstrap_halves_remaining == 0);
    assert(o.frequency_hz>.4 && o.frequency_hz<.6);
    s.peak = .6; s.phases[0] = PI; memset(s.rates, 0, sizeof(s.rates));
    assert(hip_gentle_update(&s, .01, -.01, 1040, true, &o));
    assert(s.state == HIP_GENTLE_EXIT && s.reason == HIP_GENTLE_EVENT_MISMATCH);
}

static void test_motion_smoothing(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.slew_nm_s = ASSIST_GENTLE_SLEW_NM_S;
    hip_gentle_assist_t s;
    hip_gentle_output_t o;
    for (int sign = -1; sign <= 1; sign += 2) {
        assert(hip_gentle_init(&s, &p));
        bool first = true;
        double previous = 0;
        for (unsigned i = 0; i < 40; ++i) {
            update_motion_permission(&s, sign, .01);
            apply_motion_target(&s, sign*8, .01);
            assert(fabs(s.torque-previous) <= p.slew_nm_s*.01+1e-12);
            if (s.torque != 0 && first) {
                assert(fabs(s.torque) < .1);
                first = false;
            }
            previous = s.torque;
        }
        assert(!first && s.torque == sign*8);
        for (unsigned i = 0; i < 25; ++i) {
            apply_motion_target(&s, 0, .01);
            assert(s.torque*sign >= 0 && fabs(s.torque) <= fabs(previous));
            assert(fabs(s.torque-previous) <= p.slew_nm_s*.01+1e-12);
            assert(s.target == 0);
            previous = s.torque;
        }
        assert(s.torque == 0 && !s.motion_releasing);
        /* 真正低速停顿后仍要柔和重接，不能按连续换向直接放行。 */
        for (unsigned i = 0; i < 5; ++i) update_motion_permission(&s, 0, .01);
        assert(s.motion_ramp_pending);
        for (unsigned i = 0; i < 3; ++i) {
            update_motion_permission(&s, -sign, .01);
            apply_motion_target(&s, -sign*8, .01);
        }
        assert(s.torque == 0);
        update_motion_permission(&s, -sign, .01);
        apply_motion_target(&s, -sign*8, .01);
        assert(s.torque*sign < 0 && fabs(s.torque) < .1);
        hip_gentle_reset(&s);
        s.torque = sign*8;
        apply_motion_target(&s, 0, .01);
        previous = s.torque;
        assert(s.motion_releasing);
        stop(&s, HIP_GENTLE_WEAK_PEAK, .01);
        assert(!s.motion_releasing && fabs(previous-s.torque) == .5);
        assert(!hip_gentle_update(&s, NAN, 0, 10, true, &o));
        assert(o.left_torque_nm == 0 && !s.motion_releasing);
    }
    p.motion_fall_s = .21;
    assert(!hip_gentle_init(&s, &p));
    p = *base; p.motion_rise_s = 0;
    assert(!hip_gentle_init(&s, &p));
}

static void test_long_hold_resume(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.amplitude_nm = 4; p.torque_limit_nm = 8; p.extension_torque_ratio = 1;
    p.start_angle_rad = .2; p.quiet_speed_rad_s = .1; p.quiet_duration_s = 2;
    p.hold_confirm_s = ASSIST_GENTLE_HOLD_CONFIRM_S; p.reversal_speed_rad_s = .2;
    p.inner_rise_s = .6; p.outer_rise_s = .45; p.peak_drop_ratio = .65;
    p.weak_peak_stop_enabled = true;
    const unsigned hold_samples[] = {500, 3000, 6000};
    for (int sign = -1; sign <= 1; sign += 2) {
        for (unsigned duration = 0; duration < 3; ++duration) {
            hip_gentle_assist_t s;
            hip_gentle_output_t o;
            assert(hip_gentle_init(&s, &p));
            uint32_t now = 0;
            for (unsigned i = 0; i <= 425; ++i) {
                const double q = sign*.35*sin(TAU*i*.01);
                now = i*10;
                assert(hip_gentle_update(&s, q, -q, now, true, &o));
            }
            assert(s.state == HIP_GENTLE_TRACK && s.startup_released);
            const uint8_t first_leg = s.startup_first_leg;
            bool held = false;
            double held_phase = 0, held_age = 0;
            for (unsigned i = 0; i < hold_samples[duration]; ++i) {
                now += 10;
                assert(hip_gentle_update(&s, sign*.35, -sign*.35, now, true, &o));
                if (s.state == HIP_GENTLE_HOLD && !held) {
                    held = true; held_phase = s.phase; held_age = s.time_s-s.event_time;
                }
                if (held) {
                    assert(s.state == HIP_GENTLE_HOLD && s.phase == held_phase);
                    assert(fabs((s.time_s-s.event_time)-held_age) < 1e-7);
                    assert(o.left_torque_nm == 0 && o.right_torque_nm == 0);
                }
                assert(s.history_count < HIP_GENTLE_HISTORY_CAPACITY);
            }
            assert(held);
            const hip_gentle_assist_t held_state = s;
            assert(!hip_gentle_update(&s, sign*.35, -sign*.35, now+10, false, &o));
            assert(s.state == HIP_GENTLE_WAIT && o.left_torque_nm == 0);
            s = held_state;
            assert(hip_gentle_update(&s, sign*.35, -sign*.35, now+30, true, &o));
            assert(s.state == HIP_GENTLE_WAIT && o.reason == HIP_GENTLE_TIMESTAMP_RESET && o.left_torque_nm == 0);
            s = held_state;
            /* 小抖动不能恢复出力。 */
            for (unsigned i = 0; i < 20; ++i) {
                now += 10;
                const double q = sign*(.35+.015*sin(TAU*i*.1));
                assert(hip_gentle_update(&s, q, -q, now, true, &o));
                assert(s.state == HIP_GENTLE_HOLD && o.left_torque_nm == 0);
            }
            bool resumed = false, assisted = false;
            double previous = 0;
            const unsigned fall_samples = duration == 2 ? 100U : 50U;
            for (unsigned i = 1; i <= fall_samples; ++i) {
                now += 10;
                const double q = i == fall_samples ? 0 : sign*.35*cos(PI/2*i/fall_samples);
                assert(hip_gentle_update(&s, q, -q, now, true, &o));
                if (s.state != HIP_GENTLE_HOLD && !resumed) {
                    assert(s.state == HIP_GENTLE_TRACK && o.left_torque_nm == 0);
                    resumed = true;
                }
                assert(s.state == HIP_GENTLE_HOLD || s.state == HIP_GENTLE_TRACK);
                assert(s.startup_released && s.startup_first_leg == first_leg);
                assert(o.left_torque_nm*sign <= 0 && fabs(o.left_torque_nm) <= 8);
                assert(fabs(o.left_torque_nm-previous) <= p.slew_nm_s*.01+1e-9);
                assisted |= fabs(o.left_torque_nm) > .01;
                previous = o.left_torque_nm;
            }
            assert(resumed && assisted);
            assert(s.recent_peak_count[0] == 0 && s.recent_peak_count[1] == 0);
            for (unsigned i = 1; i <= HIP_GENTLE_RESUME_SYNC_HALVES*50U+10U; ++i) {
                now += 10;
                const double q = -sign*.35*sin(TAU*i*.01);
                assert(hip_gentle_update(&s, q, -q, now, true, &o));
                assert(s.state == HIP_GENTLE_TRACK);
            }
            assert(s.resume_sync_remaining == 0 && !s.bootstrap_pending);
            const hip_gentle_assist_t tracked = s;
            s.direction = 1; s.peak = .7; s.previous_x = .1;
            s.phases[0] = 0; memset(s.rates, 0, sizeof(s.rates));
            assert(hip_gentle_update(&s, -.05, .05, now+10, true, &o));
            assert(s.state == HIP_GENTLE_EXIT && s.reason == HIP_GENTLE_EVENT_MISMATCH);
            s = tracked;
            /* 回到中立站稳仍结束步态，不能无限保留所有停止情况。 */
            for (unsigned i = 0; i < 260; ++i) {
                now += 10;
                assert(hip_gentle_update(&s, 0, 0, now, true, &o));
            }
            assert(s.state == HIP_GENTLE_WAIT && o.left_torque_nm == 0);
        }
        /* 首腿悬停可保留采集进度，但不能绕过首次两腿确认。 */
        hip_gentle_assist_t s;
        hip_gentle_output_t o;
        assert(hip_gentle_init(&s, &p));
        uint32_t now = 0;
        for (unsigned i = 0; i <= 40; ++i) {
            const double q = .65*sin(PI/2*i/40);
            now = i*10;
            assert(hip_gentle_update(&s, sign>0 ? q:0, sign<0 ? q:0, now, true, &o));
        }
        for (unsigned i = 0; i < 500; ++i) {
            now += 10;
            assert(hip_gentle_update(&s, sign>0 ? .65:0, sign<0 ? .65:0, now, true, &o));
        }
        assert(s.state == HIP_GENTLE_HOLD && s.pause_from_state == HIP_GENTLE_COLLECT);
        for (unsigned i = 1; i <= 50; ++i) {
            now += 10;
            const double q = i==50 ? 0:.65*cos(PI/2*i/50);
            assert(hip_gentle_update(&s, sign>0 ? q:0, sign<0 ? q:0, now, true, &o));
            assert(!s.startup_released && o.left_torque_nm == 0);
        }
        assert(s.state == HIP_GENTLE_START || s.state == HIP_GENTLE_TRACK);
        assert(!hip_gentle_update(&s, 0, 0, now+10, false, &o));
        assert(s.state == HIP_GENTLE_WAIT && o.left_torque_nm == 0);
    }
    hip_gentle_assist_t invalid_state;
    p.hold_confirm_s = p.low_speed_release_confirm_s;
    assert(!hip_gentle_init(&invalid_state, &p));
    p.hold_confirm_s = p.quiet_duration_s;
    assert(!hip_gentle_init(&invalid_state, &p));
}

static void test_slow_progress(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.start_angle_rad = .2;
    p.min_frequency_hz = ASSIST_GENTLE_MIN_FREQUENCY_HZ;
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_assist_t s;
        hip_gentle_output_t o;
        assert(hip_gentle_init(&s, &p));
        seed(&s, .2, sign > 0 ? PI/2 : 3*PI/2, TAU*p.min_frequency_hz);
        s.state = HIP_GENTLE_TRACK; s.direction = sign; s.peak = .2;
        assert(hip_gentle_update(&s, sign*.1, -sign*.1, 0, true, &o));
        /* 连续动作四秒不需要过零；旧事件计时已超时仍不能退出。 */
        for (unsigned i = 1; i <= 400; ++i) {
            const double x = sign*(.2+.0025*i);
            assert(hip_gentle_update(&s, x/2, -x/2, i*10, true, &o));
            assert(s.state == HIP_GENTLE_TRACK);
            assert(s.time_s-s.progress_time < p.progress_timeout_s);
        }
        assert(s.time_s-s.event_time > p.progress_timeout_s);
        /* 累计路程很大的原地往复也不能反复充当本半波的新进展。 */
        bool stopped = false;
        for (unsigned i = 401; i <= 700; ++i) {
            const double x = sign*(1.2+.04*sin(TAU*4*(i-400)*.01));
            assert(hip_gentle_update(&s, x/2, -x/2, i*10, true, &o));
            if (s.state == HIP_GENTLE_EXIT) {
                assert(s.reason == HIP_GENTLE_NO_PROGRESS);
                stopped = true;
                break;
            }
        }
        assert(stopped);
    }
}

static void test_interrupted_pulse(const hip_gentle_parameters_t *base)
{
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_assist_t s;
        assert(hip_gentle_init(&s, base));
        s.torque = sign*.05; s.completed_direction = -sign;
        s.motion_waiting = true;
        for (unsigned i = 0; i < 30; ++i) apply_motion_target(&s, sign*4, .01);
        assert(s.torque == 0 && s.paused_direction == sign && s.completed_direction == -sign);
        for (unsigned i = 0; i < 40; ++i) {
            s.previous_x = sign*.02*sin(TAU*i*.1);
            update_motion_permission(&s, sign, .01);
            apply_motion_target(&s, sign*4, .01);
            assert(s.torque == 0);
        }
        s.previous_x = sign*.06;
        apply_motion_target(&s, sign*4, .01);
        assert(s.paused_direction == 0 && s.torque == 0);
        apply_motion_target(&s, sign*4, .01);
        assert(s.torque*sign > 0 && fabs(s.torque) < .1);
        for (unsigned i = 0; i < 30; ++i) apply_motion_target(&s, sign*4, .01);
        assert(s.torque == sign*4);
        for (unsigned i = 0; i < 30; ++i) apply_motion_target(&s, 0, .01);
        assert(s.torque == 0 && s.completed_direction == sign && s.paused_direction == 0);
        s.previous_x = sign*.2;
        for (unsigned i = 0; i < 30; ++i) apply_motion_target(&s, sign*4, .01);
        assert(s.torque == 0); // 正常相位段已结束，同向波形反弹仍不能重启。
    }
}

static void test_hold_micro_motion(const hip_gentle_parameters_t *base)
{
    hip_gentle_parameters_t p = *base;
    p.start_angle_rad = .2; p.quiet_duration_s = 2;
    p.hold_confirm_s = 1; p.hold_range_rad = .08;
    for (unsigned break_window = 0; break_window < 2; ++break_window) {
        hip_gentle_assist_t s;
        hip_gentle_output_t o;
        assert(hip_gentle_init(&s, &p));
        seed(&s, .6, PI/2, TAU);
        s.phase = PI/2; s.state = HIP_GENTLE_TRACK; s.direction = 1; s.peak = .6;
        s.startup_first_leg = 1; s.startup_released = true; s.time_s = 1.5;
        assert(hip_gentle_update(&s, .3, -.3, 0, true, &o));
        s.progress_time = 0; s.progress_x = .6;
        for (unsigned i = 1; i <= 110; ++i) {
            const double q = break_window && i >= 85 ? .2 : .3+.008*sin(TAU*8*i*.01);
            const double right = break_window && i >= 85 ? -.4 : -q;
            assert(hip_gentle_update(&s, q, right, i*10, true, &o));
            if (!break_window) assert(s.state == HIP_GENTLE_TRACK || s.state == HIP_GENTLE_HOLD);
            if (i >= 80) assert(o.left_torque_nm == 0 && o.right_torque_nm == 0);
            if (break_window && i == 85) {
                assert(s.state == HIP_GENTLE_EXIT && s.reason == HIP_GENTLE_NO_PROGRESS);
                break;
            }
        }
        if (!break_window) assert(s.state == HIP_GENTLE_HOLD);
    }
    hip_gentle_assist_t s;
    p.hold_range_rad = 0;
    assert(!hip_gentle_init(&s, &p));
    p.hold_range_rad = p.start_angle_rad;
    assert(!hip_gentle_init(&s, &p));
}

int main(void)
{
    const hip_gentle_parameters_t p = {
        .amplitude_nm=2, .torque_limit_nm=8, .start_angle_rad=.5, .derivative_window_s=.04,
        .extension_torque_ratio=.5,
        .startup_lift_angle_rad=.05, .startup_lift_speed_rad_s=.2, .startup_lift_confirm_s=.05,
        .min_frequency_hz=.3, .max_frequency_hz=2.5, .quiet_speed_rad_s=.1,
        .quiet_duration_s=2, .peak_drop_ratio=.5, .weak_half_cycles=2, .progress_timeout_s=2,
        .hold_confirm_s=ASSIST_GENTLE_HOLD_CONFIRM_S, .hold_range_rad=ASSIST_GENTLE_HOLD_RANGE_RAD,
        .low_speed_release_confirm_s=.05,
        .reversal_speed_rad_s=.35, .reversal_confirm_s=.03,
        .motion_rise_s=ASSIST_GENTLE_MOTION_RISE_S, .motion_fall_s=ASSIST_GENTLE_MOTION_FALL_S,
        .tracking_error_stop_enabled=false, .tracking_error_rad=.5, .tracking_error_duration_s=.3,
        .phase_error_limit_rad=PI/2, .phase_correction_rad_s=1, .slew_nm_s=40,
        .reset_gap_s=.02, .width=1.04, .asymmetry=-.02, .response=65,
        .inner_rise_s=2, .outer_rise_s=1.5, .amplitude_gain=2.5, .phase_gain=1.6,
        .frequency_gain=3.2, .bias_gain=.5
    };
    hip_gentle_assist_t s = {0};
    hip_gentle_output_t o;
    assert(!hip_gentle_init(&s, NULL));
    assert(hip_gentle_init(&s, &p));
    s.torque = s.target = 8;
    output(&s, 0, &o);
    assert(o.left_torque_nm == 8 && o.right_torque_nm == -4);
    assert(o.target_left_torque_nm == 8 && o.target_right_torque_nm == -4);
    s.torque = s.target = -8;
    output(&s, 0, &o);
    assert(o.left_torque_nm == -4 && o.right_torque_nm == 8);
    assert(o.target_left_torque_nm == -4 && o.target_right_torque_nm == 8);
    hip_gentle_reset(&s);
    /* 40ms低速不截断当前半波；恢复即清计时，重新连续满50ms才撤销。 */
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_reset(&s);
        for (unsigned i = 0; i < 30; ++i) {
            update_motion_permission(&s, sign, .01);
            apply_motion_target(&s, sign*2, .01);
        }
        for (unsigned i = 0; i < 4; ++i) {
            update_motion_permission(&s, 0, .01);
            apply_motion_target(&s, sign*2, .01);
            assert(!s.motion_waiting && s.torque == sign*2);
        }
        update_motion_permission(&s, sign*.2, .01);
        apply_motion_target(&s, sign*2, .01);
        assert(s.quiet_time == 0 && !s.motion_waiting && s.torque == sign*2);
        for (unsigned i = 0; i < 5; ++i) {
            update_motion_permission(&s, 0, .01);
            apply_motion_target(&s, sign*2, .01);
        }
        const double duration = fmax(p.motion_fall_s, 1.5*2/p.slew_nm_s);
        assert(s.motion_waiting && s.target == 0);
        assert(fabs(s.torque-sign*2*(1-smooth_fraction(.01/duration))) < 1e-12);
    }
    /* 回零后同向不可重启；低速/中间区不解锁，正负方向完全对称。 */
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_reset(&s);
        for (unsigned i = 0; i < 30; ++i) {
            update_motion_permission(&s, sign, .01);
            apply_motion_target(&s, sign*2, .01);
        }
        assert(s.torque == sign*2);
        for (unsigned i = 0; i < 30; ++i) apply_motion_target(&s, 0, .01);
        assert(s.torque == 0 && s.completed_direction == sign);
        for (unsigned i = 0; i < 10; ++i) {
            update_motion_permission(&s, i < 3 ? 0 : i < 6 ? sign*.2 : sign, .01);
            apply_motion_target(&s, sign*2, .01);
            assert(s.torque == 0 && s.target == 0 && s.completed_direction == sign);
        }
        for (unsigned i = 0; i < 2; ++i) {
            update_motion_permission(&s, -sign, .01);
            apply_motion_target(&s, -sign*2, .01);
            assert(s.torque == 0);
        }
        update_motion_permission(&s, -sign, .01);
        apply_motion_target(&s, -sign*2, .01);
        assert(s.torque == 0); // 大目标迟到接入先从零渐入。
        update_motion_permission(&s, -sign, .01);
        apply_motion_target(&s, -sign*2, .01);
        assert(s.torque*sign < 0 && fabs(s.torque) < .1);
        stop(&s, HIP_GENTLE_QUIET, .01);
        assert(s.completed_direction == 0);
        hip_gentle_reset(&s);
        assert(s.completed_direction == 0);
    }
    for (int sign = -1; sign <= 1; sign += 2) {
        hip_gentle_reset(&s);
        assert(hip_gentle_update(&s, 0, 0, 0, true, &o));
        assert(hip_gentle_update(&s, sign*.5, 0, 10, true, &o) && o.state == HIP_GENTLE_WAIT);
        assert(hip_gentle_update(&s, sign*.51, 0, 20, true, &o) && o.state == HIP_GENTLE_COLLECT);
        assert(s.direction == sign);
    }
    hip_gentle_reset(&s);
    double previous = 0, maximum = 0;
    unsigned states = 0;
    for (unsigned i = 0; i <= 1500; ++i) {
        const double angle = i < 800 ? .35*sin(TAU*i*.01) : 0;
        assert(hip_gentle_update(&s, angle, -angle, i*10, true, &o));
        states |= 1U << o.state;
        assert(fabs(o.left_torque_nm) <= 2.000001 && fabs(o.right_torque_nm) <= 2.000001);
        if (o.left_torque_nm > 0) assert(o.right_torque_nm == -o.left_torque_nm*.5);
        else assert(o.left_torque_nm == -o.right_torque_nm*.5);
        assert(fabs(o.left_torque_nm-previous) <= .400001);
        previous = o.left_torque_nm;
        maximum = fmax(maximum, fabs(previous));
        if (i > 1200) assert(previous == 0);
    }
    assert((states & (1U << HIP_GENTLE_TRACK)) && (states & (1U << HIP_GENTLE_EXIT)) && maximum > .1);
    assert(hip_gentle_update(&s, 0, 0, 15030, true, &o) && o.reason == HIP_GENTLE_TIMESTAMP_RESET && o.left_torque_nm == 0);
    assert(!hip_gentle_update(&s, NAN, 0, 15040, true, &o) && o.left_torque_nm == 0);
    hip_gentle_reset(&s);
    s.parameters.weak_peak_stop_enabled = true;
    assert(hip_gentle_update(&s, .15, -.15, 0, true, &o));
    s.state = HIP_GENTLE_TRACK; s.direction = 1; s.peak = .3; s.last_peak = .8;
    s.torque = 2;
    seed(&s, .8, PI/2, TAU);
    assert(hip_gentle_update(&s, .1, -.1, 10, true, &o));
    assert(o.state == HIP_GENTLE_TRACK && s.weak_pending);
    assert(o.target_left_torque_nm == 0 && fabs(o.left_torque_nm-1.6) < 1e-12);
    assert(hip_gentle_update(&s, -.01, .01, 20, true, &o));
    assert(o.state == HIP_GENTLE_EXIT && o.reason == HIP_GENTLE_WEAK_PEAK);
    assert(hip_gentle_update(&s, 0, 0, UINT32_MAX-4U, true, &o));
    assert(hip_gentle_update(&s, .01, 0, 5, true, &o));
    assert(fabs(s.time_s-.01) < 1e-12 && s.timestamp_ms == 5 && s.history_count == 2);
    hip_gentle_reset(&s);
    for (unsigned i = 0; i < 1000; ++i) assert(hip_gentle_update(&s, 0, 0, i, true, &o));
    assert(s.history_count <= 42);
    for (unsigned i = 0; i < 100; ++i) {
        const double x = TAU*i/100;
        const double expected = (x < PI ? -1 : 1)*sin(x)*sin(x);
        assert(fabs(hip_gentle_waveform(x+PI/2, 1, 0)-expected) < 1e-12);
    }
    hip_gentle_parameters_t invalid = p;
    invalid.width = 0;
    assert(!hip_gentle_init(&s, &invalid));
    test_swing_tiers(&p);
    test_alternating_startup(&p);
    test_small_slow_steps(&p);
    test_weak_peak_rearm(&p);
    test_double_peak(&p);
    test_pacing_recovery(&p);
    test_motion_smoothing(&p);
    test_torque_handoff(&p);
    test_phase_paused(&p);
    test_direction_wait_pause(&p);
    test_startup_full_cycle_window(&p);
    test_long_hold_resume(&p);
    test_hold_micro_motion(&p);
    test_slow_progress(&p);
    test_interrupted_pulse(&p);
    puts("PASS C hip gentle: direction-wait pause, bounded startup cycle, unified resume, timers, faults");
    return 0;
}
#endif
