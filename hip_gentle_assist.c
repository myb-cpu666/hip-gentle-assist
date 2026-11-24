#include "hip_gentle_assist.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#define HIP_PI            3.14159265358979323846f
#define HIP_TWO_PI        6.28318530717958647692f
#define HIP_HALF_PI       1.57079632679489661923f
#define HIP_THREE_HALF_PI 4.71238898038468985769f
#define HIP_RAD_TO_DEG    57.2957795130823208768f
#define HIP_DEG_TO_RAD    0.01745329251994329577f

/* ==================== 基础数学与圆周插值工具 ==================== */

static float clampf(float x, float lo, float hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

static void ramp_to(float *cur, float target, float step)
{
    if (*cur < target) {
        *cur += step;
        if (*cur >= target) {
            *cur = target;
        }
    } else if (*cur > target) {
        *cur -= step;
        if (*cur <= target) {
            *cur = target;
        }
    }
}

static float wrap_2pi(float x)
{
    float r = fmodf(x, HIP_TWO_PI);
    if (r < 0.0f) {
        r += HIP_TWO_PI;
    }
    return r == HIP_TWO_PI ? 0.0f : r; // 极小负数加 2π 可能舍入到上边界
}

static float wrap_pm_pi(float x)
{
    return wrap_2pi(x + HIP_PI) - HIP_PI;
}

static float circular_lerp(float a, float b, float alpha)
{
    return wrap_2pi(a + alpha * wrap_pm_pi(b - a));
}

static float asin_safe(float u)
{
    return asinf(clampf(u, -1.0f, 1.0f));
}

/* §4.1 S₁ = 9/8·(sin φ + sin 3φ / 9) = 1.5 sin φ − 0.5 sin³φ。 */
static float afo_reconstruct(const hip_gentle_afo_t *afo, float basis[5])
{
    float sn = sinf(afo->phi[0]);
    basis[0] = sn * (1.5f - 0.5f * sn * sn);
    float estimate = afo->alpha0 + afo->alpha[0] * basis[0];
    for (unsigned k = 1; k < 5; ++k) {
        basis[k] = sinf(afo->phi[k]);
        estimate += afo->alpha[k] * basis[k];
    }
    return estimate;
}

static void afo_refresh_rates(hip_gentle_afo_t *a, const hip_gentle_afo_gains_t *g, float signal)
{
    float basis[5];
    a->y_hat = afo_reconstruct(a, basis);
    a->phase_out = wrap_2pi(a->phi[0] + a->phase_offset_rad);
    float error = signal - a->y_hat;
    float cs = cosf(a->phi[0]);
    float slope = 1.5f * cs * cs * cs;
    a->d_alpha0 = g->eta_0 * error;
    a->d_omega = g->gamma_omega * error * slope;
    a->d_alpha[0] = g->eta_alpha * error * basis[0];
    a->d_phi[0] = fmaxf(0.0f, a->omega + g->nu_phi * error * slope);
    for (unsigned k = 1; k < 5; ++k) {
        a->d_alpha[k] = 0.5f * g->eta_alpha * error * basis[k];
        a->d_phi[k] = fmaxf(0.0f, (k+1) * a->omega + g->nu_phi * error * cosf(a->phi[k]));
    }
}

static void afo_reset(hip_gentle_afo_t *a)
{
    float offset = a->phase_offset_rad;
    *a = (hip_gentle_afo_t){.omega = HIP_PI, .phase_offset_rad = offset};
}

static void afo_hot_start(hip_gentle_assist_t *s, float omega, float peak_rad,
                          float fundamental_phase, float signal_rad)
{
    hip_gentle_afo_t *a = &s->afo;
    afo_reset(a);
    a->enabled = true;
    a->omega = a->omega_filt = clampf(omega, HIP_GENTLE_AFO_MIN_OMEGA_RAD_S, HIP_GENTLE_AFO_MAX_OMEGA_RAD_S);
    a->alpha[0] = peak_rad;
    a->phi[0] = fundamental_phase;
    a->velocity_filt_rad_s = a->omega * peak_rad;
    afo_refresh_rates(a, &s->afo_gains, signal_rad);
}

static void afo_step(hip_gentle_afo_t *a, const hip_gentle_afo_gains_t *g, float signal, float dt)
{
    if (!a->enabled) {
        a->phase_out = a->omega_filt = a->velocity_filt_rad_s = 0.0f;
        return;
    }
    const float h = dt / HIP_GENTLE_AFO_SUBSTEPS, weight = h / HIP_GENTLE_AFO_FILTER_TAU_S;
    for (unsigned i = 0; i < HIP_GENTLE_AFO_SUBSTEPS; ++i) {
        a->omega = clampf(a->omega + h * a->d_omega, HIP_GENTLE_AFO_MIN_OMEGA_RAD_S, HIP_GENTLE_AFO_MAX_OMEGA_RAD_S);
        a->alpha0 += h * a->d_alpha0;
        for (unsigned k=0; k<5; ++k) {
            a->alpha[k] += h * a->d_alpha[k];
            a->phi[k] = wrap_2pi(a->phi[k] + h * a->d_phi[k]);
        }
        afo_refresh_rates(a, g, signal);
        a->omega_filt += weight * (a->omega - a->omega_filt);
        a->velocity_filt_rad_s += weight * (a->omega*a->alpha[0] - a->velocity_filt_rad_s);
    }
}

static void afo_sync_phase(hip_gentle_assist_t *s, float target, float signal_rad)
{
    hip_gentle_afo_t *a = &s->afo;
    float delta = wrap_pm_pi(target - a->phase_out);
    for (unsigned k=0; k<5; ++k)
        a->phi[k] = wrap_2pi(a->phi[k] + (k+1)*delta);
    afo_refresh_rates(a, &s->afo_gains, signal_rad);
}

static void fsm_reset_to_idle(hip_gentle_fsm_t *fsm, hip_gentle_afo_t *afo, double now)
{
    afo_reset(afo);
    memset(fsm, 0, sizeof(*fsm));
    fsm->state = HIP_GENTLE_IDLE;
    fsm->evt = fsm->prev_evt = 8;
    fsm->evt_ts[8] = fsm->state_enter_ts = now;
    fsm->evt_ts_valid[8] = true;
    fsm->half_period_s = 0.5f;
    fsm->static_l_since = fsm->static_r_since = -1.0;
    fsm->still_l_since = fsm->still_r_since = -1.0;
    fsm->scale_timeout = fsm->scale_env = fsm->scale_energy = 1.0f;
    fsm->dyn_to_coeff = HIP_GENTLE_DYN_TO_FAST_COEFF;
}

static void fsm_enter_cooldown(hip_gentle_fsm_t *fsm, hip_gentle_afo_t *afo, double now)
{
    afo_reset(afo);
    fsm->state = HIP_GENTLE_COOLDOWN;
    fsm->static_l_since = fsm->static_r_since = -1.0; // 冷却用更窄的停稳带重新计时
    fsm->still_l_since = fsm->still_r_since = -1.0;
    fsm->peak_fade_s = 0.0f;
    fsm->phase_override_valid = false;
    fsm->dir_sign = 0.0f;
    fsm->recent_peak_rad = 0.0f;
    fsm->phase_amp_filt_rad = 0.0f;
    fsm->scale_energy = 1.0f;
    fsm->abrupt_shrink_resonance = false;
    fsm->collapse_halves = 0U;
    fsm->early_peak_armed = false;
    fsm->early_peak_sign = 0;
    fsm->early_peak_deg = 0.0f;
    fsm->early_peak_ref_deg = 0.0f;
    fsm->last_peak_deg[0] = fsm->last_peak_deg[1] = 0.0f;
    fsm->cooldown_min_s = 0.0f;
    fsm->state_enter_ts = now;
}

/* 8→0/4 是跨越起步阈值，不能用作真实过零点。直接在角差采样间插值，
 * 特征点后续完成 3→4 / 7→0 时再采信，容纳一帧跨过阈值和零点的情况。 */
static void observe_zero(hip_gentle_fsm_t *fsm, float y, double now, float dt)
{
    if (!fsm->prev_y_valid) return;
    float prev = fsm->prev_y;
    if (!((prev <= 0.0f && y > 0.0f) || (prev >= 0.0f && y < 0.0f))) return;
    unsigned side = y < 0.0f;
    fsm->observed_zero_ts[side] = now - dt * (double)(y / (y - prev));
    fsm->observed_zero_valid[side] = true;
}

static float observe_half_period(hip_gentle_fsm_t *fsm, int sign)
{
    unsigned side = sign < 0;
    fsm->abrupt_shrink_resonance = false;
    if (!fsm->observed_zero_valid[side]) return 0.0f;
    float measured = 0.0f;
    double cross = fsm->observed_zero_ts[side];
    if (fsm->half_start_sign == -sign && cross > fsm->half_start_ts) {
        measured = (float)(cross - fsm->half_start_ts);
        /* 过正零点刚完成负半波，过负零点刚完成正半波。 */
        float peak = sign > 0 ? -fsm->seed_peak_neg : fsm->seed_peak_pos;
        unsigned ref_idx = sign > 0 ? 1U : 0U;
        float *reference = &fsm->reference_peak_deg[ref_idx];
        float eff_ref = fmaxf(*reference, fsm->last_peak_deg[ref_idx]);
        bool shrinking = *reference > HIP_GENTLE_START_THRESHOLD_DEG &&
            peak < HIP_GENTLE_SHRINK_PEAK_MAX_DEG && peak < HIP_GENTLE_SHRINK_PEAK_RATIO * *reference;
        bool contracting = measured < HIP_GENTLE_SHORT_HALF_S ||
            measured < HIP_GENTLE_CONTRACT_HALF_RATIO * fsm->reference_half_s;
        fsm->abrupt_shrink_resonance = *reference > HIP_GENTLE_START_THRESHOLD_DEG &&
            fsm->reference_half_s > 0.0f &&
            peak < HIP_GENTLE_ABRUPT_SHRINK_PEAK_DEG &&
            peak < HIP_GENTLE_ABRUPT_SHRINK_PEAK_RATIO * *reference &&
            measured < HIP_GENTLE_ABRUPT_SHRINK_HALF_RATIO * fsm->reference_half_s;
        if (shrinking && contracting)
            ++fsm->short_halves;
        else
            fsm->short_halves = 0;
        if (eff_ref >= HIP_GENTLE_EARLY_STOP_MIN_REF_DEG &&
            peak < HIP_GENTLE_EARLY_STOP_MAX_PEAK_DEG &&
            peak < HIP_GENTLE_COLLAPSE_PEAK_RATIO * eff_ref)
            ++fsm->collapse_halves;
        else
            fsm->collapse_halves = 0U;
        if (!fsm->short_halves && !fsm->collapse_halves && peak > HIP_GENTLE_START_THRESHOLD_DEG) {
            *reference = *reference > 0.0f ? *reference + HIP_GENTLE_REF_ADAPT_ALPHA * (peak - *reference) : peak;
            fsm->reference_half_s = fsm->reference_half_s > 0.0f
                ? fsm->reference_half_s + HIP_GENTLE_REF_ADAPT_ALPHA * (measured - fsm->reference_half_s) : measured;
        }

    }
    fsm->half_start_ts = cross;
    fsm->half_start_sign = sign;
    fsm->observed_zero_valid[side] = false;
    return measured;
}

static void fsm_update_dy_ma(hip_gentle_fsm_t *fsm, float y_hip_deg, float dt)
{
    if (dt <= 0.0f) return;
    if (!fsm->prev_y_valid) {
        fsm->prev_y = y_hip_deg;
        fsm->prev_y_valid = true;
        fsm->prev_dy_ma = fsm->dy_ma;
        return;
    }

    float dy_raw = (y_hip_deg - fsm->prev_y) / dt;
    fsm->ma_buf[fsm->ma_idx] = dy_raw;
    fsm->ma_idx = (uint8_t)((fsm->ma_idx + 1U) % (unsigned int)HIP_GENTLE_DY_WINDOW);

    float sum = 0.0f;
    for (unsigned int i = 0U; i < (unsigned int)HIP_GENTLE_DY_WINDOW; i++) {
        sum += fsm->ma_buf[i];
    }
    fsm->prev_dy_ma = fsm->dy_ma;
    fsm->dy_ma = sum / (float)HIP_GENTLE_DY_WINDOW;
}

static void update_quiet_since(double *since, float speed_deg_s, float limit_deg_s, double now)
{
    if (fabsf(speed_deg_s) > limit_deg_s) *since = -1.0;
    else if (*since < 0.0) *since = now;
}

static bool both_since_for(double l_since, double r_since, double now, float duration_s)
{
    return l_since >= 0.0 && r_since >= 0.0 &&
           now - fmax(l_since, r_since) + 1e-6 >= duration_s;
}

static void fsm_update_static_and_dynto(hip_gentle_fsm_t *fsm,
                                        float left_deg_s,
                                        float right_deg_s,
                                        double now, float dt)
{
    float limit_rad_s = fsm->state == HIP_GENTLE_COOLDOWN ? HIP_GENTLE_SETTLED_SPEED_RAD_S : HIP_GENTLE_STATIC_SPEED_RAD_S;
    float limit_deg_s = limit_rad_s * HIP_RAD_TO_DEG;
    update_quiet_since(&fsm->static_l_since, left_deg_s, limit_deg_s, now);
    update_quiet_since(&fsm->static_r_since, right_deg_s, limit_deg_s, now);
    float still_deg_s = HIP_GENTLE_STILL_SPEED_RAD_S * HIP_RAD_TO_DEG;
    update_quiet_since(&fsm->still_l_since, left_deg_s, still_deg_s, now);
    update_quiet_since(&fsm->still_r_since, right_deg_s, still_deg_s, now);

    /* 第 5.2 节动态超时系数 k_to 更新 */
    float v_act = fmaxf(fmaxf(fabsf(left_deg_s), fabsf(right_deg_s)), fabsf(fsm->dy_ma));
    bool at_peak = fsm->evt == 2 || fsm->evt == 6;
    float target_kto = (at_peak || v_act >= HIP_GENTLE_DYN_TO_SPEED_DEG_S)
        ? HIP_GENTLE_DYN_TO_FAST_COEFF : HIP_GENTLE_DYN_TO_SLOW_COEFF;
    /* 达峰帧先切入保持档，避免把峰前的低速系数锁进整段回程。 */
    if (at_peak && fsm->evt != fsm->prev_evt)
        fsm->dyn_to_coeff = target_kto;
    float alpha_kto = dt / (dt + HIP_GENTLE_DYN_TO_TAU_S);
    fsm->dyn_to_coeff += alpha_kto * (target_kto - fsm->dyn_to_coeff);
}

static bool fsm_both_static_for(const hip_gentle_fsm_t *f, double now, float duration_s)
{
    return both_since_for(f->static_l_since, f->static_r_since, now, duration_s);
}

static bool fsm_both_still_for(const hip_gentle_fsm_t *f, double now, float duration_s)
{
    return both_since_for(f->still_l_since, f->still_r_since, now, duration_s);
}

static bool fsm_narrow_stopped(const hip_gentle_fsm_t *f, const hip_gentle_afo_t *a,
                               float y_deg, double now)
{
    float amp_deg = fabsf(a->alpha[0]) * HIP_RAD_TO_DEG;
    bool near_mid = fabsf(y_deg) <= HIP_GENTLE_NARROW_STOP_DEG &&
                    fabsf(y_deg) <= HIP_GENTLE_STOP_MID_AMP_RATIO * amp_deg;
    bool high_cadence = a->enabled && a->omega_filt >= HIP_GENTLE_FAST_STOP_OMEGA_RAD_S;
    float fast_s = clampf(HIP_GENTLE_NARROW_STOP_HALF_RATIO * f->half_period_s,
                          HIP_GENTLE_NARROW_STOP_MIN_S, HIP_GENTLE_NARROW_STOP_STATIC_S);
    if (high_cadence && fsm_both_still_for(f, now, fast_s))
        return true;
    if (fabsf(y_deg) <= HIP_GENTLE_NARROW_STOP_DEG)
        return fsm_both_still_for(f, now, near_mid ? fast_s : HIP_GENTLE_NARROW_STOP_STATIC_S);
    return false;
}

static bool fsm_track_peak_budget(hip_gentle_fsm_t *fsm,
                                  const hip_gentle_afo_t *afo,
                                  float signed_y,
                                  float signed_dy,
                                  float signed_prev_dy,
                                  float signed_rev_vel_m,
                                  float thr_mag,
                                  float margin_mag,
                                  float rev_vel_thresh,
                                  int sign,
                                  float *max_dy_mag,
                                  float *peak_mag)
{
    if (signed_y > *peak_mag) {
        *peak_mag = signed_y;
    }

    if (fsm->rapid_channel == (sign>0 ? HIP_GENTLE_RAPID_POSITIVE : HIP_GENTLE_RAPID_NEGATIVE))
        return signed_y < 0.0f && signed_rev_vel_m > rev_vel_thresh;

    if (fsm->geo_decel_init == 0.0f) {
        if (fsm->reacquiring_speed && signed_dy > *max_dy_mag) {
            *max_dy_mag = signed_dy;
        }
        float half_v = 0.5f * (*max_dy_mag);
        if (half_v > 0.0f && signed_dy < half_v) {
            fsm->decel_acc = half_v;
            fsm->geo_decel_init = half_v;
            fsm->geo_phi_anchor = afo->phase_out;
            fsm->reacquiring_speed = false;
        }
    } else if (signed_dy > 0.60f * (*max_dy_mag)) {
        /* 双支撑期的暂时减速结束，重新追踪本半步的最大外摆速度。 */
        fsm->decel_acc = fsm->geo_decel_init = 0.0f;
        fsm->reacquiring_speed = true;
        *max_dy_mag = fmaxf(*max_dy_mag, signed_dy);
    } else if (signed_prev_dy > signed_dy) {
        fsm->decel_acc -= signed_prev_dy - signed_dy;
        if (fsm->decel_acc == 0.0f) fsm->decel_acc = -1.0e-4f;
        if (fsm->decel_acc < 0.0f) return true;
    }

    /* 舒缓慢步允许在真实转向且摆幅足够时确认峰。 */
    if (signed_dy <= 0.0f && *peak_mag > (thr_mag + margin_mag)) {
        return true;
    }
    return false;
}

static bool fsm_update_landmarks(hip_gentle_fsm_t *fsm,
                                 const hip_gentle_afo_t *afo,
                                 float y_hip_deg,
                                 float left_deg_s,
                                 float right_deg_s,
                                 double now)
{
    float thr_pos = HIP_GENTLE_START_THRESHOLD_DEG;
    float thr_neg = -HIP_GENTLE_START_THRESHOLD_DEG;
    float margin = HIP_GENTLE_PEAK_MARGIN_DEG;
    float rev_vel = HIP_GENTLE_RAPID_REVERSAL_DEG_S;

    fsm->prev_evt = fsm->evt;

    switch (fsm->evt) {
    case 8: /* 未初始化 / 静止复位 */
        if (y_hip_deg > thr_pos && fsm->prev_y <= thr_pos) {
            fsm->evt = 0;
            fsm->rapid_channel = HIP_GENTLE_RAPID_OFF;


            fsm->max_dy_pos = fsm->dy_ma;
        } else if (y_hip_deg < thr_neg && fsm->prev_y >= thr_neg) {
            fsm->evt = 4;
            fsm->rapid_channel = HIP_GENTLE_RAPID_OFF;


            fsm->min_dy_neg = fsm->dy_ma;
        }
        break;

    case 0: /* 过正阈值起步 [0, +thr] */
        if (y_hip_deg > thr_pos) {
            fsm->seed_peak_pos = fmaxf(thr_pos, y_hip_deg);
            fsm->decel_acc = 0.0f;
            fsm->geo_decel_init = 0.0f;
            fsm->reacquiring_speed = false;
            fsm->evt = 1;


        } else {
            fsm->max_dy_pos = fmaxf(fsm->max_dy_pos, fsm->dy_ma);
        }
        break;

    case 1: /* 正半周上升沿追踪 */
        if (fsm_track_peak_budget(fsm, afo, y_hip_deg, fsm->dy_ma, fsm->prev_dy_ma, right_deg_s,
                                  thr_pos, margin, rev_vel, 1,
                                  &fsm->max_dy_pos, &fsm->seed_peak_pos)) {
            fsm->evt = 2;


            fsm->peak_pos_seen = true;
        }
        break;

    case 2: /* 正波峰确认，继续更新最高值直至低于 +thr */
        if (y_hip_deg < thr_pos) {
            fsm->evt = 3;
            fsm->rapid_channel = HIP_GENTLE_RAPID_OFF;


        } else {
            fsm->seed_peak_pos = fmaxf(fsm->seed_peak_pos, y_hip_deg);
        }
        break;

    case 3: /* 正峰后回落过阈区间 [+thr, 0] */
        if (y_hip_deg < 0.0f) {
            fsm->min_dy_neg = fsm->dy_ma;
            fsm->evt = 4;
            fsm->rapid_channel = HIP_GENTLE_RAPID_OFF;


        }
        break;

    case 4: /* 过负阈值起步 [0, -thr] */
        if (y_hip_deg < thr_neg) {
            fsm->seed_peak_neg = fminf(thr_neg, y_hip_deg);
            fsm->decel_acc = 0.0f;
            fsm->geo_decel_init = 0.0f;
            fsm->reacquiring_speed = false;
            fsm->evt = 5;


        } else {
            fsm->min_dy_neg = fminf(fsm->min_dy_neg, fsm->dy_ma);
        }
        break;

    case 5: { /* 负半周下降沿追踪（正半周镜像） */
        float peak_mag = -fsm->seed_peak_neg;
        float max_dy_mag = -fsm->min_dy_neg;
        bool reached = fsm_track_peak_budget(fsm, afo, -y_hip_deg, -fsm->dy_ma, -fsm->prev_dy_ma,
                                             left_deg_s, thr_pos, margin, rev_vel, -1,
                                             &max_dy_mag, &peak_mag);
        fsm->seed_peak_neg = -peak_mag;
        fsm->min_dy_neg = -max_dy_mag;
        if (reached) {
            fsm->evt = 6;


            fsm->peak_neg_seen = true;
        }
        break;
    }

    case 6: /* 负波谷确认，继续更新最低值直至回升高于 -thr */
        if (y_hip_deg > thr_neg) {
            fsm->evt = 7;
            fsm->rapid_channel = HIP_GENTLE_RAPID_OFF;


        } else {
            fsm->seed_peak_neg = fminf(fsm->seed_peak_neg, y_hip_deg);
        }
        break;

    case 7: /* 负谷后回升过阈区间 [-thr, 0] */
        if (y_hip_deg > 0.0f) {
            fsm->max_dy_pos = fsm->dy_ma;
            fsm->evt = 0;
            fsm->rapid_channel = HIP_GENTLE_RAPID_OFF;


        }
        break;

    default:
        fsm->evt = 8;
        break;
    }

    if (fsm->evt == fsm->prev_evt || fsm->prev_evt == 8 || fsm->evt == 8) return false;
    fsm->evt_ts[fsm->evt] = now;
    fsm->evt_ts_valid[fsm->evt] = true;
    return true;
}

static float fsm_calc_half_period(float dt_quarter, float thr_mag, float peak_mag)
{
    float safe_peak = fmaxf(peak_mag, thr_mag + 1e-4f);
    float denom = fmaxf(HIP_HALF_PI - asin_safe(thr_mag / safe_peak), 0.01f);
    float dt_clamped = fmaxf(dt_quarter, 0.02f);
    return clampf((dt_clamped * HIP_PI) / denom,
                  HIP_GENTLE_MIN_HALF_PERIOD_S,
                  HIP_GENTLE_MAX_HALF_PERIOD_S);
}

/* 已观测或镜像补填的两个特征必须按时间先后排列；首步之前的合成时间可为负。 */
static float landmark_span(const hip_gentle_fsm_t *f, unsigned from, unsigned to)
{
    from &= 7U;
    to &= 7U;
    if (!f->evt_ts_valid[from] || !f->evt_ts_valid[to] || f->evt_ts[to] <= f->evt_ts[from])
        return 0.0f;
    return (float)(f->evt_ts[to] - f->evt_ts[from]);
}

static void backfill_first_half(hip_gentle_fsm_t *f, unsigned base)
{
    float edge = landmark_span(f, base+3, base+4);
    float top = landmark_span(f, base+1, base+3);
    if (edge <= 0 || top <= 0) return;
    double stamp = f->evt_ts[(base+1) & 7U];
    for (unsigned i=0; i<4; ++i) {
        unsigned slot = (base+8-i) & 7U;
        stamp -= i<2 ? edge : 0.5*top;
        f->evt_ts[slot] = stamp;
        f->evt_ts_valid[slot] = true;
    }
}

static float bounded_anchor(float omega, float half_period)
{
    float measured = clampf(HIP_PI / half_period,
                            HIP_GENTLE_AFO_ANCHOR_MIN_RATIO * omega,
                            HIP_GENTLE_AFO_ANCHOR_MAX_RATIO * omega);
    return clampf(0.5f*(omega + measured), HIP_GENTLE_AFO_MIN_OMEGA_RAD_S, HIP_GENTLE_AFO_MAX_OMEGA_RAD_S);
}

static void recover_timeout(hip_gentle_fsm_t *f, float dt)
{
    f->scale_timeout += dt/(HIP_GENTLE_TIMEOUT_RECOVER_TAU_S+dt) * (1.0f-f->scale_timeout);
}

static bool interval_expired(hip_gentle_fsm_t *f, float elapsed, float limit, float dt)
{
    if (elapsed > limit) {
        f->scale_timeout = 1.0f;
        return true;
    }
    if (elapsed > (2.0f/3.0f)*limit) {
        float factor = 1.0f - 2.0f/(1.0f+(1000.0f/3.0f)*limit);
        f->scale_timeout *= powf(fmaxf(0.0f, factor), dt / 0.002f);
    } else {
        recover_timeout(f, dt);
    }
    return false;
}

static bool peak_expired(hip_gentle_fsm_t *f, float elapsed, float dt)
{
    if (elapsed >= f->peak_fade_s + HIP_GENTLE_PEAK_FADE_S) return true;
    if (elapsed > f->peak_fade_s)
        f->scale_timeout = clampf(1.0f-(elapsed-f->peak_fade_s)/HIP_GENTLE_PEAK_FADE_S, 0.0f, 1.0f);
    else
        recover_timeout(f, dt);
    return false;
}

static float amplitude_confidence(float amplitude_rad, float min_rad, float full_rad)
{
    float u = clampf((amplitude_rad - min_rad) / fmaxf(full_rad - min_rad, 1e-4f), 0.0f, 1.0f);
    return 0.5f * (1.0f - cosf(HIP_PI * u));
}

/* 速度仅来自双髋角差；锁相因子、幅值置信度和退缩包络各自有独立时间常数。 */
static void update_supervisors(hip_gentle_assist_t *s, float y_deg, float dt)
{
    hip_gentle_fsm_t *f = &s->fsm;
    float target = 1.0f, tau = HIP_GENTLE_ENV_RECOVER_TAU_S;
    float reverse = 0.0f, reference = 0.0f;
    if (f->evt == 7 || f->evt == 0) {
        f->rise_ref_deg_s = fmaxf(f->rise_ref_deg_s, f->dy_ma);
        reverse = -f->dy_ma;
        reference = f->rise_ref_deg_s;
    } else if (f->evt == 3 || f->evt == 4) {
        f->fall_ref_deg_s = fminf(f->fall_ref_deg_s, f->dy_ma);
        reverse = f->dy_ma;
        reference = -f->fall_ref_deg_s;
    }
    if (fabsf(y_deg) <= HIP_GENTLE_START_THRESHOLD_DEG && reference > 1e-5f && reverse > reference) {
        float ratio = reverse/reference;
        target = ratio >= 2.0f ? 0.0f : sinf(HIP_HALF_PI*ratio);
        tau = HIP_GENTLE_ENV_FAST_TAU_S;
    }
    f->scale_env += dt/(tau+dt) * (target-f->scale_env);
    f->lockin_gain += dt/(HIP_GENTLE_ENV_RECOVER_TAU_S+dt) * (1.0f-f->lockin_gain);
    float amp_afo = fmaxf(0.0f, s->afo.alpha[0]);
    float abs_y_rad = fabsf(y_deg) * HIP_DEG_TO_RAD;
    if (abs_y_rad > f->recent_peak_rad)
        f->recent_peak_rad = abs_y_rad;
    float eff_amp_rad = amp_afo;
    if (f->recent_peak_rad > 0.0f) {
        eff_amp_rad = (f->state == HIP_GENTLE_LOCK_IN)
            ? fmaxf(amp_afo, f->recent_peak_rad)
            : fminf(amp_afo, f->recent_peak_rad);
    }
    float confidence = amplitude_confidence(eff_amp_rad,
                                            s->parameters.amp_scale_min_rad,
                                            s->parameters.amp_scale_full_rad);
    float amp_tau = confidence > f->amp_confidence ? HIP_GENTLE_AMP_CONF_RISE_TAU_S : HIP_GENTLE_AMP_CONF_FALL_TAU_S;
    f->amp_confidence += dt / (amp_tau + dt) * (confidence - f->amp_confidence);

    /* 相平面瞬时步幅能量比：正常行走与叉腿停顿保持 ~1.0，停步减速时在最后半步即平滑减幅 */
    if (f->state == HIP_GENTLE_STEADY) {
        float w_eff = fmaxf(HIP_GENTLE_AFO_MIN_OMEGA_RAD_S,
                            s->afo.omega_filt > 1.0f ? s->afo.omega_filt : s->afo.omega);
        float v_term = (f->dy_ma * HIP_DEG_TO_RAD) / w_eff;
        float r_inst = sqrtf(abs_y_rad * abs_y_rad + v_term * v_term);
        float r_tau = r_inst > f->phase_amp_filt_rad
            ? HIP_GENTLE_ENERGY_AMP_RISE_TAU_S : HIP_GENTLE_ENERGY_AMP_FALL_TAU_S;
        f->phase_amp_filt_rad += dt / (r_tau + dt) * (r_inst - f->phase_amp_filt_rad);
        unsigned cur_side = f->evt < 4 ? 0U : 1U;
        float side_ref_rad = fmaxf(f->reference_peak_deg[cur_side], f->last_peak_deg[cur_side]) * HIP_DEG_TO_RAD;
        float ref_r = fmaxf(fmaxf(amp_afo, side_ref_rad), HIP_GENTLE_ENERGY_MIN_REF_DEG * HIP_DEG_TO_RAD);
        float rho = f->phase_amp_filt_rad / ref_r;
        float u_rho = clampf((rho - HIP_GENTLE_ENERGY_RATIO_LO) /
                             fmaxf(HIP_GENTLE_ENERGY_RATIO_HI - HIP_GENTLE_ENERGY_RATIO_LO, 1e-4f),
                             0.0f, 1.0f);
        float tgt_energy = HIP_GENTLE_ENERGY_SCALE_FLOOR +
                           (1.0f - HIP_GENTLE_ENERGY_SCALE_FLOOR) * (u_rho * u_rho * (3.0f - 2.0f * u_rho));
        float e_tau = tgt_energy > f->scale_energy
            ? HIP_GENTLE_ENERGY_SCALE_RISE_TAU_S : HIP_GENTLE_ENERGY_SCALE_FALL_TAU_S;
        f->scale_energy += dt / (e_tau + dt) * (tgt_energy - f->scale_energy);
    } else {
        f->phase_amp_filt_rad = fmaxf(f->recent_peak_rad, abs_y_rad);
        f->scale_energy = 1.0f;
    }
}

static void fsm_enter_lock_in(hip_gentle_assist_t *s, double now, float quarter_s,
                              float peak_deg, float y_deg, int sign)
{
    hip_gentle_fsm_t *f = &s->fsm;
    f->half_period_s = fsm_calc_half_period(quarter_s, HIP_GENTLE_START_THRESHOLD_DEG, peak_deg);
    f->state = HIP_GENTLE_LOCK_IN;
    f->state_enter_ts = now;
    f->reason = HIP_GENTLE_EVENT_SEED;
    f->phase_override_valid = true;
    f->fsm_phase = sign > 0 ? HIP_HALF_PI : HIP_THREE_HALF_PI;
    f->dir_sign = sign * 0.5f;
    f->lockin_gain = 0.5f;
    f->recent_peak_rad = fmaxf(0.0f, peak_deg) * HIP_DEG_TO_RAD;
    f->amp_confidence = amplitude_confidence(f->recent_peak_rad,
                                             s->parameters.amp_scale_min_rad,
                                             s->parameters.amp_scale_full_rad);
    f->phase_amp_filt_rad = f->recent_peak_rad;
    f->abrupt_shrink_resonance = false;
    f->collapse_halves = 0U;
    f->scale_timeout = f->scale_env = f->scale_energy = 1.0f;
    float mode_scale = s->parameters.extended_timing ? 1.5f : 1.0f;
    f->peak_fade_s = (2.0f/3.0f)*clampf(f->half_period_s*mode_scale, 0.5f, 1000.0f);
    afo_hot_start(s, HIP_PI/f->half_period_s, peak_deg*HIP_DEG_TO_RAD, f->fsm_phase, y_deg*HIP_DEG_TO_RAD);
    f->phi_ref = s->afo.phase_out;
    ++f->step_count;
}

static void fsm_handle_lock_in(hip_gentle_assist_t *s, float y_deg, double now, float dt)
{
    hip_gentle_fsm_t *f = &s->fsm;
    unsigned base = f->dir_sign > 0 ? 0U : 4U;
    int sign = base == 0 ? 1 : -1;
    int crest = (int)((base+2)&7U), release = (int)((base+3)&7U), zero = (int)((base+4)&7U);
    float peak_phase = sign > 0 ? HIP_HALF_PI : HIP_THREE_HALF_PI;
    float peak_deg = sign > 0 ? f->seed_peak_pos : -f->seed_peak_neg;
    float relative_y = sign*y_deg;
    float angle_phase = wrap_2pi(peak_phase + HIP_HALF_PI -
                           asin_safe(clampf(relative_y/fmaxf(peak_deg, 1e-3f), 0.0f, 1.0f)));
    if (fsm_both_static_for(f, now, HIP_GENTLE_ACTIVE_STATIC_S) ||
        fsm_narrow_stopped(f, &s->afo, y_deg, now)) {
        fsm_reset_to_idle(f, &s->afo, now);
        f->reason = HIP_GENTLE_QUIET;
        return;
    }
    if (f->evt == crest) {
        float elapsed = (float)(now - f->evt_ts[crest]);
        if (peak_expired(f, elapsed, dt)) {
            fsm_reset_to_idle(f, &s->afo, now);
            f->reason = HIP_GENTLE_NO_PROGRESS;
            return;
        }
        afo_hot_start(s, HIP_PI/f->half_period_s, peak_deg*HIP_DEG_TO_RAD, peak_phase, y_deg*HIP_DEG_TO_RAD);
        float quarter = 0.5f*f->half_period_s;
        f->fsm_phase = quarter > 0.001f
            ? wrap_2pi(peak_phase + HIP_HALF_PI*clampf(elapsed/quarter, 0.0f, 1.0f))
            : angle_phase;
    } else if (f->evt == release) {
        if (f->prev_evt == crest) {
            if (sign > 0) f->fall_ref_deg_s = f->dy_ma;
            else f->rise_ref_deg_s = f->dy_ma;
            float span = landmark_span(f, base+1, base+3);
            if (span > 0) f->half_period_s = fsm_calc_half_period(0.5f*span, HIP_GENTLE_START_THRESHOLD_DEG, peak_deg);
            f->peak_fade_s = 0.0f;
        }
        float limit = 0.5f*f->half_period_s*(s->parameters.extended_timing ? 4.0f : 1.0f);
        if (interval_expired(f, (float)(now-f->evt_ts[release]), fmaxf(HIP_GENTLE_MIN_INTERVAL_LIMIT_S, limit), dt)) {
            fsm_reset_to_idle(f, &s->afo, now);
            f->reason = HIP_GENTLE_NO_PROGRESS;
            return;
        }
        f->fsm_phase = angle_phase;
    } else if (f->evt == zero && f->prev_evt == release) {
        /* 保护所用的真实过零记录与四槽合成历史分开；合成时间不计为真实步幅。 */
        observe_half_period(f, -sign);
        f->half_period_s = clampf(2.0f*(float)(now-f->evt_ts[crest]), HIP_GENTLE_MIN_HALF_PERIOD_S, HIP_GENTLE_MAX_HALF_PERIOD_S);
        float omega = bounded_anchor(s->afo.omega, f->half_period_s);
        afo_hot_start(s, omega, peak_deg*HIP_DEG_TO_RAD, sign>0 ? HIP_PI : 0.0f, y_deg*HIP_DEG_TO_RAD);
        backfill_first_half(f, base);
        f->state = HIP_GENTLE_STEADY;
        f->reason = HIP_GENTLE_HALF_ANCHORED;
        f->state_enter_ts = now;
        f->dir_sign = 0.0f;
        f->phase_override_valid = false;
        f->scale_timeout = f->scale_energy = 1.0f;
        f->phase_amp_filt_rad = fmaxf(0.0f, peak_deg) * HIP_DEG_TO_RAD;
        f->collapse_halves = 0U;
        f->phi_ref = s->afo.phase_out;
    } else {
        fsm_reset_to_idle(f, &s->afo, now);
        f->reason = HIP_GENTLE_EVENT_MISMATCH;
        return;
    }
    update_supervisors(s, y_deg, dt);
}

/* 正半周用上一负半周的时长，负半周的三个槽位整体平移四格。 */
static float opposite_reference(const hip_gentle_fsm_t *f, unsigned base, bool peak_hold)
{
    float top = landmark_span(f, base+5, base+7);
    float edges = landmark_span(f, base+7, base+1);
    if (top <= 0 || edges <= 0) return f->half_period_s;
    float crest_return = landmark_span(f, base+6, base+7);
    if (peak_hold && f->evt_ts_valid[2] && f->evt_ts_valid[6] && crest_return > 0)
        return fmaxf(crest_return, 0.5f*edges);
    return 0.5f*fmaxf(top, edges);
}

static float steady_interval_limit(const hip_gentle_fsm_t *f, bool extended)
{
    unsigned base = f->evt < 4 ? 0U : 4U;
    float window;
    float min_limit;
    switch (f->evt & 3) {
    case 0:
        window = 2.0f*landmark_span(f, base+3, base+5);
        min_limit = HIP_GENTLE_MIN_NEUTRAL_LIMIT_S;
        break;
    case 1:
        window = 2.0f*opposite_reference(f, base, false);
        min_limit = HIP_GENTLE_MIN_INTERVAL_LIMIT_S;
        break;
    default: /* 3 / 7 */
        window = landmark_span(f, base+7, base+1);
        if (f->return_rebounded) window = fmaxf(window, f->reference_half_s);
        window *= (extended ? 4.0f : 7.0f)*0.5f;
        min_limit = HIP_GENTLE_MIN_RETURN_LIMIT_S;
        break;
    }
    if (window <= 0.0f) window = f->half_period_s;
    /* 11-21 慢步的外摆仍在推进时，0.3 倍会先于峰值耗尽。外摆至少给对侧实测的同段时间。 */
    float coeff = (f->evt & 3) == 1 ? fmaxf(1.0f, f->dyn_to_coeff) : f->dyn_to_coeff;
    return fmaxf(min_limit, window*coeff);
}

static void fsm_handle_steady(hip_gentle_assist_t *s, float y_deg, double now, float dt, bool changed)
{
    hip_gentle_fsm_t *f = &s->fsm;
    hip_gentle_afo_t *a = &s->afo;
    float signal = y_deg * HIP_DEG_TO_RAD;
    unsigned base = f->evt < 4 ? 0U : 4U;
    int part = f->evt & 3;
    if (changed) f->return_rebounded = false;
    if (part == 3 && (base == 0 ? y_deg : -y_deg) > HIP_GENTLE_START_THRESHOLD_DEG + HIP_GENTLE_PEAK_MARGIN_DEG)
        f->return_rebounded = true;
    if (changed) {
        if (part == 1) {
            f->early_peak_armed = false;
            ++f->step_count;
        } else if (part == 2) {
            float span = landmark_span(f, base+6, base+2);
            if (span > 0) f->half_period_s = clampf(span, HIP_GENTLE_MIN_HALF_PERIOD_S, HIP_GENTLE_MAX_HALF_PERIOD_S);
            float peak_deg = base == 0U ? f->seed_peak_pos : -f->seed_peak_neg;
            f->recent_peak_rad = fmaxf(0.0f, peak_deg) * HIP_DEG_TO_RAD;
            unsigned ref_side = base == 0U ? 0U : 1U;
            float ref_peak = fmaxf(f->reference_peak_deg[ref_side], f->last_peak_deg[ref_side]);
            float expected_phase = base == 0U ? HIP_HALF_PI : HIP_THREE_HALF_PI;
            float phase_err_deg = wrap_pm_pi(a->phi[0] - expected_phase) * HIP_RAD_TO_DEG;
            bool early_lead = phase_err_deg < -HIP_GENTLE_EARLY_STOP_PHASE_LEAD_DEG ||
                              (peak_deg < HIP_GENTLE_EARLY_STOP_SEVERE_RATIO * ref_peak &&
                               phase_err_deg < -HIP_GENTLE_EARLY_STOP_SEVERE_LEAD_DEG);
            if (s->parameters.early_stop_enable &&
                ref_peak >= HIP_GENTLE_EARLY_STOP_MIN_REF_DEG &&
                peak_deg < HIP_GENTLE_EARLY_STOP_MAX_PEAK_DEG &&
                peak_deg < HIP_GENTLE_EARLY_STOP_PEAK_RATIO * ref_peak &&
                early_lead) {
                f->early_peak_armed = true;
                f->early_peak_sign = base == 0U ? 1 : -1;
                f->early_peak_deg = peak_deg;
                f->early_peak_ref_deg = ref_peak;
            } else {
                f->early_peak_armed = false;
                if (peak_deg > HIP_GENTLE_START_THRESHOLD_DEG)
                    f->last_peak_deg[ref_side] = peak_deg;
            }
            if (s->parameters.sync_peak_amplitude && f->peak_pos_seen && f->peak_neg_seen) {
                float dc = 0.5f*(f->seed_peak_pos + f->seed_peak_neg);
                a->alpha0 = dc*HIP_DEG_TO_RAD;
                a->alpha[0] = (f->seed_peak_pos-dc)*HIP_DEG_TO_RAD;
                afo_refresh_rates(a, &s->afo_gains, signal);
            }
            /* 异步短边会低估慢腿回程，至少保留刚测得峰谷间隔的一半。 */
            float ref = fmaxf(opposite_reference(f, base, true), 0.5f*f->half_period_s);
            if (s->parameters.extended_timing) ref *= 1.5f;
            f->peak_fade_s = fmaxf(HIP_GENTLE_MIN_PEAK_FADE_S,
                                   (2.0f/3.0f)*clampf(ref*f->dyn_to_coeff, 0.5f, 1000.0f));
        } else if (part == 3) {
            if (base == 0) f->fall_ref_deg_s = f->dy_ma;
            else f->rise_ref_deg_s = f->dy_ma;
            float body = landmark_span(f, base, base+3), lead = landmark_span(f, base, base+1);
            if (body > 0 && lead > 0)
                f->half_period_s = clampf(body+lead, HIP_GENTLE_MIN_HALF_PERIOD_S, HIP_GENTLE_MAX_HALF_PERIOD_S);
        } else if (f->prev_evt == ((f->evt+7)&7)) {
            f->early_peak_armed = false;
            float measured = observe_half_period(f, base == 0 ? 1 : -1);
            if (measured > 0) {
                f->half_period_s = clampf(measured, HIP_GENTLE_MIN_HALF_PERIOD_S, HIP_GENTLE_MAX_HALF_PERIOD_S);
                a->omega = bounded_anchor(a->omega, f->half_period_s);
                afo_refresh_rates(a, &s->afo_gains, signal);
                if (f->abrupt_shrink_resonance ||
                    f->short_halves >= HIP_GENTLE_SHRINK_CONFIRM_HALVES ||
                    f->collapse_halves >= HIP_GENTLE_COLLAPSE_CONFIRM_HALVES) {
                    fsm_enter_cooldown(f, a, now);
                    f->reason = HIP_GENTLE_RESONANCE;
                    return;
                }
            }
        }
    }

    if (s->parameters.geometric_phase_pull && part == 1 && f->decel_acc > 0 && f->geo_decel_init > 1e-6f) {
        float progress = clampf(1-f->decel_acc/f->geo_decel_init, 0, 1);
        float target = circular_lerp(f->geo_phi_anchor, base == 0 ? HIP_HALF_PI : HIP_THREE_HALF_PI, progress);
        afo_sync_phase(s, circular_lerp(a->phase_out, target, HIP_GENTLE_GEO_PHASE_PULL_ALPHA), signal);
    }
    if (fabsf(wrap_pm_pi(a->phase_out-f->phi_ref)) > 0.05f) f->phi_ref = a->phase_out;
    update_supervisors(s, y_deg, dt);

    bool early_stop_hit = false;
    if (f->early_peak_armed) {
        float signed_y = (float)f->early_peak_sign * y_deg;
        if (signed_y > HIP_GENTLE_EARLY_STOP_PEAK_RATIO * f->early_peak_ref_deg ||
            signed_y >= HIP_GENTLE_EARLY_STOP_MAX_PEAK_DEG) {
            f->early_peak_armed = false;
        } else if (signed_y < f->early_peak_deg - HIP_GENTLE_EARLY_STOP_FALL_DEG ||
                   fsm_both_still_for(f, now, HIP_GENTLE_EARLY_STOP_STILL_S)) {
            early_stop_hit = true;
        }
    }

    hip_gentle_reason_t cause = HIP_GENTLE_WAITING_CROSSING;
    bool early_stop_cooldown = false;
    f->high_cadence_elapsed = a->omega_filt > HIP_GENTLE_RESONANCE_OMEGA_RAD_S ? f->high_cadence_elapsed+dt : 0;
    f->low_cadence_elapsed = s->parameters.low_cadence_exit && a->omega_filt <= HIP_GENTLE_LOW_CADENCE_OMEGA_RAD_S ? f->low_cadence_elapsed+dt : 0;
    if (fsm_both_static_for(f, now, HIP_GENTLE_ACTIVE_STATIC_S) ||
        fsm_narrow_stopped(f, a, y_deg, now)) cause = HIP_GENTLE_QUIET;
    else if (early_stop_hit) {
        cause = HIP_GENTLE_QUIET;
        early_stop_cooldown = true;
    }
    else if (f->high_cadence_elapsed + 1e-6f >= HIP_GENTLE_RESONANCE_CONFIRM_S) cause = HIP_GENTLE_RESONANCE;
    else if (f->low_cadence_elapsed + 1e-6f >= HIP_GENTLE_LOW_CADENCE_CONFIRM_S) cause = HIP_GENTLE_LOW_CADENCE;
    else if (f->evt_ts_valid[f->evt]) {
        float elapsed = (float)(now-f->evt_ts[f->evt]);
        bool expired = part == 2 ? peak_expired(f, elapsed, dt)
                      : interval_expired(f, elapsed, steady_interval_limit(f, s->parameters.extended_timing), dt);
        if (expired) cause = HIP_GENTLE_NO_PROGRESS;
    }
    if (cause != HIP_GENTLE_WAITING_CROSSING) {
        fsm_enter_cooldown(f, a, now);
        f->reason = cause;
        if (early_stop_cooldown)
            f->cooldown_min_s = HIP_GENTLE_EARLY_STOP_COOLDOWN_S;
    }
}

static void fsm_step(hip_gentle_assist_t *ctrl,
                     float y_hip_deg,
                     float left_deg_s,
                     float right_deg_s,
                     double now,
                     float dt)
{
    hip_gentle_fsm_t *fsm = &ctrl->fsm;
    hip_gentle_afo_t *afo = &ctrl->afo;

    bool evt_changed = fsm->state != HIP_GENTLE_COOLDOWN &&
                       fsm_update_landmarks(fsm, afo, y_hip_deg, left_deg_s, right_deg_s, now);
    fsm_update_static_and_dynto(fsm, left_deg_s, right_deg_s, now, dt);
    if (fsm->evt != 2 && fsm->evt != 6)
        fsm->peak_fade_s = 0.0f;

    switch (fsm->state) {
    case HIP_GENTLE_IDLE:
        if (afo->enabled) afo_reset(afo);
        fsm->phase_override_valid = false;
        fsm->dir_sign = 0.0f;
        if (ctrl->enabled && evt_changed && (fsm->evt == 1 || fsm->evt == 5)) {
            unsigned side = fsm->evt == 5;
            if (fsm->observed_zero_valid[side]) {
                fsm->half_start_ts = fsm->observed_zero_ts[side];
                fsm->half_start_sign = side ? -1 : 1;
                fsm->observed_zero_valid[side] = false;
            }
            fsm->state = HIP_GENTLE_FIRST_STEP;
            fsm->reason = HIP_GENTLE_WAITING_FIRST_PEAK;
            fsm->rise_ref_deg_s = fabsf(fsm->dy_ma);
            fsm->fall_ref_deg_s = -fsm->rise_ref_deg_s;
            fsm->state_enter_ts = now;
            fsm->low_cadence_elapsed = 0.0f;
            fsm->scale_env = 1.0f;
            fsm->scale_timeout = 1.0f;
            fsm->peak_pos_seen = false;
            fsm->peak_neg_seen = false;
        }
        break;

    case HIP_GENTLE_FIRST_STEP:
        {
            bool crest = evt_changed && (fsm->evt == 2 || fsm->evt == 6);
            bool never_moved = fsm->static_l_since >= 0 && fsm->static_r_since >= 0 &&
                fmax(fsm->static_l_since, fsm->static_r_since) < fsm->state_enter_ts;
            bool quiet = fsm_both_static_for(fsm, now, HIP_GENTLE_FIRST_STEP_STATIC_S) || (crest && never_moved);
            if (quiet || interval_expired(fsm, (float)(now-fsm->state_enter_ts), HIP_GENTLE_FIRST_STEP_TIMEOUT_S, dt)) {
                fsm_reset_to_idle(fsm, afo, now);
                fsm->reason = quiet ? HIP_GENTLE_QUIET : HIP_GENTLE_NO_PROGRESS;
            } else if (crest) {
                int sign = fsm->evt == 2 ? 1 : -1;
                float peak = sign > 0 ? fsm->seed_peak_pos : -fsm->seed_peak_neg;
                fsm_enter_lock_in(ctrl, now, landmark_span(fsm, fsm->evt-1, fsm->evt), peak, y_hip_deg, sign);
            }
        }
        break;

    case HIP_GENTLE_LOCK_IN:
        fsm_handle_lock_in(ctrl, y_hip_deg, now, dt);
        break;

    case HIP_GENTLE_STEADY:
        fsm_handle_steady(ctrl, y_hip_deg, now, dt, evt_changed);
        break;

    case HIP_GENTLE_COOLDOWN: {
        /* 共振锁定等双腿真正停稳，不能用下一次弹性过零重新起步；普通退出也须等卸力归零再回 IDLE。 */
        bool at_rest = fsm_both_static_for(fsm, now, HIP_GENTLE_COOLDOWN_STATIC_S);
        bool settled_neutral = at_rest && fabsf(y_hip_deg) <= HIP_GENTLE_IDLE_OFFSET_LIMIT_RAD * HIP_RAD_TO_DEG;
        bool torque_released = ctrl->torque.release_rate_nm_s == 0.0f &&
                               fabsf(ctrl->torque.carrier_nm) <= 1e-6f;
        bool hold_done = (float)(now - fsm->state_enter_ts) + 1e-6f >= fsm->cooldown_min_s;
        afo->enabled = false;
        fsm->phase_override_valid = false;
        fsm->dir_sign = 0.0f;
        if (torque_released &&
            (settled_neutral ||
             (fsm->reason != HIP_GENTLE_RESONANCE && hold_done &&
              fabsf(y_hip_deg) < HIP_GENTLE_COOLDOWN_THRESHOLD_DEG)))
            fsm_reset_to_idle(fsm, afo, now);
        break;
    }

    default:
        fsm_reset_to_idle(fsm, afo, now);
        break;
    }

    fsm->prev_y = y_hip_deg;
    fsm->prev_y_valid = true;
}

#define PI 3.14159265358979323846
#define TAU (2.0*PI)

double hip_gentle_waveform(double phase, double c, double d)
{
    if (!isfinite(phase) || !isfinite(c) || !isfinite(d) ||
        c < 0.1 || c > 2.0 || d <= -1.0 || d >= 1.0) return NAN;
    const double cc = cos(PI/c);
    const double arg = c < 1.0
        ? (d >= 0.0 ? 2*d/(1+d)-1 : -2*d/(1-d)-1)
        : (d >= 0.0 ? (cc+d)/(1+d) : (cc-d)/(1-d));
    const double arc = asin(fmax(-1.0, fmin(1.0, arg)));
    const double delta = d >= 0.0 ? c/2*(PI/2-arc)-PI/2 : c/2*(arc-PI/2)+PI/2;
    double x = fmod(phase-PI/2+delta, TAU);
    if (x < 0) x += TAU;
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

static void update_torque(hip_gentle_assist_t *s, float y_deg, float dt)
{
    hip_gentle_torque_t *t = &s->torque;
    const hip_gentle_parameters_t *p = &s->parameters;
    const hip_gentle_fsm_t *fsm = &s->fsm;
    ramp_to(&t->gear_ramp_nm, p->gear_limit_nm, HIP_GENTLE_GEAR_RISE_NM_S * dt);
    float z = p->zero_cross_latch ? y_deg : 0.0f;
    if (t->prev_zero_sig * z <= 0.0f)
        t->amplitude_nm = t->gear_ramp_nm;
    t->prev_zero_sig = z;
    ramp_to(&t->resp_ramp, p->response, HIP_GENTLE_RESP_RAMP_RATE_S * dt);
    bool active = s->enabled && (fsm->state == HIP_GENTLE_LOCK_IN || fsm->state == HIP_GENTLE_STEADY);
    if (!active && t->ramp_inner > 0.0f)
        t->release_rate_nm_s = fminf(fabsf(t->carrier_nm) / HIP_GENTLE_EXIT_RAMP_S, HIP_GENTLE_TORQUE_SLEW_NM_S);
    if (!active || t->release_rate_nm_s > 0.0f) {
        t->ramp_inner = t->ramp_outer = t->target_nm = 0.0f;
        t->amp_scale_valid = false;
        t->prev_wave_shape = 0.0f;
        ramp_to(&t->carrier_nm, 0.0f, t->release_rate_nm_s * dt);
        if (fabsf(t->carrier_nm) <= 1e-6f)
            t->carrier_nm = t->release_rate_nm_s = 0.0f;
        return;
    }
    float phase = fsm->phase_override_valid ? fsm->fsm_phase : s->afo.phase_out;
    float offset = clampf(t->resp_ramp / 100.0f - 0.92f, -3.2f, 1.6f);
    float shape = (float)hip_gentle_waveform(phase + offset, p->width, p->asymmetry);
    float raw = t->amplitude_nm * shape;
    ramp_to(&t->ramp_inner, 1.0f, dt / HIP_GENTLE_INNER_RAMP_S);
    ramp_to(&t->ramp_outer, 1.0f, dt / HIP_GENTLE_OUTER_RAMP_S);
    float inner = raw * t->ramp_inner;
    if (fsm->phase_override_valid && inner * fsm->dir_sign > 0.0f)
        inner = 0.0f;
    float carrier = inner * t->ramp_outer;
    if (!t->amp_scale_valid || !p->amp_zero_cross_latch || t->prev_wave_shape * shape <= 0.0f) {
        t->amp_scale_latched = fsm->amp_confidence;
        t->amp_scale_valid = true;
    }
    t->prev_wave_shape = shape;
    if (p->supervision_scale) {
        float amp_scale = p->amp_zero_cross_latch ? t->amp_scale_latched : fsm->amp_confidence;
        carrier *= fsm->scale_timeout * amp_scale * fsm->scale_energy; // 本机保留小幅、相平面能量塌缩与超时保护
        if (p->reversal_envelope) carrier *= fsm->scale_env;
    }
    t->target_nm = clampf(carrier, -p->gear_limit_nm, p->gear_limit_nm);
    ramp_to(&t->carrier_nm, t->target_nm, HIP_GENTLE_TORQUE_SLEW_NM_S * dt);
}

void hip_gentle_default_parameters(hip_gentle_parameters_t *p)
{
    if (!p) return;
    *p = (hip_gentle_parameters_t){
        .gear_limit_nm = 0.0f,
        .response = HIP_GENTLE_DEFAULT_RESPONSE,
        .width = HIP_GENTLE_DEFAULT_WIDTH,
        .asymmetry = HIP_GENTLE_DEFAULT_ASYMMETRY,
        .zero_cross_latch = HIP_GENTLE_DEFAULT_ZERO_CROSS_LATCH,
        .supervision_scale = HIP_GENTLE_DEFAULT_SUPERVISION_SCALE,
        .sync_peak_amplitude = HIP_GENTLE_DEFAULT_SYNC_PEAK_AMP,
        .low_cadence_exit = HIP_GENTLE_DEFAULT_LOW_CADENCE_EXIT,
        .afo_phase_offset_rad = HIP_GENTLE_DEFAULT_PHASE_OFFSET_RAD,
        .extended_timing = HIP_GENTLE_DEFAULT_EXTENDED_TIMING,
        .geometric_phase_pull = HIP_GENTLE_DEFAULT_GEO_PHASE_PULL,
        .reversal_envelope = HIP_GENTLE_DEFAULT_REVERSAL_ENV,
        .amp_scale_min_rad = HIP_GENTLE_AMP_SCALE_MIN_DEG * HIP_DEG_TO_RAD,
        .amp_scale_full_rad = HIP_GENTLE_AMP_SCALE_FULL_DEG * HIP_DEG_TO_RAD,
        .amp_zero_cross_latch = HIP_GENTLE_AMP_ZERO_CROSS_LATCH,
        .early_stop_enable = HIP_GENTLE_EARLY_STOP_ENABLE,
    };
}

void hip_gentle_reset(hip_gentle_assist_t *s)
{
    if (!s || !s->initialized) return;
    hip_gentle_parameters_t p = s->parameters;
    hip_gentle_afo_gains_t gains = s->afo_gains;
    bool enabled = s->enabled, hip_enabled = s->hip_enabled;
    memset(s, 0, sizeof(*s));
    s->parameters = p;
    s->afo_gains = gains;
    s->afo.phase_offset_rad = wrap_2pi(p.afo_phase_offset_rad);
    s->enabled = enabled;
    s->hip_enabled = hip_enabled;
    s->initialized = true;
    fsm_reset_to_idle(&s->fsm, &s->afo, 0.0);
    s->torque.gear_ramp_nm = s->torque.amplitude_nm = p.gear_limit_nm;
    s->torque.resp_ramp = HIP_GENTLE_RESP_RAMP_INIT;
    s->fsm.reason = HIP_GENTLE_WAITING_CROSSING;
}

bool hip_gentle_init(hip_gentle_assist_t *s, const hip_gentle_parameters_t *params)
{
    if (!s) return false;
    hip_gentle_parameters_t p;
    if (params) p = *params;
    else hip_gentle_default_parameters(&p);
    memset(s, 0, sizeof(*s));
    if (!(isfinite(p.gear_limit_nm) && p.gear_limit_nm >= 0.0f &&
          p.response >= 0.0f && p.response <= 100.0f && p.width >= 0.1f && p.width <= 2.0f &&
          p.asymmetry > -1.0f && p.asymmetry < 1.0f && isfinite(p.afo_phase_offset_rad) &&
          isfinite(p.amp_scale_min_rad) && isfinite(p.amp_scale_full_rad) &&
          p.amp_scale_min_rad >= 0.0f && p.amp_scale_full_rad > p.amp_scale_min_rad)) return false;
    s->parameters = p;
    s->initialized = s->enabled = s->hip_enabled = true;
    hip_gentle_set_afo_gains(s, 1.0f, NULL);
    hip_gentle_reset(s);
    return true;
}

bool hip_gentle_set_afo_gains(hip_gentle_assist_t *s, float global_scale,
                               const float trim_percent[4])
{
    if (!s || !s->initialized || !isfinite(global_scale)) return false;
    static const float base[4] = {
        HIP_GENTLE_AFO_GAIN_ETA_ALPHA, HIP_GENTLE_AFO_GAIN_NU_PHI,
        HIP_GENTLE_AFO_GAIN_GAMMA_OMEGA, HIP_GENTLE_AFO_GAIN_ETA_0
    }; // eta_alpha, nu_phi, gamma_omega, eta_0
    float gains[4];
    for (unsigned i = 0; i < 4; ++i) {
        float trim = trim_percent ? trim_percent[i] : 0.0f;
        gains[i] = base[i] * global_scale * (1.0f + 0.01f * trim);
        if (!isfinite(gains[i])) return false;
    }
    s->afo_gains = (hip_gentle_afo_gains_t){.eta_alpha = gains[0], .nu_phi = gains[1],
                                            .gamma_omega = gains[2], .eta_0 = gains[3]};
    if (s->afo.enabled) afo_refresh_rates(&s->afo, &s->afo_gains, s->fsm.prev_y * HIP_DEG_TO_RAD);
    return true;
}

bool hip_gentle_set_rapid_channel(hip_gentle_assist_t *s, hip_gentle_rapid_channel_t channel)
{
    if (!s || !s->initialized || !s->enabled) return false;
    int evt = s->fsm.evt;
    if (channel != HIP_GENTLE_RAPID_OFF &&
        (s->fsm.state == HIP_GENTLE_COOLDOWN ||
         !((channel == HIP_GENTLE_RAPID_POSITIVE && (evt == 0 || evt == 1)) ||
           (channel == HIP_GENTLE_RAPID_NEGATIVE && (evt == 4 || evt == 5))))) return false;
    s->fsm.rapid_channel = channel;
    return true;
}

void hip_gentle_set_enable(hip_gentle_assist_t *s, bool enabled, bool hip_enabled)
{
    if (!s || !s->initialized) return;
    if (s->enabled != enabled) hip_gentle_reset(s);
    s->enabled = enabled;
    s->hip_enabled = hip_enabled;
}

void hip_gentle_set_response_speed(hip_gentle_assist_t *s, float response)
{
    if (s && s->initialized && isfinite(response))
        s->parameters.response = clampf(response, 0.0f, 100.0f);
}

bool hip_gentle_set_gear_limit_nm(hip_gentle_assist_t *s, float gear_limit_nm)
{
    if (!s || !s->initialized) return false;
    bool valid = isfinite(gear_limit_nm) && gear_limit_nm >= 0.0f;
    if (valid && gear_limit_nm == s->parameters.gear_limit_nm) return true;
    float limit = valid ? gear_limit_nm : 0.0f;
    s->parameters.gear_limit_nm = limit;
    hip_gentle_torque_t *t = &s->torque;

    /* 首帧前直接选档；运行中只允许升档等待斜坡或过零，降档同时收紧已锁存值和卸力余量。 */
    if (!s->have_sample) {
        t->gear_ramp_nm = t->amplitude_nm = limit;
    } else {
        t->gear_ramp_nm = fminf(t->gear_ramp_nm, limit);
        t->amplitude_nm = fminf(t->amplitude_nm, limit);
    }
    t->target_nm = clampf(t->target_nm, -limit, limit);
    t->carrier_nm = clampf(t->carrier_nm, -limit, limit);
    if (limit == 0.0f) t->release_rate_nm_s = 0.0f;
    return valid;
}

static bool reject_sample(hip_gentle_assist_t *s, hip_gentle_reason_t reason, hip_gentle_output_t *o)
{
    hip_gentle_reset(s);
    s->fsm.reason = o->reason = reason;
    return false;
}

static bool update(hip_gentle_assist_t *s, double l, double r, double vl, double vr,
                   bool has_velocity, uint32_t ms, bool valid, const float *frame_limit_nm,
                   hip_gentle_output_t *o)
{
    if (!o) return false;
    memset(o, 0, sizeof(*o));
    if (!s || !s->initialized) return false;
    if (frame_limit_nm && !hip_gentle_set_gear_limit_nm(s, *frame_limit_nm))
        return reject_sample(s, HIP_GENTLE_INVALID_GEAR_LIMIT, o);
    float signal_rad = (float)(l - r);
    float y = (signal_rad - s->idle_offset_rad) * HIP_RAD_TO_DEG;
    if (!valid || !isfinite((float)l) || !isfinite((float)r) || !isfinite(y) ||
        (has_velocity && (!isfinite((float)(vl * HIP_RAD_TO_DEG)) ||
                          !isfinite((float)(vr * HIP_RAD_TO_DEG)))))
        return reject_sample(s, HIP_GENTLE_INVALID_ANGLE, o);
    uint32_t elapsed = ms - s->timestamp_ms;
    if (s->have_sample && (!elapsed || elapsed > HIP_GENTLE_RESET_GAP_MS)) {
        hip_gentle_reset(s);
        s->fsm.reason = HIP_GENTLE_TIMESTAMP_RESET;
    }
    if (!s->have_sample) {
        if (has_velocity) {
            s->left_speed = (float)vl;
            s->right_speed = (float)vr;
        }
        s->fsm.prev_y = (signal_rad - s->idle_offset_rad) * HIP_RAD_TO_DEG;
        s->fsm.prev_y_valid = true;
        s->have_sample = true;
    } else if (s->enabled) {
        float sample_dt = (float)elapsed * 0.001f;
        float dt = clampf(sample_dt, 0.002f, 0.050f);
        float raw_l = has_velocity ? (float)vl : ((float)l - s->previous_left) / sample_dt;
        float raw_r = has_velocity ? (float)vr : ((float)r - s->previous_right) / sample_dt;
        float beta = dt / (dt + HIP_GENTLE_STATIC_VEL_FILTER_S);
        s->left_speed += beta * (raw_l - s->left_speed);
        s->right_speed += beta * (raw_r - s->right_speed);
        /* 抬腿超时也能回到 IDLE。明显分腿的姿态不能充当站立零位。 */
        bool near_upright = fabsf(signal_rad) <= HIP_GENTLE_IDLE_OFFSET_LIMIT_RAD + HIP_GENTLE_START_THRESHOLD_DEG * HIP_DEG_TO_RAD;
        if (s->fsm.state == HIP_GENTLE_IDLE && near_upright && fabsf(raw_l) <= HIP_GENTLE_SETTLED_SPEED_RAD_S &&
            fabsf(raw_r) <= HIP_GENTLE_SETTLED_SPEED_RAD_S && fabsf(s->left_speed) <= HIP_GENTLE_SETTLED_SPEED_RAD_S &&
            fabsf(s->right_speed) <= HIP_GENTLE_SETTLED_SPEED_RAD_S) {
            float old = s->idle_offset_rad;
            s->idle_offset_rad = clampf(old + dt / (HIP_GENTLE_IDLE_OFFSET_TAU_S + dt) * (signal_rad - old),
                                        -HIP_GENTLE_IDLE_OFFSET_LIMIT_RAD, HIP_GENTLE_IDLE_OFFSET_LIMIT_RAD);
            /* 校准量的变化不应被微分成运动。 */
            s->fsm.prev_y -= (s->idle_offset_rad - old) * HIP_RAD_TO_DEG;
        }
        signal_rad -= s->idle_offset_rad;
        y = signal_rad * HIP_RAD_TO_DEG;
        fsm_update_dy_ma(&s->fsm, y, dt);
        if (!isfinite(s->fsm.dy_ma) || !isfinite(s->left_speed) || !isfinite(s->right_speed))
            return reject_sample(s, HIP_GENTLE_INVALID_ANGLE, o);
        s->gait_time_s += dt;
        observe_zero(&s->fsm, y, s->gait_time_s, dt);
        float left_deg_s = s->left_speed * HIP_RAD_TO_DEG;
        float right_deg_s = s->right_speed * HIP_RAD_TO_DEG;
        float afo_dt = dt;
        if (s->afo.enabled) {
            float amp_deg = fabsf(s->afo.alpha[0]) * HIP_RAD_TO_DEG;
            bool afo_near_mid = (fabsf(y) <= HIP_GENTLE_NARROW_STOP_DEG &&
                                 fabsf(y) <= HIP_GENTLE_STOP_MID_AMP_RATIO * amp_deg) ||
                                (s->afo.omega_filt >= HIP_GENTLE_FAST_STOP_OMEGA_RAD_S);
            if (afo_near_mid) {
                float vmax_deg_s = fmaxf(fabsf(raw_l), fabsf(raw_r)) * HIP_RAD_TO_DEG;
                float u_mv = clampf((vmax_deg_s - HIP_GENTLE_AFO_FREEZE_MIN_DEG_S) /
                                    fmaxf(HIP_GENTLE_AFO_FREEZE_FULL_DEG_S - HIP_GENTLE_AFO_FREEZE_MIN_DEG_S, 1e-4f),
                                    0.0f, 1.0f);
                afo_dt = dt * (u_mv * u_mv * (3.0f - 2.0f * u_mv));
            }
        }
        afo_step(&s->afo, &s->afo_gains, signal_rad, afo_dt);
        fsm_step(s, y, left_deg_s, right_deg_s, s->gait_time_s, dt);
        update_torque(s, y, dt);
        if (!isfinite(s->afo.phase_out) || !isfinite(s->afo.y_hat) || !isfinite(s->afo.omega) ||
            !isfinite(s->afo.omega_filt) || !isfinite(s->afo.velocity_filt_rad_s) ||
            !isfinite(s->torque.carrier_nm) || fabsf(s->torque.carrier_nm) > s->parameters.gear_limit_nm)
            return reject_sample(s, HIP_GENTLE_INVALID_TRACKER, o);
    }
    s->timestamp_ms = ms;
    s->previous_left = (float)l;
    s->previous_right = (float)r;
    o->state = s->fsm.state;
    o->reason = s->fsm.reason;
    o->phase_rad = s->fsm.phase_override_valid ? s->fsm.fsm_phase : s->afo.phase_out;
    o->frequency_hz = s->afo.enabled ? s->afo.omega / HIP_TWO_PI : 0.0;
    o->filtered_frequency_hz = s->afo.omega_filt / HIP_TWO_PI;
    o->characteristic_velocity_rad_s = s->afo.velocity_filt_rad_s;
    o->signal_rad = l - r;
    o->reconstructed_rad = s->afo.y_hat;
    o->swing_amplitude_rad = fabsf(s->afo.alpha[0]);
    o->gate = s->torque.ramp_inner * s->torque.ramp_outer;
    o->adaptive_amplitude_nm = s->torque.amplitude_nm;
    o->carrier_nm = s->torque.carrier_nm;
    o->left_torque_nm = s->hip_enabled ? o->carrier_nm : 0.0;
    o->right_torque_nm = -o->left_torque_nm;
    o->target_left_torque_nm = s->hip_enabled ? s->torque.target_nm : 0.0;
    o->target_right_torque_nm = -o->target_left_torque_nm;
    return true;
}

bool hip_gentle_step(hip_gentle_assist_t *s, double l, double r, double vl, double vr,
                     uint32_t ms, bool valid, float gear_limit_nm, hip_gentle_output_t *o)
{
    return update(s, l, r, vl, vr, true, ms, valid, &gear_limit_nm, o);
}

bool hip_gentle_update(hip_gentle_assist_t *s, double l, double r, uint32_t ms, bool valid,
                       hip_gentle_output_t *o)
{
    return update(s, l, r, 0.0, 0.0, false, ms, valid, NULL, o);
}

bool hip_gentle_update_with_velocity(hip_gentle_assist_t *s, double l, double r, double vl,
                                     double vr, uint32_t ms, bool valid, hip_gentle_output_t *o)
{
    return update(s, l, r, vl, vr, true, ms, valid, NULL, o);
}
