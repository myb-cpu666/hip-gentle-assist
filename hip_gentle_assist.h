/*
 * 髋关节舒缓助力模式使用说明
 * 输入左右髋角度，内部跟踪角度差的相位，以公共波形生成双髋力矩（N·m）。
 * 公共波形经限幅、限速后左右反向输出；默认伸展与屈曲等幅。
 * 首腿动作仅用于学习相位，确认另一腿开始抬起后才从零渐入助力。
 * 起步定频允许在进入TRACK后的两个有效完整半周期内纠偏，之后恢复正常失配检查。
 *
 * 力矩设置：大摆幅上限由 ASSIST_GENTLE_ADAPTIVE_MAX_NM 设置，默认8 N·m；
 * 小摆幅取上限的 ASSIST_GENTLE_LOW_AMPLITUDE_PERCENT，默认50%（4 N·m）。
 * 判档使用角度差的AFO第一谐波幅值，不是瞬时角度差或速度；0.60 rad升档、0.50 rad降档，
 * 中间保持原档位，幅值平滑过渡。幅值是波形上限，并非每一步都必须达到的实际力矩。
 *
 * 暂停与恢复：稳定跟踪、启动渐入完成、非弱峰候选且仍处于抬腿姿态时，
 * 零力矩下若下一回落方向未确认，相对低速达30ms或该方向等待达50ms，进入短暂停。
 * 暂停时冻结预测相位，实际回落确认或真实过零后，在零力矩帧对齐相位，再平滑接入。
 * 短暂停保留两腿确认和幅值档位；不增加踱步、下楼等场景识别。
 * 已确认的长HOLD可零力矩不限时等待；落腿恢复后有限重校准。中立站稳等停止条件仍有效。
 * 普通停顿平滑卸力；正式EXIT沿限速归零，故障/无效输入复位清零。
 * 2秒无进展超时不是半步完成期限，也不是已确认长HOLD的等待上限。
 *
 * 本项目通过 assist_algorithm.c 接入，不要另建一套状态重复驱动电机。
 * 参数宏在初始化时复制到状态中；修改后需重新编译、烧录，运行中不直接改内部状态。
 * 电机使能、故障及反馈时效由上层检查，本模块仅接收角度、时间戳和有效性标志。
 *
 * 调用示例（放在使用本模块的 .c 文件中）：
 *
 * #include "hip_gentle_assist.h"
 * #include "esp_timer.h"
 *
 * static hip_gentle_assist_t phase_state;
 * static const hip_gentle_parameters_t phase_parameters = {
 *     .amplitude_nm = ASSIST_GENTLE_AMPLITUDE_NM,
 *     .torque_limit_nm = ASSIST_GENTLE_TORQUE_LIMIT_NM,
 *     .extension_torque_ratio = ASSIST_GENTLE_EXTENSION_TORQUE_RATIO,
 *     .adaptive_enabled = ASSIST_GENTLE_ADAPTIVE_ENABLED,
 *     .adaptive_max_nm = ASSIST_GENTLE_ADAPTIVE_MAX_NM,
 *     .swing_up_threshold_rad = ASSIST_GENTLE_SWING_UP_THRESHOLD_RAD,
 *     .swing_down_threshold_rad = ASSIST_GENTLE_SWING_DOWN_THRESHOLD_RAD,
 *     .amplitude_increase_nm_s = ASSIST_GENTLE_AMPLITUDE_INCREASE_NM_S,
 *     .amplitude_decrease_nm_s = ASSIST_GENTLE_AMPLITUDE_DECREASE_NM_S,
 *     .start_angle_rad = ASSIST_GENTLE_START_ANGLE_RAD,
 *     .startup_lift_angle_rad = ASSIST_GENTLE_STARTUP_LIFT_ANGLE_RAD,
 *     .startup_lift_speed_rad_s = ASSIST_GENTLE_STARTUP_LIFT_SPEED_RAD_S,
 *     .startup_lift_confirm_s = ASSIST_GENTLE_STARTUP_LIFT_CONFIRM_S,
 *     .derivative_window_s = ASSIST_GENTLE_DERIVATIVE_WINDOW_S,
 *     .min_frequency_hz = ASSIST_GENTLE_MIN_FREQUENCY_HZ,
 *     .max_frequency_hz = ASSIST_GENTLE_MAX_FREQUENCY_HZ,
 *     .quiet_speed_rad_s = ASSIST_GENTLE_QUIET_SPEED_RAD_S,
 *     .quiet_duration_s = ASSIST_GENTLE_QUIET_DURATION_S,
 *     .hold_confirm_s = ASSIST_GENTLE_HOLD_CONFIRM_S,
 *     .hold_range_rad = ASSIST_GENTLE_HOLD_RANGE_RAD,
 *     .low_speed_release_confirm_s = ASSIST_GENTLE_LOW_SPEED_RELEASE_CONFIRM_S,
 *     .reversal_speed_rad_s = ASSIST_GENTLE_REVERSAL_SPEED_RAD_S,
 *     .reversal_confirm_s = ASSIST_GENTLE_REVERSAL_CONFIRM_S,
 *     .peak_drop_ratio = ASSIST_GENTLE_PEAK_DROP_RATIO,
 *     .weak_half_cycles = ASSIST_GENTLE_WEAK_HALF_CYCLES,
 *     .weak_peak_stop_enabled = ASSIST_GENTLE_WEAK_PEAK_STOP_ENABLED,
 *     .progress_timeout_s = ASSIST_GENTLE_PROGRESS_TIMEOUT_S,
 *     .tracking_error_stop_enabled = ASSIST_GENTLE_TRACKING_ERROR_STOP_ENABLED,
 *     .tracking_error_rad = ASSIST_GENTLE_TRACKING_ERROR_RAD,
 *     .tracking_error_duration_s = ASSIST_GENTLE_TRACKING_ERROR_DURATION_S,
 *     .phase_error_limit_rad = ASSIST_GENTLE_EVENT_ERROR_LIMIT_RAD,
 *     .phase_correction_rad_s = ASSIST_GENTLE_CORRECTION_RAD_S,
 *     .slew_nm_s = ASSIST_GENTLE_SLEW_NM_S,
 *     .motion_rise_s = ASSIST_GENTLE_MOTION_RISE_S,
 *     .motion_fall_s = ASSIST_GENTLE_MOTION_FALL_S,
 *     .reset_gap_s = ASSIST_GENTLE_RESET_GAP_S,
 *     .width = ASSIST_GENTLE_WIDTH, .asymmetry = ASSIST_GENTLE_ASYMMETRY,
 *     .response = ASSIST_GENTLE_RESPONSE,
 *     .inner_rise_s = ASSIST_GENTLE_INNER_RISE_S,
 *     .outer_rise_s = ASSIST_GENTLE_OUTER_RISE_S,
 *     .amplitude_gain = ASSIST_GENTLE_AMPLITUDE_GAIN,
 *     .phase_gain = ASSIST_GENTLE_PHASE_GAIN,
 *     .frequency_gain = ASSIST_GENTLE_FREQUENCY_GAIN,
 *     .bias_gain = ASSIST_GENTLE_BIAS_GAIN,
 * };
 *
 * // 上电初始化一次，返回 false 时不进入助力控制。
 * bool hip_control_init(void)
 * {
 *     return hip_gentle_init(&phase_state, &phase_parameters);
 * }
 *
 * // 助力开启后每 10ms 调用，包括暂停期间；valid 应包含上层对反馈和电机状态的检查。
 * bool hip_control_step(double left_rad, double right_rad, bool valid,
 *                       hip_gentle_output_t *output)
 * {
 *     uint32_t now_ms = (uint32_t)(esp_timer_get_time() / 1000ULL);
 *     return hip_gentle_update(&phase_state, left_rad, right_rad,
 *                             now_ms, valid, output);
 * }
 *
 * // 上层检查 step 返回值后读取 output.left_torque_nm / right_torque_nm；
 * // 返回 false 时停止助力，不继续发送上一帧力矩。电机发送和使能由上层负责。
 * // 停止助力时调用 hip_gentle_reset(&phase_state)，并由上层处理电机停机；
 * // reset 仅清空算法状态，不负责发送零力矩或执行平滑卸力。
 */
#ifndef HIP_GENTLE_ASSIST_H
#define HIP_GENTLE_ASSIST_H

#include <stdbool.h>
#include <stdint.h>

/* ============================================================================
 * 舒缓模式专用参数：浮点常量及内部数学状态使用 double。
 * 假设编码器左右均以屈曲为正，输入左角度减右角度；公共波形限速后左右反号输出。
 * 持续低速确认后目标归零，短时保留跟踪；确认实际运动方向且指令回零后才允许换向。
 * ============================================================================ */
#define ASSIST_GENTLE_ADAPTIVE_MAX_NM 8.0     // 大摆幅目标幅值/绝对限幅，N·m
#define ASSIST_GENTLE_LOW_AMPLITUDE_PERCENT 50 // 小摆幅占上限的整数百分比，0～100；0关闭舒缓输出
#if ASSIST_GENTLE_LOW_AMPLITUDE_PERCENT < 0 || ASSIST_GENTLE_LOW_AMPLITUDE_PERCENT > 100
#error "Gentle low amplitude percentage must be an integer from 0 to 100"
#endif
#define ASSIST_GENTLE_AMPLITUDE_NM (ASSIST_GENTLE_ADAPTIVE_MAX_NM * ASSIST_GENTLE_LOW_AMPLITUDE_PERCENT / 100.0) // 派生值，不单独设置
#define ASSIST_GENTLE_TORQUE_LIMIT_NM ASSIST_GENTLE_ADAPTIVE_MAX_NM // 绝对输出限幅，与自适应最大值共用
#define ASSIST_GENTLE_EXTENSION_TORQUE_RATIO 1.0 // 伸展相对屈曲的幅值比例，0～1；0关闭伸展，1等幅
#define ASSIST_GENTLE_START_ANGLE_RAD 0.20    // 左右角度差向外越过±0.20rad开始学习，约11.5°
#define ASSIST_GENTLE_STARTUP_LIFT_ANGLE_RAD 0.05 // 单腿连续屈曲段最小抬升量，rad，约2.9°
#define ASSIST_GENTLE_STARTUP_LIFT_SPEED_RAD_S 0.12 // 单腿屈曲确认速度，rad/s
#define ASSIST_GENTLE_STARTUP_LIFT_CONFIRM_S 0.05 // 连续屈曲确认时间，s；确认第二条腿后才渐入出力
#define ASSIST_GENTLE_FEEDBACK_MAX_AGE_MS 100U // 编码器反馈有效期，ms
#define ASSIST_GENTLE_DERIVATIVE_WINDOW_S 0.04 // 约四个 10ms 间隔；固定缓冲支持不超过 61ms
#define ASSIST_GENTLE_MIN_FREQUENCY_HZ 0.15  // 完整步态周期频率下限，Hz；对应半周期约3.33s，仍受进展及失配检查约束
#define ASSIST_GENTLE_MAX_FREQUENCY_HZ 3.50  // 跟踪频率上限，Hz

/* ============================================================================
 * 启停与失锁保护：静止确认、弱半波、无进展退出及跟踪误差检测。
 * ============================================================================ */
#define ASSIST_GENTLE_QUIET_SPEED_RAD_S 0.10 // 角度差速度低速阈值，rad/s，连续确认后卸力；也用于首峰回落确认
#define ASSIST_GENTLE_LOW_SPEED_RELEASE_CONFIRM_S 0.05 // 连续低速撤力确认，s；也用作零力矩下下一回落方向等待的短暂停门槛
#define ASSIST_GENTLE_QUIET_DURATION_S 2.0  // 非悬停低速且无进展的退出时间，也用于站稳解除防震荡门槛，s
#define ASSIST_GENTLE_HOLD_CONFIRM_S 1.0    // 抬腿姿态在允许范围内持续稳定的确认时间，s
#define ASSIST_GENTLE_HOLD_RANGE_RAD 0.08   // 候选窗口内每腿及角度差的最大变化范围，rad，允许姿态微调
#define ASSIST_GENTLE_REVERSAL_SPEED_RAD_S 0.20 // 允许重新出力/换向的角度差速度幅值，rad/s，应高于低速阈值
#define ASSIST_GENTLE_REVERSAL_CONFIRM_S 0.03 // 同方向连续确认时间，s；也用作零力矩下相对低速触发短暂停的门槛
#define ASSIST_GENTLE_PEAK_DROP_RATIO 0.50  // 弱峰相对门槛=最近2个同方向有效半波最大峰值×此比例；调低会减弱停止抑振
#define ASSIST_GENTLE_WEAK_HALF_CYCLES 2U   // 连续弱半波退出数量
#define ASSIST_GENTLE_WEAK_PEAK_STOP_ENABLED 1 // 弱峰回摆先卸力保留步态，过零仍弱才退出；0关闭此项
#if ASSIST_GENTLE_WEAK_PEAK_STOP_ENABLED != 0 && ASSIST_GENTLE_WEAK_PEAK_STOP_ENABLED != 1
#error "Weak peak stop switch must be 0 or 1"
#endif
#define ASSIST_GENTLE_PROGRESS_TIMEOUT_S 2.0 // 连续无有效角度进展的退出时间，s；不是半步完成期限
#define ASSIST_GENTLE_TRACKING_ERROR_STOP_ENABLED 0 // 0 关闭重建误差超时退出，1 开启
#if ASSIST_GENTLE_TRACKING_ERROR_STOP_ENABLED != 0 && ASSIST_GENTLE_TRACKING_ERROR_STOP_ENABLED != 1
#error "Phase tracking error stop switch must be 0 or 1"
#endif
#define ASSIST_GENTLE_TRACKING_ERROR_RAD 0.50 // 重建绝对误差阈值，rad
#define ASSIST_GENTLE_TRACKING_ERROR_DURATION_S 0.30 // 误差连续超限时间，s
#define ASSIST_GENTLE_EVENT_ERROR_LIMIT_RAD 1.5707963267948966 // 事件相位失配阈值，rad
#define ASSIST_GENTLE_CORRECTION_RAD_S 1.0   // 持续跟踪的小相位偏差限速修正
#define ASSIST_GENTLE_SLEW_NM_S 50.0         // 实际 dt 限速，10ms 最多变化 0.5 N·m
#define ASSIST_GENTLE_MOTION_RISE_S 0.16     // 首次/停顿/大目标或方向等待后接入的局部渐入，s；已确认方向的近零换向直接跟随
#define ASSIST_GENTLE_MOTION_FALL_S 0.08     // 普通停顿的平滑卸力时间下限，s；连续相位回零和正式退出沿用限速
#define ASSIST_GENTLE_RESET_GAP_S 0.02       // 实际控制周期超过 20ms 即复位清零

/* ============================================================================
 * 波形、渐入和自适应振荡器：两级渐入并行相乘。
 * ============================================================================ */
#define ASSIST_GENTLE_WIDTH 1.04
#define ASSIST_GENTLE_ASYMMETRY (-0.02)
#define ASSIST_GENTLE_RESPONSE 65.0          // 相位偏移 response/100-0.92 rad
#define ASSIST_GENTLE_INNER_RISE_S 0.6       // 内层渐入时间，s；按实际 dt/0.6 推进
#define ASSIST_GENTLE_OUTER_RISE_S 0.45      // 外层渐入时间，s；两层并行相乘，加快首步接入
#define ASSIST_GENTLE_AMPLITUDE_GAIN 2.5
#define ASSIST_GENTLE_PHASE_GAIN 1.6
#define ASSIST_GENTLE_FREQUENCY_GAIN 3.2
#define ASSIST_GENTLE_BIAS_GAIN 0.5

/* ============================================================================
 * 角摆幅两档调幅：与网页角摆幅同源，取fabs(AFO第一谐波幅值)，单位rad。
 * 小摆幅=上限×百分比，默认4/8 N·m；升降档门槛不同，区间内保持原档位。
 * 首次过渡、悬停恢复再校准、弱峰候选及起步渐入期间回到低档。
 * ============================================================================ */
#define ASSIST_GENTLE_ADAPTIVE_ENABLED 1      // 0 固定基准幅值，1 稳定跟踪后按角摆幅分档
#define ASSIST_GENTLE_SWING_UP_THRESHOLD_RAD 0.60   // 低档达到此角摆幅时升档，rad
#define ASSIST_GENTLE_SWING_DOWN_THRESHOLD_RAD 0.50 // 高档降到此角摆幅时降档，rad；须小于升档门槛
#define ASSIST_GENTLE_AMPLITUDE_INCREASE_NM_S 3.0 // 幅值包络增加速度，N·m/s
#define ASSIST_GENTLE_AMPLITUDE_DECREASE_NM_S 6.0 // 幅值包络降低速度；停止不等待此包络
#if ASSIST_GENTLE_ADAPTIVE_ENABLED != 0 && ASSIST_GENTLE_ADAPTIVE_ENABLED != 1
#error "Adaptive amplitude switch must be 0 or 1"
#endif

/* 舒缓模式独立副本，源自 hip_phase_assist.c/.h；后续舒缓助力在本模块调整，不含电机通信。
 * 固定缓冲、无堆分配；实际运行耗时须实机测量。 */
#define HIP_GENTLE_HISTORY_CAPACITY 64U
#define HIP_GENTLE_HARMONICS 5U
#define HIP_GENTLE_PEAK_HISTORY_CAPACITY 2U // 每个方向各自保存，总计4个半波
#define HIP_GENTLE_STARTUP_SYNC_HALVES 2U // 首次校准观察完整左右周期，避免只验证较快一侧
#define HIP_GENTLE_RESUME_SYNC_HALVES 4U // 悬停恢复后用两个完整周期重新确认频率

typedef enum {
    HIP_GENTLE_WAIT, HIP_GENTLE_COLLECT, HIP_GENTLE_START, HIP_GENTLE_TRACK, HIP_GENTLE_EXIT,
    HIP_GENTLE_HOLD
} hip_gentle_state_id_t;

typedef enum {
    HIP_GENTLE_WAITING_CROSSING, HIP_GENTLE_INVALID_ANGLE, HIP_GENTLE_TIMESTAMP_RESET,
    HIP_GENTLE_WAITING_NEW_CROSSING, HIP_GENTLE_QUIET, HIP_GENTLE_NO_PROGRESS,
    HIP_GENTLE_INCOMPLETE_HALF_WAVE, HIP_GENTLE_INVALID_FIRST_PEAK,
    HIP_GENTLE_SEED_FREQUENCY_RANGE, HIP_GENTLE_EVENT_SEED, HIP_GENTLE_INVALID_TRACKER,
    HIP_GENTLE_AMPLITUDE_DROP, HIP_GENTLE_STARTUP_RESEED, HIP_GENTLE_EVENT_MISMATCH,
    HIP_GENTLE_TRACKING_ERROR, HIP_GENTLE_HISTORY_OVERFLOW, HIP_GENTLE_WEAK_PEAK,
    HIP_GENTLE_HOLDING, HIP_GENTLE_RESUMED
} hip_gentle_reason_t;

typedef struct {
    double amplitude_nm, torque_limit_nm, start_angle_rad, derivative_window_s;
    double startup_lift_angle_rad, startup_lift_speed_rad_s, startup_lift_confirm_s;
    double extension_torque_ratio;
    double min_frequency_hz, max_frequency_hz;
    double quiet_speed_rad_s, quiet_duration_s, peak_drop_ratio;
    double hold_confirm_s, hold_range_rad;
    double low_speed_release_confirm_s;
    double reversal_speed_rad_s, reversal_confirm_s;
    uint32_t weak_half_cycles;
    bool weak_peak_stop_enabled;
    double progress_timeout_s;
    bool tracking_error_stop_enabled;
    double tracking_error_rad, tracking_error_duration_s;
    double phase_error_limit_rad, phase_correction_rad_s;
    double slew_nm_s, reset_gap_s;
    double motion_rise_s, motion_fall_s;
    double width, asymmetry, response, inner_rise_s, outer_rise_s;
    double amplitude_gain, phase_gain, frequency_gain, bias_gain;
    bool adaptive_enabled;
    double adaptive_max_nm, swing_up_threshold_rad, swing_down_threshold_rad;
    double amplitude_increase_nm_s, amplitude_decrease_nm_s;
} hip_gentle_parameters_t;

typedef struct {
    double left_torque_nm, right_torque_nm;
    hip_gentle_state_id_t state;
    hip_gentle_reason_t reason;
    double phase_rad, frequency_hz, signal_rad, reconstructed_rad, gate;
    double swing_amplitude_rad; // 角度差第一谐波幅值绝对值，不是米制步长
    double target_left_torque_nm;  // 限速前目标，正常停止当帧即清零
    double target_right_torque_nm;
    double adaptive_amplitude_nm; // 诊断：当前平滑波形幅值，N·m，非反馈力矩
} hip_gentle_output_t;

typedef struct {
    double time_s, left_rad, right_rad;
} hip_gentle_history_t;

typedef struct {
    hip_gentle_parameters_t parameters;
    bool initialized, has_timestamp, bootstrap_pending;
    bool startup_reference_valid, startup_released, startup_second_confirmed;
    bool collect_wait_cycle, collect_has_zero;
    double collect_zero_time;
    uint8_t startup_first_leg; // 0未确认，1左腿，2右腿；正式退出/复位后重新识别
    double startup_min_angle[2], startup_lift_time[2];
    uint32_t timestamp_ms, history_head, history_count, weak_count;
    hip_gentle_state_id_t state;
    hip_gentle_state_id_t pause_from_state;
    hip_gentle_reason_t reason;
    int direction;
    int torque_direction, motion_direction; // 公共力矩已允许方向、实测角度差运动候选方向
    int completed_direction; // 已结束的相位出力方向；临时卸力不视作结束
    int paused_direction; // 临时卸力不消费当前相位段；同向恢复需有新的实际位移
    double paused_x;
    bool motion_waiting; // 低速卸力后的许可等待，不等于正式退出或暂停整个振荡器
    bool motion_releasing;
    bool motion_ramp_pending;
    int requested_direction;
    double direction_wait_s; // 等实际方向许可的时间，不含正常跨零卸力
    bool phase_paused, pause_used_this_half; // 短停与长HOLD共用暂停状态；每半波限制短停重入
    double motion_entry_elapsed, motion_release_elapsed, motion_release_duration, motion_release_torque;
    double torque_rate_nm_s, motion_release_rate; // 卸力衔接保留已有的回零速度
    double time_s, previous_x, peak, last_peak, cross_time, peak_time, event_time;
    double progress_time, progress_x;
    bool progress_returning;
    bool weak_pending; // 局部弱峰回摆：仅卸力，等待同半波继续增大或确认结束
    double hold_quiet_time, pause_anchor_x;
    double hold_low[3], hold_high[3]; // 左腿、右腿、角度差的候选稳定范围
    bool pause_zeroed;
    bool relearn_peaks_after_hold;
    bool resume_cross_seen;
    uint8_t resume_sync_remaining, bootstrap_halves_remaining;
    double resume_cross_time;
    double recent_peaks[2][HIP_GENTLE_PEAK_HISTORY_CAPACITY];
    uint32_t recent_peak_head[2], recent_peak_count[2];
    double restart_angle_rad, cross_angle_rad; // 弱峰退出后抑制残余摆动；记录本次实际越阈值角度
    double quiet_time, bad_time, enabled_time, phase_correction;
    double motion_confirm_time;
    double target, torque, phase, gate, estimate, bias, omega;
    double adaptive_amplitude_nm;
    bool adaptive_high; // 当前目标档位，与平滑过渡中的实际幅值分开
    double amplitudes[HIP_GENTLE_HARMONICS], phases[HIP_GENTLE_HARMONICS], rates[12];
    hip_gentle_history_t history[HIP_GENTLE_HISTORY_CAPACITY];
} hip_gentle_assist_t;

bool hip_gentle_init(hip_gentle_assist_t *state, const hip_gentle_parameters_t *parameters);
void hip_gentle_reset(hip_gentle_assist_t *state);
/* 100Hz 每帧一次，timestamp_ms 应为实际单调时钟。无效数据返回 false 且清零。
 * 周期重复/倒退/超期则复位并返回有效的零输出，等待新的越阈值事件。
 * 两侧输入须零位正确、同以屈曲为正。左右反号输出为本机适配设计。 */
bool hip_gentle_update(hip_gentle_assist_t *state, double left_rad, double right_rad,
                      uint32_t timestamp_ms, bool valid, hip_gentle_output_t *output);
/* 单位幅值形状；无效参数返回 NAN。 */
double hip_gentle_waveform(double phase, double width, double asymmetry);

#endif
