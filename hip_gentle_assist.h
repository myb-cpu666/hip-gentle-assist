#ifndef HIP_GENTLE_ASSIST_H
#define HIP_GENTLE_ASSIST_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 100 Hz 髋舒缓模式。两腿均以屈曲为正，输入 rad、rad/s，输出解剖系 N·m。
 * 调用 hip_gentle_init(&state, NULL) 初始化一次。之后在控制任务中每 10 ms 调用 hip_gentle_step，
 * 同时传入双髋反馈、本地单调毫秒时间戳和本帧每侧最大力矩（N·m）。
 * FIRST_STEP 观察首峰。LOCK_IN 和 STEADY 生成波形；正常退出按当前力矩线性卸力，硬失能立即清零。
 * 升档幅值默认连续渐增，zero_cross_latch 可选角差过零锁存；降档立即按新上限约束双髋。
 */
/* ==================== 一、采样与零偏校准 ==================== */
#define ASSIST_GENTLE_FEEDBACK_MAX_AGE_MS  100U  /* 编码器反馈最大允许过期时间 (ms) */
#define HIP_GENTLE_RESET_GAP_MS            50U   /* 两帧间隔超过此阈值触发时间戳复位 (ms) */
#define HIP_GENTLE_DY_WINDOW               4U    /* 角差导数滑窗点数（100 Hz 下约 40 ms） */
#define HIP_GENTLE_STATIC_VEL_FILTER_S     0.010f /* 两侧角速度一阶低通时间常数 (s) */
#define HIP_GENTLE_IDLE_OFFSET_LIMIT_RAD   0.15f /* 站立静止零偏自适应校准最大允许偏置 (rad) */
#define HIP_GENTLE_IDLE_OFFSET_TAU_S       0.50f /* 站立静止零偏慢跟踪时间常数 (s) */

/* ==================== 二、默认助力与波形参数 ==================== */
#define HIP_GENTLE_DEFAULT_RESPONSE          0.0f   /* 默认响应超前/滞后调节值 [0, 100]（越大越超前，越小越偏向屈曲托举、伸展柔和） */
#define HIP_GENTLE_DEFAULT_WIDTH             1.04f   /* 升余弦波形占空/脉宽参数 [0.1, 2.0] */
#define HIP_GENTLE_DEFAULT_ASYMMETRY         (0.00f)/* 升余弦波形非对称偏置 (-1, 1) */
#define HIP_GENTLE_DEFAULT_PHASE_OFFSET_RAD  -0.15f    /* AFO 输出相位附加偏置 (rad) */
#define HIP_GENTLE_DEFAULT_ZERO_CROSS_LATCH  false   /* true=角差过零时锁存升档幅值，降档总是立即限幅 */
#define HIP_GENTLE_DEFAULT_SUPERVISION_SCALE true    /* 默认开启小幅与超时监督缩放 */
#define HIP_GENTLE_DEFAULT_SYNC_PEAK_AMP     false   /* 是否在 STEADY 波峰处硬同步 AFO 直流/幅值 */
#define HIP_GENTLE_DEFAULT_LOW_CADENCE_EXIT  false   /* 默认关闭低步频连续退出保护（避免叉腿停顿后落腿续步被误断） */
#define HIP_GENTLE_DEFAULT_EXTENDED_TIMING   false   /* 默认关闭扩展半步超时倍率 */
#define HIP_GENTLE_DEFAULT_GEO_PHASE_PULL    false   /* 默认关闭减速段几何相位柔性牵引 */
#define HIP_GENTLE_DEFAULT_REVERSAL_ENV      false   /* 默认关闭角差退缩包络乘法调制（仅观测） */

/* ==================== 三、软启动、斜坡与卸力限速 ==================== */
#define HIP_GENTLE_GEAR_RISE_NM_S         (20.0f / 3.0f) /* 升档时波形幅值增加速率 (N·m/s) */
#define HIP_GENTLE_RESP_RAMP_INIT          HIP_GENTLE_DEFAULT_RESPONSE /* 复位时响应参数斜坡初始值（与默认响应对齐，避免起步相位漂移） */
#define HIP_GENTLE_RESP_RAMP_RATE_S        20.0f /* 响应参数变化追踪斜坡速率 (1/s) */
#define HIP_GENTLE_INNER_RAMP_S            2.0f  /* 级联内层软启动总时长 (s) */
#define HIP_GENTLE_OUTER_RAMP_S            1.5f  /* 级联外层软启动总时长 (s) */
#define HIP_GENTLE_EXIT_RAMP_S             0.08f /* 正常停步/退步线性卸力参考时间 (s) */
#define HIP_GENTLE_TORQUE_SLEW_NM_S        60.0f /* 最终输出力矩变化率硬限幅 (N·m/s) */

/* ==================== 四、步幅-力矩上限映射（S 曲线） ==================== */
#define HIP_GENTLE_AMP_SCALE_MIN_DEG       6.0f  /* 步幅上限缩放起始角差幅值 (deg，<=此值为 0) */
#define HIP_GENTLE_AMP_SCALE_FULL_DEG      45.0f /* 步幅上限缩放置满角差幅值 (deg，>=此值为 1) */
#define HIP_GENTLE_AMP_CONF_RISE_TAU_S     0.15f /* 步幅增大时上限置信度慢升滤波时间常数 (s) */
#define HIP_GENTLE_AMP_CONF_FALL_TAU_S     0.06f /* 步幅减小时上限置信度快降滤波时间常数 (s) */
#define HIP_GENTLE_AMP_ZERO_CROSS_LATCH    true  /* true=仅在波形过零点锁存步幅上限（单半步保形），false=逐帧实时缩放 */
#define HIP_GENTLE_ENERGY_MIN_REF_DEG      10.0f /* 相平面步幅能量比参考幅值下限 (deg) */
#define HIP_GENTLE_ENERGY_RATIO_HI         0.65f /* 相平面步幅能量比开始平滑衰减门槛（>=此值不衰减） */
#define HIP_GENTLE_ENERGY_RATIO_LO         0.25f /* 相平面步幅能量比达到最低倍率门槛 */
#define HIP_GENTLE_ENERGY_SCALE_FLOOR      0.20f /* 相平面步幅能量塌缩时的最低平滑缩放倍率 */
#define HIP_GENTLE_ENERGY_AMP_RISE_TAU_S   0.08f /* 瞬时相平面幅值上升低通时间常数 (s) */
#define HIP_GENTLE_ENERGY_AMP_FALL_TAU_S   0.04f /* 瞬时相平面幅值下降快跟时间常数 (s) */
#define HIP_GENTLE_ENERGY_SCALE_RISE_TAU_S 0.08f /* 能量缩放倍率恢复时间常数 (s) */
#define HIP_GENTLE_ENERGY_SCALE_FALL_TAU_S 0.04f /* 能量缩放倍率下降时间常数 (s) */


/* ==================== 五、五谐波 AFO 振荡器参数 ==================== */
#define HIP_GENTLE_AFO_SUBSTEPS            5U     /* 单帧控制周期内部积分细分步数 */
#define HIP_GENTLE_AFO_GAIN_ETA_ALPHA      2.5f   /* AFO 基波幅值自适应学习基准增益 */
#define HIP_GENTLE_AFO_GAIN_NU_PHI         6.0f   /* AFO 相位自适应学习基准增益（弧度制） */
#define HIP_GENTLE_AFO_GAIN_GAMMA_OMEGA    5.0f   /* AFO 基频自适应学习基准增益（弧度制） */
#define HIP_GENTLE_AFO_GAIN_ETA_0          0.5f   /* AFO 直流偏置自适应学习基准增益 */
#define HIP_GENTLE_AFO_MIN_OMEGA_RAD_S     2.0f   /* AFO 角频率安全下限 (rad/s，约 38 步/分) */
#define HIP_GENTLE_AFO_MAX_OMEGA_RAD_S     12.5f  /* AFO 角频率安全上限 (rad/s) */
#define HIP_GENTLE_AFO_FILTER_TAU_S        0.126f /* 步频与特征速度低通滤波时间常数 (s) */
#define HIP_GENTLE_AFO_ANCHOR_MIN_RATIO    0.75f  /* 过零重锚允许单步降频最小倍率 */
#define HIP_GENTLE_AFO_ANCHOR_MAX_RATIO    1.30f  /* 过零重锚允许单步升频最大倍率 */
#define HIP_GENTLE_GEO_PHASE_PULL_ALPHA    0.02f  /* 减速段可选几何相位柔性牵引权重 */

/* ==================== 六、八特征点与五状态机时序门限 ==================== */
#define HIP_GENTLE_START_THRESHOLD_DEG     2.5f   /* 离开中立位的角差起步阈值 (deg) */
#define HIP_GENTLE_COOLDOWN_THRESHOLD_DEG  0.5f   /* 普通退出 COOLDOWN 回落回 IDLE 角差阈值 (deg) */
#define HIP_GENTLE_PEAK_MARGIN_DEG         0.5f   /* 慢步转峰判据超阈裕度 (deg) */
#define HIP_GENTLE_RAPID_REVERSAL_DEG_S    500.0f /* 急换向保护通道触发单侧屈曲角速度门限 (deg/s) */
#define HIP_GENTLE_STATIC_SPEED_RAD_S      0.70f  /* 步态静止检测速度死区 (rad/s，约 40.1 deg/s) */
#define HIP_GENTLE_SETTLED_SPEED_RAD_S     0.25f  /* 站立零偏校准与冷却停稳确认速度死区 (rad/s) */
#define HIP_GENTLE_FIRST_STEP_STATIC_S     1.0f   /* FIRST_STEP 双侧静止超时退出时长 (s) */
#define HIP_GENTLE_FIRST_STEP_TIMEOUT_S    3.0f   /* FIRST_STEP 未达峰谷总超时退出时长 (s) */
#define HIP_GENTLE_ACTIVE_STATIC_S         1.5f   /* LOCK_IN/STEADY 叉腿宽站姿双侧静止停步时长 (s) */
#define HIP_GENTLE_COOLDOWN_STATIC_S       0.5f   /* COOLDOWN 共振退出后双侧停稳恢复时长 (s) */
#define HIP_GENTLE_DYN_TO_SPEED_DEG_S      100.0f /* 动态超时系数高速/低速切换活动速度门限 (deg/s) */
#define HIP_GENTLE_DYN_TO_FAST_COEFF       3.0f   /* 正常活动或峰/谷保持时的动态超时倍率 */
#define HIP_GENTLE_DYN_TO_SLOW_COEFF       0.3f   /* 低速停滞收缩时的动态超时收紧倍率 */
#define HIP_GENTLE_DYN_TO_TAU_S            0.10f  /* 动态超时系数 k_to 一阶低通时间常数 (s) */
#define HIP_GENTLE_MIN_INTERVAL_LIMIT_S    0.40f  /* 外摆区间 (evt=1/5) 超时监督最小保底时限 (s) */
#define HIP_GENTLE_MIN_NEUTRAL_LIMIT_S     0.12f  /* 并腿过零区 (evt=0/4) 超时监督最小保底时限 (s) */
#define HIP_GENTLE_MIN_RETURN_LIMIT_S      0.15f  /* 峰后回程区 (evt=3/7) 超时监督最小保底时限 (s) */
#define HIP_GENTLE_STILL_SPEED_RAD_S       0.35f  /* 真静止检测速度死区 (rad/s，约 20.1 deg/s) */
#define HIP_GENTLE_NARROW_STOP_DEG         15.0f  /* 并腿窄站姿快退判定角差上限 (deg) */
#define HIP_GENTLE_NARROW_STOP_STATIC_S    0.22f  /* 并腿窄站姿双侧真静止快退确认时长上限 (s) */
#define HIP_GENTLE_NARROW_STOP_MIN_S       0.06f  /* 并腿中段/快踏步自适应静止快退确认时长下限 (s) */
#define HIP_GENTLE_NARROW_STOP_HALF_RATIO  0.25f  /* 并腿中段/快踏步静止快退时长占当前半周期比例 */
#define HIP_GENTLE_STOP_MID_AMP_RATIO      0.45f  /* 判定步幅中段并腿区相对 AFO 幅值的角差比例上限（避开慢走波峰） */
#define HIP_GENTLE_FAST_STOP_OMEGA_RAD_S   7.5f   /* 高步频快踏步急停保护角频率门槛 (rad/s，约 143 步/分) */
#define HIP_GENTLE_AFO_FREEZE_MIN_DEG_S    2.0f   /* 中段/快踏步 AFO 随腿速冻结下限 (deg/s，<=此值 AFO 立即停转) */
#define HIP_GENTLE_AFO_FREEZE_FULL_DEG_S   12.0f  /* 中段/快踏步 AFO 全速积分上限 (deg/s，>=此值 AFO 正常积分) */
#define HIP_GENTLE_PEAK_FADE_S             0.20f  /* 峰/谷超时后预淡出持续窗口 (s) */
#define HIP_GENTLE_MIN_PEAK_FADE_S         1.30f  /* STEADY 峰/谷保持预淡出前最小保底时长 (s，加 0.20 s 淡出约 1.50 s) */
#define HIP_GENTLE_MIN_HALF_PERIOD_S       0.25f  /* 状态机步态半周期下限钳位 (s) */
#define HIP_GENTLE_MAX_HALF_PERIOD_S       3.0f   /* 状态机步态半周期上限钳位 (s) */
#define HIP_GENTLE_TIMEOUT_RECOVER_TAU_S   0.10f  /* 未超时阶段 scale_timeout 回升时间常数 (s) */
#define HIP_GENTLE_ENV_FAST_TAU_S          0.074f /* 角差快速退缩时 scale_env 衰减时间常数 (s) */
#define HIP_GENTLE_ENV_RECOVER_TAU_S       0.25f  /* lockin_gain 与 scale_env 恢复时间常数 (s) */

/* ==================== 七、异常步频与停步共振防护 ==================== */
#define HIP_GENTLE_LOW_CADENCE_OMEGA_RAD_S    2.2f  /* 稳态低步频退出判定滤波角频率上限 (rad/s) */
#define HIP_GENTLE_LOW_CADENCE_CONFIRM_S      0.25f /* 稳态低步频连续确认退出时长 (s) */
#define HIP_GENTLE_RESONANCE_OMEGA_RAD_S      11.5f /* 异常高频共振冷却滤波角频率门槛 (rad/s) */
#define HIP_GENTLE_RESONANCE_CONFIRM_S        0.18f /* 异常高频连续确认进入共振冷却时长 (s) */
#define HIP_GENTLE_SHORT_HALF_S               0.30f /* 短半步绝对时长门槛 (s) */
#define HIP_GENTLE_SHRINK_PEAK_MAX_DEG        26.0f /* 两半步缩幅缩时保护触发峰值上限 (deg) */
#define HIP_GENTLE_SHRINK_PEAK_RATIO          0.80f /* 两半步缩幅判定相对同侧参考峰值比例 */
#define HIP_GENTLE_CONTRACT_HALF_RATIO        0.75f /* 两半步缩时判定相对参考半周期比例 */
#define HIP_GENTLE_SHRINK_CONFIRM_HALVES      2U    /* 连续缩幅缩时触发共振冷却的半步数门槛 */
#define HIP_GENTLE_ABRUPT_SHRINK_PEAK_DEG     20.0f /* 单半步突发坍缩直切共振的峰值上限 (deg) */
#define HIP_GENTLE_ABRUPT_SHRINK_PEAK_RATIO   0.62f /* 单半步突发坍缩相对同侧参考峰值比例 */
#define HIP_GENTLE_ABRUPT_SHRINK_HALF_RATIO   0.65f /* 单半步突发坍缩相对参考半周期比例 */
#define HIP_GENTLE_REF_ADAPT_ALPHA            0.15f /* 正常步幅与半周期参考值的指数更新系数 */
#define HIP_GENTLE_EARLY_STOP_ENABLE          true  /* true=开启早峰回落停步快退（识别半步刹车后收腿） */
#define HIP_GENTLE_EARLY_STOP_MIN_REF_DEG     28.0f /* 启用早峰停步判定的同侧参考步幅下限 (deg) */
#define HIP_GENTLE_EARLY_STOP_MAX_PEAK_DEG    33.0f /* 触发早峰挂起的本次峰值角差上限 (deg) */
#define HIP_GENTLE_EARLY_STOP_PEAK_RATIO      0.65f /* 触发早峰挂起相对同侧参考峰值的比例上限 */
#define HIP_GENTLE_EARLY_STOP_PHASE_LEAD_DEG  18.0f /* 触峰相对 90/270 deg 的最小提前相位角 (deg) */
#define HIP_GENTLE_EARLY_STOP_SEVERE_RATIO    0.45f /* 单半步严重缩幅（相对参考峰值）放宽相位提前门槛的比例上限 */
#define HIP_GENTLE_EARLY_STOP_SEVERE_LEAD_DEG 5.0f  /* 单半步严重缩幅时所需的最小提前相位角 (deg) */
#define HIP_GENTLE_COLLAPSE_PEAK_RATIO        0.55f /* 连续两半步大幅缩幅（不要求缩时）的相对参考峰值比例上限 */
#define HIP_GENTLE_COLLAPSE_CONFIRM_HALVES    2U    /* 连续大幅缩幅触发共振冷却的半步数门槛 */
#define HIP_GENTLE_EARLY_STOP_FALL_DEG        6.0f  /* 早峰挂起后确认收腿退出所需的角差回落量 (deg) */
#define HIP_GENTLE_EARLY_STOP_STILL_S         0.10f /* 早峰挂起后双腿真静止确认停步退出时长 (s) */
#define HIP_GENTLE_EARLY_STOP_COOLDOWN_S      0.50f /* 早峰停步退出后的最短冷却防反抽时长 (s) */

typedef enum {
    HIP_GENTLE_IDLE, HIP_GENTLE_FIRST_STEP, HIP_GENTLE_LOCK_IN, HIP_GENTLE_STEADY,
    HIP_GENTLE_COOLDOWN
} hip_gentle_state_id_t;

typedef enum {
    HIP_GENTLE_WAITING_CROSSING = 0, HIP_GENTLE_INVALID_ANGLE = 1, HIP_GENTLE_TIMESTAMP_RESET = 2,
    HIP_GENTLE_QUIET = 4, HIP_GENTLE_NO_PROGRESS = 5, HIP_GENTLE_EVENT_SEED = 9,
    HIP_GENTLE_INVALID_TRACKER = 10, HIP_GENTLE_EVENT_MISMATCH = 13,
    HIP_GENTLE_WAITING_FIRST_PEAK = 19, HIP_GENTLE_HALF_ANCHORED, HIP_GENTLE_LOW_CADENCE,
    HIP_GENTLE_RESONANCE, HIP_GENTLE_INVALID_GEAR_LIMIT
} hip_gentle_reason_t;

typedef enum {
    HIP_GENTLE_RAPID_OFF = 0, HIP_GENTLE_RAPID_POSITIVE = 1, HIP_GENTLE_RAPID_NEGATIVE = 2
} hip_gentle_rapid_channel_t;

typedef struct {
    float response, width, asymmetry;
    float gear_limit_nm; // 当前档位的每侧力矩上限 (N·m)，须为非负有限值
    bool zero_cross_latch, supervision_scale, sync_peak_amplitude, low_cadence_exit;
    float afo_phase_offset_rad; // §4.1 输出相位偏置，默认 0
    bool extended_timing, geometric_phase_pull; // §5.2 扩展时限和可选几何牵引
    bool reversal_envelope; // 可选角差退缩包络；舒缓默认只观测
    float amp_scale_min_rad, amp_scale_full_rad; // 步幅力矩上限 S 曲线起止角差幅值 (rad)
    bool amp_zero_cross_latch; // true=波形过零点锁存步幅上限，保持单半步曲线形状不变
    bool early_stop_enable; // true=开启早峰回落停步快退（半步刹车回落时提前退出并短暂冷却）
} hip_gentle_parameters_t;

typedef struct {
    double left_torque_nm, right_torque_nm;
    hip_gentle_state_id_t state;
    hip_gentle_reason_t reason;
    double phase_rad, frequency_hz, signal_rad, reconstructed_rad, gate; // signal_rad 为原始角差，重构值已去零偏
    double swing_amplitude_rad; // 角度差第一谐波幅值绝对值，不是米制步长
    double target_left_torque_nm;  // 波形经首步守卫和监督后的髋侧目标
    double target_right_torque_nm;
    double adaptive_amplitude_nm; // 当前锁存的波形幅值，N·m
    double carrier_nm; // 髋通道关闭时仍可读取有符号载体
    double filtered_frequency_hz, characteristic_velocity_rad_s; // §4.2 低通后的 ω/(2π)、ωα₁
} hip_gentle_output_t;

typedef struct {
    float eta_alpha, nu_phi, gamma_omega, eta_0;
} hip_gentle_afo_gains_t;

typedef struct {
    bool enabled;
    float omega, omega_filt, velocity_filt_rad_s;
    float d_omega, d_alpha0, d_alpha[5], d_phi[5]; // 上一积分子步的导数
    float alpha0, alpha[5], phi[5], phase_offset_rad, phase_out, y_hat; // 角度、幅值、重构均为 rad
} hip_gentle_afo_t;

typedef struct {
    hip_gentle_state_id_t state;
    hip_gentle_reason_t reason;
    int32_t evt;
    int32_t prev_evt;
    double evt_ts[9];
    bool evt_ts_valid[9];
    double state_enter_ts;
    double static_l_since, static_r_since; // <0 表示当前超出各腿静止带
    double still_l_since, still_r_since; // <0 表示当前超出双侧真静止带
    float prev_y;
    float ma_buf[HIP_GENTLE_DY_WINDOW];
    uint8_t ma_idx;
    float dy_ma;
    float prev_dy_ma;
    float max_dy_pos;
    float min_dy_neg;
    float decel_acc;
    float geo_decel_init;
    float geo_phi_anchor;
    bool return_rebounded; // 回落越阈后重新明显外摆，沿用实测过零观察时限
    bool reacquiring_speed; // 双支撑凹陷解除后允许重新建立本半步速度参考
    hip_gentle_rapid_channel_t rapid_channel; // 跨回阈值或过零后清除
    float phi_ref;
    float seed_peak_pos;
    float seed_peak_neg;
    bool peak_pos_seen;
    bool peak_neg_seen;
    float half_period_s, peak_fade_s;
    float reference_half_s;
    float reference_peak_deg[2];  // 分侧记住稳定半步摆幅，避免把左右不对称当作收步
    unsigned short_halves;
    float high_cadence_elapsed;
    float dir_sign;
    bool phase_override_valid;
    float fsm_phase;
    float scale_timeout;
    float scale_env, lockin_gain;
    float rise_ref_deg_s, fall_ref_deg_s;
    float dyn_to_coeff;
    float amp_confidence;
    float recent_peak_rad; // 最近实测半步摆幅峰值 (rad)，用于步幅上限快降慢升
    float phase_amp_filt_rad, scale_energy; // 相平面瞬时步幅低通估计 (rad) 与步幅能量塌缩平滑倍率
    bool abrupt_shrink_resonance; // 单半步突发大幅缩幅+缩时共振标志
    uint8_t collapse_halves; // 连续大幅缩幅（不要求缩时）的半步计数
    float low_cadence_elapsed;
    uint32_t step_count;
    bool prev_y_valid;
    double observed_zero_ts[2], half_start_ts;
    bool observed_zero_valid[2];
    int half_start_sign;
    bool early_peak_armed; // 检测到明显早峰+缩幅后挂起，等待回落确认停步
    int early_peak_sign; // 挂起半步方向 (+1 正半步 / -1 负半步)
    float early_peak_deg, early_peak_ref_deg; // 挂起时实测峰值与同侧参考步幅 (deg)
    float last_peak_deg[2]; // 分侧记录最近一次正常 STEADY 峰值 (deg)，防止短距步行 reference_peak_deg 滞后
    float cooldown_min_s; // 普通冷却态允许过零回 IDLE 前的最短保持时长 (s)
} hip_gentle_fsm_t;

typedef struct {
    float gear_ramp_nm, prev_zero_sig, resp_ramp;
    float ramp_inner, ramp_outer, amplitude_nm, carrier_nm, target_nm, release_rate_nm_s;
    float amp_scale_latched, prev_wave_shape;
    bool amp_scale_valid;
} hip_gentle_torque_t;

typedef struct {
    hip_gentle_parameters_t parameters;
    bool initialized, enabled, hip_enabled, have_sample;
    uint32_t timestamp_ms;
    double gait_time_s;
    float previous_left, previous_right, left_speed, right_speed;
    float idle_offset_rad; // 仅在 IDLE 双髋静止时更新；AFO 和事件检测共用校正输入
    hip_gentle_afo_t afo;
    hip_gentle_afo_gains_t afo_gains;
    hip_gentle_fsm_t fsm;
    hip_gentle_torque_t torque;
} hip_gentle_assist_t;

void hip_gentle_default_parameters(hip_gentle_parameters_t *parameters);
bool hip_gentle_init(hip_gentle_assist_t *state, const hip_gentle_parameters_t *parameters);
void hip_gentle_reset(hip_gentle_assist_t *state);
/* 控制任务在两帧之间调节。四个百分比依次对应 eta_alpha、nu_phi、gamma_omega、eta_0。
 * trim_percent 为 NULL 时各通道微调为 0。无效输入返回 false，保留上次增益。 */
bool hip_gentle_set_afo_gains(hip_gentle_assist_t *state, float global_scale,
                               const float trim_percent[4]);
/* 为当前外摆半周选择一次急换向通道，传 OFF 撤销。 */
bool hip_gentle_set_rapid_channel(hip_gentle_assist_t *state, hip_gentle_rapid_channel_t channel);
void hip_gentle_set_enable(hip_gentle_assist_t *state, bool enabled, bool hip_enabled);
void hip_gentle_set_response_speed(hip_gentle_assist_t *state, float response);
/* 每帧入口：gear_limit_nm 立即作用于本帧。传 0 时两侧零出力；无效上限使本帧失败并清零。 */
bool hip_gentle_step(hip_gentle_assist_t *state, double left_rad, double right_rad,
                     double left_rad_s, double right_rad_s, uint32_t timestamp_ms,
                     bool valid, float gear_limit_nm, hip_gentle_output_t *output);

/* 分开设档和更新的调用方式也可使用以下三个接口。
 * 初始化后或两帧之间换档；不重置步态。0 表示零力矩。负值或非有限值清零并返回 false。
 * 下一帧起每侧指令绝对值不超过本档上限，正常卸力期间同样生效。reset 保留所选档位。 */
bool hip_gentle_set_gear_limit_nm(hip_gentle_assist_t *state, float gear_limit_nm);
bool hip_gentle_update_with_velocity(hip_gentle_assist_t *state, double left_rad, double right_rad,
                                     double left_rad_s, double right_rad_s, uint32_t timestamp_ms,
                                     bool valid, hip_gentle_output_t *output);
bool hip_gentle_update(hip_gentle_assist_t *state, double left_rad, double right_rad,
                       uint32_t timestamp_ms, bool valid, hip_gentle_output_t *output);
/* 归一化波形，phase 为包含响应偏移的步态相位。 */
double hip_gentle_waveform(double phase, double width, double asymmetry);

#ifdef __cplusplus
}
#endif
#endif
