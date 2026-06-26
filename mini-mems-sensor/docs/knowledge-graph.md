# MEMS Sensor Knowledge Graph — L1-L9 Coverage

## L1: Definitions (Complete)

| # | Concept | C Implementation | Lean Definition |
|---|---------|-----------------|-----------------|
| 1 | MEMS Sensor Types (12 types) | `MemsSensorType` enum | `MemsSensorType` inductive |
| 2 | Transduction Mechanisms (8 types) | `MemsTransduceType` enum | `MemsTransduceType` inductive |
| 3 | 3D Vector | `MemsVector3` struct | `MemsVec3` structure |
| 4 | Quaternion | `MemsQuaternion` struct | `MemsQuat` structure |
| 5 | Rotation Matrix 3×3 | `MemsMatrix33` struct | — |
| 6 | MEMS Sensor Specification | `MemsSpec` struct (25 fields) | `MemsSpec` structure |
| 7 | Sensor Timestamped Sample | `MemsSample` struct | — |
| 8 | IMU Configuration | `MemsImuConfig` struct | — |
| 9 | Fabrication Process Types | `MemsFabricationType` enum | — |
| 10 | Seismic Mass Parameters | `MemsSeismicMass` struct | `SpringMassDamper` structure |
| 11 | Spring-Mass-Damper State | `MemsSpringMassDamper` struct | — |
| 12 | Noise Process Types | `MemsNoiseProcess` enum | — |
| 13 | Sensor Axis Enum | `MemsAxis` enum | — |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Capacitive transduction | `mems_capacitive_sensitivity()`, `mems_capacitance_parallel_plate()` |
| 2 | Piezoresistive/piezoelectric transduction | `MemsTransduceType` with full enumeration |
| 3 | Spring-Mass-Damper ODE | `mems_smd_step()`, `MemsSpringMassDamper` |
| 4 | Natural frequency | `mems_resonant_freq()`, `mems_natural_freq_rad()` |
| 5 | Quality factor (Q) | `mems_quality_factor()` |
| 6 | Damping ratio (ζ) | `mems_damping_ratio()` |
| 7 | Squeeze-film damping | `mems_squeeze_film_damping()` |
| 8 | Couette damping | `mems_couette_damping()` |
| 9 | Thermoelastic damping (TED) | `mems_thermoelastic_q_inv()` |
| 10 | Electrostatic pull-in | `mems_pullin_voltage()` |
| 11 | MEMS Fabrication overview | `MemsFabricationType` with 6 processes |
| 12 | Capacitive readout (differential) | `mems_capacitive_displacement()` |

## L3: Mathematical Structures (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Mass-spring-damper ODE solution | `mems_smd_step()` (semi-implicit Euler), harmonic amplitude/phase |
| 2 | Transfer function H(s) | `mems_freq_response_mag()`, `mems_freq_response_phase()` |
| 3 | Laplace domain ↔ time domain | Documented in `mems_dynamics.c` |
| 4 | Power Spectral Density (PSD) | `MemsPsdEstimate`, `mems_welch_psd()` (Welch's method) |
| 5 | Allan Variance (IEEE 1139) | `mems_allan_variance()`, `mems_allan_analysis()` |
| 6 | Allan Deviation noise slopes | `mems_allan_fit_noise()` (all 5 slopes) |
| 7 | Quaternion kinematics | `mems_quat_multiply()`, `mems_quat_rotate()` |
| 8 | Euler angle conversion | `mems_quat_to_euler()`, `mems_euler_to_quat()` |
| 9 | Cross-axis sensitivity matrix | 3×3 scale+misalignment matrices |
| 10 | Coriolis coupling matrix | `mems_coriolis_force()`, gyro scale factor |
| 11 | 3×3 matrix operations | `mat3_det()`, `mat3_inv()`, `mat3_mul()`, `mat3_transpose()` |

## L4: Fundamental Laws (Complete)

| # | Law | C Implementation | Lean Theorem |
|---|-----|-----------------|--------------|
| 1 | Newton's 2nd Law (F=ma) | `mems_smd_step()` core | `newton_second_law_proportional` |
| 2 | Hooke's Law (F=kx) | `mems_static_displacement()` | `static_displacement_k_decreases` |
| 3 | Coriolis Force (F_c=2m·v×Ω) | `mems_coriolis_force()` | `coriolis_force_proportional` |
| 4 | Thermo-mechanical Noise Limit (Gabrielson 1993) | `mems_thermomech_nea()`, `mems_thermomech_ner()` | `nea_temperature_monotone`, `nea_q_decreases` |
| 5 | Johnson-Nyquist Noise | `mems_johnson_noise_density()`, `_rms()` | — |
| 6 | Shot Noise (Schottky 1918) | `mems_shot_noise_density()` | — |
| 7 | Stoney's Formula (1909) | `mems_stoney_stress()`, `mems_stoney_curvature()` | `stoney_stress_curvature_inverse` |
| 8 | Fluctuation-Dissipation Theorem | Documented in `mems_noise.c` | — |
| 9 | Resonance: |H(jω_n)| = 1/(2ζ) | `mems_freq_response_mag(1.0, zeta)` = 1/(2ζ) | `resonance_gain_ge_one` |
| 10 | Barometric Formula (ISA) | `mems_baro_to_altitude()` | `baro_altitude_monotone` |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | 6-Position Accelerometer Calibration | `mems_calib_accel_6pos()` |
| 2 | Rate Table Gyro Calibration (LS) | `mems_calib_gyro_rate()` |
| 3 | Magnetometer Ellipsoid Fitting | `mems_calib_mag_ellipsoid()` |
| 4 | Least-Squares Axis Alignment | `mems_alignment_matrix()` |
| 5 | Polynomial Temperature Compensation | `mems_fit_temp_comp()` (up to 4th order, Gaussian elimination) |
| 6 | Complementary Filter (SO(3)) | `mems_comp_filter_update()` |
| 7 | Mahony Explicit Complementary Filter | `mems_mahony_update()` (proven stability) |
| 8 | Madgwick Gradient Descent AHRS (6+9 axis) | `mems_madgwick_update_imu()`, `_marg()` |
| 9 | 1D Kalman Filter | `mems_kalman1d_predict()`, `_update()` |
| 10 | Allan Variance Full Analysis | `mems_allan_analysis()` + coefficient fitting |
| 11 | Welch's PSD Estimation (FFT-based) | `mems_welch_psd()` |
| 12 | White Noise Generation (Box-Muller) | `mems_generate_white_noise()` |
| 13 | 1/f Flicker Noise Generation | `mems_generate_flicker_noise()` |
| 14 | Xorshift64* PRNG | `mems_rand64()` |

## L6: Canonical Problems (Complete)

| # | Problem | Solution |
|---|---------|----------|
| 1 | AHRS — Attitude & Heading Reference | Mahony, Madgwick, Comp Filter with accel+gyro+mag |
| 2 | Roll/Pitch from Accelerometer | `mems_accel_to_roll_pitch()` |
| 3 | Tilt-Compensated Compass | `mems_mag_to_yaw()` with declination |
| 4 | IMU Dead Reckoning | `mems_dead_reckon()` |
| 5 | IMU Calibration Pipeline (example) | `examples/example_imu_calib.c` |
| 6 | Tilt Sensing (example) | `examples/example_accel_tilt.c` |
| 7 | Sensor Fusion Comparison (example) | `examples/example_sensor_fusion.c` |

## L7: Applications (Complete)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | Smartphone Orientation (iPhone/Android) | `mems_phone_orientation_update()` — Keywords: iPhone, Samsung, Huawei |
| 2 | Drone Stabilization (DJI/ArduPilot) | `mems_drone_attitude_update()` — Keywords: Quadrotor, DJI, ArduPilot |
| 3 | Automotive Airbag (Bosch/Continental) | `mems_airbag_crash_check()` — Keywords: Airbag, Bosch, NCAP |
| 4 | Automotive ESC (Electronic Stability) | `mems_esc_check()` — Keywords: ESC, Bosch ESP, Toyota, Tesla |
| 5 | Wearable Step Counter (Fitbit/Apple Watch) | `mems_pedometer_init()`, `_update()` — Keywords: Fitbit, Apple Watch, Garmin |
| 6 | Industrial Vibration Monitoring (ISO 10816) | `mems_iso_10816_zone()`, `mems_bearing_fault_freqs()` — Keywords: ISO, SKF, Bently Nevada |

## L8: Advanced Topics (Complete)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | GPS-Denied Navigation (ZUPT-aided INS) | `mems_ins_zupt_update()` — Keywords: GPS-denied, indoor, tunnel |
| 2 | Zero-Velocity Detection (ZUPT) | `mems_detect_stationary()` — Foxlin 2005 method |
| 3 | MEMS Micromirror for Structured Light | `mems_mirror_optical_angle()`, `mems_mirror_resonant_freq()` — Keywords: TI DLP, ST LBS |
| 4 | Noise Propagation through 2nd-Order Systems | `mems_noise_propagate_2nd_order()` |

## L9: Research Frontiers (Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | Standard Quantum Limit in MEMS/NEMS | Documented as SQL theorem in Lean (`sql_trace_proportional`) |
| 2 | Chip-Scale Atomic MEMS Sensors | Referenced in course alignment docs |
| 3 | Quantum MEMS (superposition in mechanical resonators) | Referenced in course-tree.md |
