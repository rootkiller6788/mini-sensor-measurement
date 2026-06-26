# Knowledge Graph — mini-sensor-principle

## L1: Definitions — COMPLETE

| Term | Struct/Enum | Location |
|------|------------|----------|
| Sensitivity S = dY/dX | sensor_specs_t.sensitivity | sensor_defs.h |
| Resolution | sensor_specs_t.resolution | sensor_defs.h |
| Accuracy / Precision | sensor_specs_t (all error fields) | sensor_defs.h |
| Hysteresis | sensor_specs_t.hysteresis | sensor_defs.h |
| Linearity (independent/endpoint/terminal) | sensor_linearity_*() | sensor_transfer.h |
| Drift / Offset | sensor_specs_t (offset_drift, gain_drift) | sensor_defs.h |
| Span / Full Scale | sensor_specs_t.full_scale_* | sensor_defs.h |
| Noise Equivalent Power (NEP) | nep_compute() | sensor_noise.h |
| Specific Detectivity D* | detectivity_d_star() | sensor_noise.h |
| Response Time / Bandwidth | sensor_specs_t.bandwidth | sensor_defs.h |
| Gauge Factor GF = (ΔR/R₀)/ε | bridge_strain_*() | sensor_bridge.h |
| Seebeck Coefficient | tc_sensitivity() | sensor_transfer.h |
| Temperature Coefficient α | rtd_cvd_model_t.alpha | sensor_transfer.h |
| Thermal Time Constant | sensor_1st_order_t.tau | sensor_transfer.h |
| Sensor Classification (passive/active) | sensor_energy_type_t | sensor_defs.h |
| Noise Mechanism Types | noise_mechanism_t | sensor_defs.h |

## L2: Core Concepts — COMPLETE

| Concept | Implementation | Location |
|---------|---------------|----------|
| Polynomial Transfer Function | poly_transfer_*() | sensor_transfer.c |
| Thermistor Steinhart-Hart / β-model | thermistor_model_*() | sensor_transfer.c |
| RTD Callendar-Van Dusen (IEC 60751) | rtd_cvd_*() | sensor_transfer.c |
| Thermocouple ITS-90 + CJC | tc_model_*() | sensor_transfer.c |
| 1st-Order Dynamic Response | sensor_1st_order_*() | sensor_transfer.c |
| 2nd-Order Dynamic Response | sensor_2nd_order_*() | sensor_transfer.c |
| Lookup Table Transfer Function | sensor_lut_*() | sensor_transfer.c |
| Wheatstone Bridge (quarter/half/full) | bridge_*() | sensor_bridge.c |
| Lead Wire Error (2/3/4-wire) | bridge_lead_wire_error_2wire() | sensor_bridge.c |
| Bridge Thermal EMF | bridge_thermal_emf_offset() | sensor_bridge.c |
| Self-Heating Effects | bridge_self_heating_*() | sensor_bridge.c |
| Instrumentation Amplifier Gain | amp_instrumentation_gain() | sensor_signal_cond.c |
| Anti-Alias Filter Design | anti_alias_*() | sensor_signal_cond.c |
| 4-20 mA Current Loop | current_loop_*() | sensor_signal_cond.c |
| Sensor Excitation (V/I/AC) | excitation_*() | sensor_signal_cond.c |
| Common-Mode Rejection (CMRR) | bridge_amplifier_cmrr_required() | sensor_bridge.c |

## L3: Mathematical Structures — COMPLETE

| Structure | Implementation | Location |
|-----------|---------------|----------|
| Polynomial (Horner evaluation) | poly_transfer_evaluate() | sensor_transfer.c |
| Newton-Raphson Inverse | poly_transfer_inverse() | sensor_transfer.c |
| Normal Equations (Linear LS) | calibration_linear_ls() | sensor_calibration.c |
| Vandermonde + Gaussian Elim | calibration_polynomial_ls() | sensor_calibration.c |
| Weighted LS Normal Equations | calibration_weighted_linear_ls() | sensor_calibration.c |
| 1/f Noise PSD Integration | flicker_noise_*() | sensor_noise.c |
| Allan Variance (Overlapping) | allan_variance_compute_overlapping() | sensor_noise.c |
| Noise Source RSS Combination | noise_rss_combine() | sensor_noise.c |
| Bridge Thevenin Equivalent | bridge_thevenin_resistance() | sensor_bridge.c |
| Cross-Sensitivity Matrix Inversion | calibration_cross_sensitivity_matrix() | sensor_calibration.c |
| ENOB from Noise | enob_from_noise() | sensor_noise.c |
| Dynamic Error (ramp, sinusoidal) | sensor_dynamic_*_error_*() | sensor_transfer.c |
| Uncertainty Budget (GUM) | sensor_uncertainty_budget_compute() | sensor_defs.c |
| Piecewise Linear Interpolation | sensor_lut_evaluate() | sensor_transfer.c |

## L4: Fundamental Laws — COMPLETE

| Law / Theorem | Implementation | Formalized (Lean) |
|---------------|---------------|-------------------|
| Johnson-Nyquist Thermal Noise (1928) | johnson_noise_voltage() | YES — `johnson_noise_scaling` |
| Shot Noise (Schottky, 1918) | shot_noise_current() | — |
| Seebeck Effect (1821) | tc_emf_from_T() | — |
| Steinhart-Hart Equation (1968) | thermistor_temp_from_R() | YES — `steinhart_hart_monotonic` |
| Callendar-Van Dusen (1887) | rtd_R_from_T() | — |
| Friis Noise Cascade Formula (1944) | noise_factor_cascade() | — |
| Wheatstone Bridge Balance (1843) | bridge_output_voltage_exact() | YES — `BridgeBalanced` |
| First-Order ODE Step Response | sensor_1st_order_step_response() | YES — `first_order_monotonic` |
| Second-Order Damped Response | sensor_2nd_order_step_response() | — |
| Fluctuation-Dissipation Theorem | johnson_noise_power_available() | — |
| Allan Variance Power-Law Dependence | allan_noise_coefficients_identify() | — |
| GUM Uncertainty Combination | sensor_uncertainty_budget_compute() | — |
| Rose Detection Criterion (1948) | rose_criterion_check() | — |
| Noise Additivity (Uncorrelated) | noise_rss_combine() | YES — `noise_additive_monotonic` |
| Temperature Coefficient of Resistance | rtd_cvd_init_iec751() | YES — `rtd_monotonic` |

## L5: Algorithms/Methods — COMPLETE

| Algorithm | Function | Location |
|-----------|----------|----------|
| Linear Least-Squares | calibration_linear_ls() | sensor_calibration.c |
| Polynomial LS (Vandermonde) | calibration_polynomial_ls() | sensor_calibration.c |
| Weighted Least-Squares | calibration_weighted_linear_ls() | sensor_calibration.c |
| Newton-Raphson (Polynomial Inverse) | poly_transfer_inverse() | sensor_transfer.c |
| Newton-Raphson (TC Temperature) | tc_T_from_emf() | sensor_transfer.c |
| Newton-Raphson (Thermistor R→T) | thermistor_R_from_temp() | sensor_transfer.c |
| Gaussian Elimination | solve_linear_system() [static] | sensor_calibration.c |
| Gauss-Jordan (Matrix Inversion) | calibration_cross_sensitivity_matrix() | sensor_calibration.c |
| Temperature Compensation | calibration_temp_compensation() | sensor_calibration.c |
| Cross-Sensitivity Compensation | calibration_cross_sensitivity_comp() | sensor_calibration.c |
| Auto-Zero / CDS | calibration_auto_zero() | sensor_calibration.c |
| Moving Average Filter | signal_moving_average() | sensor_signal_cond.c |
| Exponential (IIR) Filter | signal_exponential_filter() | sensor_signal_cond.c |
| Median Filter | signal_median_filter() | sensor_signal_cond.c |
| Allan Variance (Overlapping) | allan_variance_compute_overlapping() | sensor_noise.c |
| Noise Coefficient Identification | allan_noise_coefficients_identify() | sensor_noise.c |
| Piecewise Linearization (3-point) | calibration_piecewise_3point() | sensor_calibration.c |
| Anti-Alias Filter Order Selection | anti_alias_filter_order_required() | sensor_signal_cond.c |
| Thermistor Linearization (Rp) | thermistor_linearization_Rp() | sensor_transfer.c |

## L6: Canonical Problems — COMPLETE

| Problem | Example/Implementation |
|---------|----------------------|
| Thermocouple Measurement with CJC | ex_thermocouple_cjc.c |
| Strain Gauge Bridge Measurement | ex_strain_gauge_bridge.c |
| MEMS Accelerometer Dynamic Response | ex_accelerometer_sim.c |
| RTD Temperature Measurement | test_rtd_cvd in test_sensor_core.c |
| Thermistor Temperature Sensing | test_thermistor_beta in test_sensor_core.c |
| Wheatstone Bridge Nonlinearity Correction | bridge_strain_nonlinearity_corrected() |
| Sensor Calibration (multi-point polynomial) | calibration_polynomial_ls() |
| Dynamic Error in Ramp Tracking | sensor_dynamic_ramp_error_1st/2nd_order() |
| Noise Floor Characterization (Allan Var) | allan_variance_compute_overlapping() |
| 4-20 mA Industrial Loop | current_loop_*() in sensor_signal_cond.c |

## L7: Applications — COMPLETE (3+ applications)

| Application | Keywords | Evidence |
|-------------|----------|----------|
| Industrial Temperature Monitoring | thermocouple, RTD, 4-20mA | ex_thermocouple_cjc.c, current_loop_*() |
| Automotive / Aerospace Strain Gauging | Boeing, Detroit, F-35 | ex_strain_gauge_bridge.c |
| MEMS Inertial (GPS/IMU) | Accelerometer, gyro, quadrotor | ex_accelerometer_sim.c |
| Structural Health Monitoring | strain-to-stress, Hooke's Law | bridge_strain_from_bridge_output() |
| Industrial Process Control (ISO 5167) | √ΔP flow measurement, 4-20mA | linearize_sqrt_law() |

## L8: Advanced Topics — PARTIAL (2/5 topics)

| Topic | Evidence | Status |
|-------|----------|--------|
| Sensor Noise Figure Cascade (Friis) | noise_factor_cascade() | Complete |
| Allan Variance Multi-Source Decomposition | allan_noise_coefficients_identify() | Complete |
| MEMS Sensor Fusion | docs only | Partial |
| Self-Calibrating Sensors | docs only | Partial |
| Time-Varying Sensor Models | docs only | Partial |

## L9: Research Frontiers — PARTIAL

| Topic | Evidence | Status |
|-------|----------|--------|
| Smart Sensors / IoT Sensors | docs only | Partial |
| Energy-Harvesting Sensors | docs only | Partial |
| Quantum Sensor Noise Limits | docs (BLIP D* limit) | Partial |
