# Knowledge Graph — mini-smart-sensor

Nine-level knowledge coverage map for Smart Sensor Measurement, Calibration, and Fusion.

## L1: Definitions (Complete)

| Item | C Implementation | Lean Formalization |
|------|-----------------|-------------------|
| Transducer types (14 types) | `ss_transducer_t` enum | `TransducerType` inductive |
| Measurand categories (17 types) | `ss_measurand_t` enum | `Measurand` inductive |
| Communication interfaces (15 types) | `ss_interface_t` enum | — |
| Operating modes (8 states) | `ss_operating_mode_t` enum | `OperatingMode` inductive |
| TEDS template types | `ss_teds_template_t` enum | — |
| Sensor specification | `ss_sensor_spec_t` struct | `SensorSpec` structure |
| Accuracy specification | `ss_accuracy_spec_t` struct | — |
| Sensor runtime state | `ss_sensor_state_t` struct | `SensorState` structure |
| TEDS record (IEEE 1451.4) | `ss_teds_t` struct | — |
| Measurement record | `ss_measurement_t` struct | `Measurement` structure |
| Measurement quality flags | `SS_FLAG_*` macros | — |
| Analog front-end config | `ss_analog_frontend_t` struct | — |
| Digital processing config | `ss_digital_processing_t` struct | — |
| Complete sensor config | `ss_sensor_config_t` struct | — |
| Multi-channel sensor node | `ss_node_t` struct | `SensorNode` structure |
| Linear calibration model | `ss_linear_model_t` struct | `LinearCalibration` structure |
| Polynomial calibration model | `ss_polynomial_model_t` struct | — |
| Steinhart-Hart thermistor model | `ss_steinhart_hart_t` struct | — |
| Thermocouple model (ITS-90) | `ss_thermocouple_model_t` struct | — |
| Allan variance structure | `ss_allan_variance_t` struct | — |
| 1D Kalman filter state | `ss_kalman_1d_t` struct | `KalmanState1D` structure |
| 3D IMU Kalman filter state | `ss_kalman_3d_imu_t` struct | — |
| Complementary filter state | `ss_complementary_filter_t` struct | — |

## L2: Core Concepts (Complete)

| Concept | Implementation |
|---------|---------------|
| Sensor transduction principles | `ss_sensor_spec_init()` with type-specific defaults |
| Signal conditioning chain | `ss_analog_frontend_init()`, ADC/amplifier/filter config |
| Anti-aliasing filter design | `ss_analog_frontend_t`, Butterworth/Bessel/Chebyshev |
| ADC oversampling for ENOB | `adc_oversampling_ratio` field |
| Digital post-processing chain | `ss_digital_processing_init()`, moving avg/median/IIR |
| Sensor state machine | `ss_sensor_state_set_mode()` with transition validation |
| Measurement ring buffer | `ss_measurement_buffer_*()` operations |
| Multi-channel node management | `ss_node_init()`, `ss_node_add_channel()` |
| IEEE 1451.4 TEDS | `ss_teds_t` with calibration data |
| Running statistics (Welford) | `ss_sensor_state_update_stats()` |
| Measurement quality assessment | `ss_measurement_t.flags` with 11 quality flags |
| Calibration validity tracking | `ss_sensor_state_t.calibration_valid` |

## L3: Mathematical Structures (Complete)

| Structure | Implementation |
|-----------|---------------|
| Wheatstone bridge circuit model | `ss_cond_wheatstone_voltage()` |
| Linear sensor model (y = ax + b) | `ss_linear_model_t` + `ss_cal_linear_apply()` |
| Polynomial sensor model (Horner) | `ss_polynomial_model_t` + `ss_cal_polynomial_apply()` |
| Steinhart-Hart equation | `ss_cal_steinhart_hart_temp()` |
| Beta parameter thermistor model | `ss_cal_thermistor_beta()` |
| ITS-90 thermocouple polynomials | `ss_cal_thermocouple_to_temp()` |
| Generic transfer function models | `ss_cond_transfer_function()` (5 model types) |
| Sensitivity (derivative) analysis | `ss_cond_sensitivity_analysis()` |
| ADC quantization noise model | `ss_cond_adc_quantization_noise()` |
| ADC ENOB model | `ss_cond_adc_enob()` |
| Dynamic range formula | `ss_cond_dynamic_range_db()` |
| SNR fundamentals | `ss_cond_snr_db()`, `ss_cond_min_detectable_signal()` |
| Kalman state-space model | `ss_kalman_1d_t`, `ss_kalman_3d_imu_t` |
| Complementary filter transfer function | `ss_complementary_filter_t` |
| Inverse-variance weighting | `ss_fusion_weighted_average()` |
| Mahalanobis distance (1D) | `ss_fusion_mahalanobis_1d()` |
| Allan variance statistical model | `ss_allan_variance_t` + computation |

## L4: Fundamental Laws (Complete)

| Law/Theorem | C Verification | Lean Statement |
|-------------|---------------|----------------|
| Wheatstone bridge balance condition | `test_wheatstone_voltage()` | `wheatstoneBalanced`, `wheatstone_balance_symmetric` |
| Gauss-Markov theorem (OLS is BLUE) | `test_linear_regression()` | `LinearCalibration` structure |
| Steinhart-Hart thermistor equation | `test_steinhart_hart_3pt()` | — |
| Inverse-variance weighting optimality | `test_weighted_average()` | `fuseTwoSensors`, `fusedVariance` |
| Kalman filter optimality (MVUE) | `test_kalman_1d_step()` | `KalmanState1D`, `kalman1DUpdate` |
| Shannon-Nyquist sampling theorem | `ss_analog_frontend_t` anti-aliasing filter | — |
| Propagation of uncertainty (GUM) | `ss_cal_combined_uncertainty()` | — |
| ADC quantization noise (Bennett) | `test_adc_quantization_noise()` | `adcQuantizationNoise` |
| Ideal gas law (TPMS pressure comp.) | `test_tpms_compensate()` | `compensatePressure`, `no_compensation_when_equal_temp` |
| Ohm's law for resistive sensors | `test_strain_from_bridge()` | — |
| SNR definition and MDS | `test_snr_db()`, `test_min_detectable_signal()` | `signalToNoiseRatioDB` |
| Conservation of energy (harvesting) | `ss_advanced_max_duty_cycle()` | — |

## L5: Algorithms/Methods (Complete)

| Algorithm | Implementation | Complexity |
|-----------|---------------|------------|
| Linear regression (OLS) | `ss_cal_linear_regression()` | O(n) |
| Polynomial regression (normal eq.) | `ss_cal_polynomial_regression()` | O(n·d² + d³) |
| Gaussian elimination w/ pivoting | `gauss_solve()` | O(d³) |
| Horner's method (polynomial eval) | `ss_cal_polynomial_apply()` | O(d) |
| Newton-Raphson root finding | `ss_cal_polynomial_inverse()` | O(d·iter) |
| Steinhart-Hart 3-point (Cramer) | `ss_cal_steinhart_hart_3pt()` | O(1) |
| Beta parameter 2-point | `ss_cal_thermistor_beta()` | O(1) |
| ITS-90 thermocouple inversion | `ss_cal_thermocouple_to_temp()` | O(d) |
| Two-point calibration | `ss_cal_two_point()` | O(1) |
| Binary search + linear interpolation | `ss_cal_table_lookup()`, `ss_cond_piecewise_linearize()` | O(log n) |
| Simple moving average | `ss_cond_moving_average()` | O(N) |
| Running moving average (O(1)) | `ss_cond_running_moving_avg()` | O(1) |
| Median filter (insertion sort) | `ss_cond_median_filter()` | O(N²) small N |
| IIR lowpass filter | `ss_cond_iir_lowpass()` | O(1) |
| IIR highpass filter | `ss_cond_iir_highpass()` | O(1) |
| Peak detector (exp. decay) | `ss_cond_peak_detector()` | O(1) |
| Rate-of-change estimator | `ss_cond_rate_of_change()` | O(1) |
| Schmitt trigger hysteresis | `ss_cond_threshold_hysteresis()` | O(1) |
| Temperature compensation | `ss_cond_temp_compensate()` | O(1) |
| Poly nonlinearity correction | `ss_cond_poly_nonlinearity_correct()` | O(d) |
| Kalman filter 1D | `ss_fusion_kalman_1d_step()` | O(1) |
| Kalman filter IMU (2-state) | `ss_fusion_kalman_imu_predict/update()` | O(1) |
| Complementary filter | `ss_fusion_complementary_update()` | O(1) |
| Weighted average (inverse-variance) | `ss_fusion_weighted_average()` | O(n) |
| Median fusion (qsort) | `ss_fusion_median()` | O(n·log n) |
| Trimmed mean | `ss_fusion_trimmed_mean()` | O(n·log n) |
| Exponential moving average | `ss_fusion_ema()` | O(1) |
| Cumulative moving average | `ss_fusion_cumulative_avg()` | O(1) |
| Mahalanobis outlier rejection | `ss_fusion_outlier_rejection()` | O(n) |
| Allan variance computation | `ss_advanced_allan_variance()` | O(N²/k) |
| RMS computation | `ss_advanced_signal_rms()` | O(n) |
| Crest factor | `ss_advanced_crest_factor()` | O(n) |
| Skewness computation | `ss_advanced_skewness()` | O(n) |
| Spectral centroid | `ss_advanced_spectral_centroid()` | O(n) |
| Zero-crossing rate | `ss_advanced_zero_crossing_rate()` | O(n) |

## L6: Canonical Problems (Complete)

| Problem | Implementation |
|---------|---------------|
| Thermistor calibration (Steinhart-Hart + Beta) | `example_thermistor_cal.c` |
| IMU orientation from gyro + accelerometer | `example_imu_fusion.c` |
| Multi-sensor statistical fusion with outlier rejection | `example_imu_fusion.c` Part 3 |
| Vibration severity classification (ISO 10816) | `test_vibration_severity()`, `test_vibration_zone()` |
| Tire pressure temperature compensation | `test_tpms_compensate()`, `test_tpms_deviation()` |
| Slow leak detection from pressure trend | `test_tpms_slow_leak()` |
| SpO2 estimation from PPG ratio-of-ratios | `test_spo2_estimate()` |
| PM2.5 AQI computation | `test_aqi_pm25()` |
| Multi-channel IoT sensor node | `example_smart_node.c` |

## L7: Applications (Complete)

| Application | Domain | Key Functions |
|-------------|--------|---------------|
| Industrial vibration monitoring (ISO 13374/10816) | Industry 4.0 | `ss_app_vibration_severity`, `ss_app_bpfo_frequency` |
| Automotive TPMS (ISO 21750, FMVSS 138) | Automotive | `ss_app_tpms_compensate_pressure`, `ss_app_tpms_detect_slow_leak` |
| Medical pulse oximeter signal processing | Healthcare | `ss_app_spo2_estimate`, `ss_app_heart_rate_from_ppg` |
| Environmental air quality monitoring (EPA AQI) | Smart City | `ss_app_aqi_pm25`, `ss_app_dew_point`, `ss_app_heat_index` |

## L8: Advanced Topics (Complete)

| Topic | Implementation |
|-------|---------------|
| Energy-harvesting power budget (PV, TEG) | `ss_advanced_pv_power_uw`, `ss_advanced_teg_power` |
| Duty-cycling for energy-neutral operation | `ss_advanced_avg_power_uw`, `ss_advanced_max_duty_cycle` |
| Allan variance for MEMS noise analysis | `ss_advanced_allan_variance`, `ss_advanced_allan_fit` |
| Edge-ML feature extraction (RMS, crest, skewness, ZCR) | `ss_advanced_signal_rms`, `ss_advanced_crest_factor`, etc. |
| Mahalanobis distance outlier detection | `ss_fusion_mahalanobis_1d`, `ss_fusion_outlier_rejection` |

## L9: Research Frontiers (Partial)

| Topic | Status |
|-------|--------|
| 6G RIS (Reconfigurable Intelligent Surface) smart sensors | Documented |
| Quantum sensor interfaces | Documented |
| Bio-inspired neuromorphic smart sensors | Documented |
| AI-native sensor edge computing (TinyML) | Referenced in L8 edge-ML |
| Terahertz CMOS sensor interfaces | Documented |
