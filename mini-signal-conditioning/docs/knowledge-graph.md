# Knowledge Graph - Signal Conditioning

## L1: Definitions (Complete)
| # | Definition | C Type | Status |
|---|-----------|--------|--------|
| 1 | Signal conditioning pipeline stages | `sigcond_stage_t` | Implemented |
| 2 | Filter approximation families | `filter_family_t` | Implemented |
| 3 | Active filter topologies | `filter_topology_t` | Implemented |
| 4 | Linearization methods | `linearize_method_t` | Implemented |
| 5 | Sensor excitation types | `excitation_type_t` | Implemented |
| 6 | Isolation barrier types | `isolation_type_t` | Implemented |
| 7 | Bridge configuration types | `bridge_topology_t` | Implemented |
| 8 | Error classification | `sigcond_error_type_t` | Implemented |
| 9 | Sensor input types | `sensor_input_type_t` | Implemented |
| 10 | Filter specifications | `analog_filter_spec_t`, `digital_filter_spec_t` | Implemented |
| 11 | Biquad coefficients | `filter_biquad_t` | Implemented |
| 12 | Polynomial coefficients | `polynomial_coeffs_t` | Implemented |
| 13 | Thermistor model | `steinhart_hart_params_t` | Implemented |
| 14 | RTD model | `rtd_cvd_params_t` | Implemented |
| 15 | Thermocouple NIST model | `thermocouple_nist_model_t` | Implemented |
| 16 | Isolation specification | `isolation_spec_t` | Implemented |
| 17 | Error budget | `error_budget_t` | Implemented |
| 18 | Bridge completion | `bridge_completion_t` | Implemented |
| 19 | CJC configuration | `cjc_config_t` | Implemented |
| 20 | Channel configuration | `sigcond_channel_config_t` | Implemented |

## L2: Core Concepts (Complete)
| # | Concept | Implementation | Status |
|---|---------|---------------|--------|
| 1 | Anti-aliasing filter | `filter_verify_antialias()` | Implemented |
| 2 | Ratiometric measurement | `ratiometric_divider_output()` | Implemented |
| 3 | Galvanic isolation | `isolation_withstand_from_working()` | Implemented |
| 4 | Cold junction compensation | `cjc_compensate()` | Implemented |
| 5 | Wheatstone bridge | `bridge_output_exact()` | Implemented |
| 6 | Lead wire compensation | `bridge_3wire_compensate()` | Implemented |
| 7 | Auto-zero / offset nulling | `selfcal_autozero_correct()` | Implemented |
| 8 | Sensor self-heating | `rtd_self_heating_error()` | Implemented |
| 9 | Excitation stability | `excitation_load_error()` | Implemented |
| 10 | Calibration concepts | `selfcal_two_point()` | Implemented |

## L3: Mathematical Structures (Complete)
| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Laplace domain transfer functions | Butterworth/Chebyshev/Bessel poles |
| 2 | z-domain difference equations | Bilinear transform, FIR/IIR design |
| 3 | Horner polynomial evaluation | `poly_eval_horner()` |
| 4 | Lagrange interpolation | `lut_quadratic_interp()` |
| 5 | Catmull-Rom cubic spline | `lut_cubic_spline()` |
| 6 | Howland current pump model | `howland_current_source()` |
| 7 | Wheatstone bridge transfer function | `bridge_output_exact()` |
| 8 | Barrier capacitance model | `isolation_cap_barrier_cmrr()` |
| 9 | Calibration matrix | `selfcal_multichannel_calibrate()` |

## L4: Fundamental Laws (Complete)
| # | Law | Verification |
|---|-----|--------------|
| 1 | Nyquist-Shannon sampling | `filter_verify_antialias()` |
| 2 | Law of successive temperatures | `cjc_verify_successive_law()` |
| 3 | Kirchhoff voltage law (bridge) | `bridge_is_balanced()` |
| 4 | Ohm's law (excitation) | `excitation_load_error()` |
| 5 | Seebeck effect | NIST ITS-90 implementation |
| 6 | Steinhart-Hart equation | `steinhart_hart_T()` |
| 7 | Callendar-Van Dusen equation | `rtd_cvd_R_from_T()` |
| 8 | IEC creepage/clearance | `isolation_creepage()` |
| 9 | IEC 60751 RTD standard | Pt100/Pt1000 CVD parameters |
| 10 | Bridge balance condition | `bridge_is_balanced()` theorem |

## L5: Algorithms/Methods (Complete)
19 algorithms fully implemented including:
- Butterworth/Chebyshev/Bessel/Elliptic filter design
- Sallen-Key, MFB, Twin-T component design
- Bilinear transform, Kaiser window FIR
- Horner evaluation, Newton-Raphson inversion
- LUT interpolation (linear, quadratic, cubic)
- Howland current pump, bridge optimization
- Two-point/multi-channel calibration
- NIST ITS-90 polynomial evaluation

## L6: Canonical Problems (Complete)
3 end-to-end examples solving real measurement problems.

## L7: Applications (Partial+)
Boeing 787, NASA engine testing, F-35 structural, NHS pharmaceutical.

## L8: Advanced Topics (Partial)
Adaptive filter design, multi-channel calibration matrix, hysteresis-aware linearization.

## L9: Research Frontiers (Partial)
Smart sensor self-calibration documented. AI-based fusion conditioning for future.
