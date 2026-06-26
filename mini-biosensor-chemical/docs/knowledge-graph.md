# Knowledge Graph — mini-biosensor-chemical

## L1: Definitions (COMPLETE)

| # | Knowledge Item | C Type | Lean Type | Status |
|---|---------------|--------|-----------|--------|
| 1.1 | BiosensorDescriptor (bioreceptor + transducer) | `BiosensorDescriptor` struct | `Biosensor` structure | ✅ |
| 1.2 | Bioreceptor types (enzyme, Ab, DNA, cell, MIP...) | `BioreceptorType` enum | `BioreceptorType` inductive | ✅ |
| 1.3 | Transducer types (ampero, potentio, impedi, optical...) | `TransducerType` enum | `TransducerType` inductive | ✅ |
| 1.4 | Measurement modes (end-point, kinetic, differential...) | `MeasurementMode` enum | `MeasurementMode` inductive | ✅ |
| 1.5 | Michaelis-Menten parameters (Km, Vmax, kcat) | `MichaelisMentenParams` struct | `MichaelisMentenNat` struct | ✅ |
| 1.6 | Hill parameters (Kd, n_Hill, Rmax) | `HillParams` struct | `HillModel` struct | ✅ |
| 1.7 | Langmuir parameters (k_ads, k_des, q_max, K_eq) | `LangmuirParams` struct | `LangmuirModel` struct | ✅ |
| 1.8 | Nernst parameters (E°, n, T) | `NernstParams` struct | `NernstParams` struct | ✅ |
| 1.9 | Butler-Volmer parameters (i₀, α, A, η) | `ButlerVolmerParams` struct | — | ✅ |
| 1.10 | EIS parameters (R_s, R_ct, C_dl, σ) | `EISParams` struct | — | ✅ |
| 1.11 | ISFET parameters (V_th, sensitivity, g_m) | `ISFETParams` struct | — | ✅ |
| 1.12 | SPR parameters (λ, RI sensitivity, d_Au) | `SPRParams` struct | — | ✅ |
| 1.13 | Fluorescence parameters (λ_ex, λ_em, Φ, ε, R₀) | `FluorescenceParams` struct | — | ✅ |
| 1.14 | QCM parameters (f₀, A, Sauerbrey sensitivity) | `QCMParameters` struct | — | ✅ |
| 1.15 | Calibration standards + results | `CalibrationStandards`, `LinearCalibration` | `Calibration` struct | ✅ |
| 1.16 | Sensor array (chemiresistor, e-nose) | `SensorArray`, `ChemiresistorElement` | — | ✅ |
| 1.17 | Measurement data containers | `MeasurementPoint`, `MeasurementSeries` | — | ✅ |
| 1.18 | 4PL/5PL logistic models | `FourPLLogisticFit`, `FivePLLogisticFit` | — | ✅ |

## L2: Core Concepts (COMPLETE)

| # | Knowledge Item | Implementation |
|---|---------------|----------------|
| 2.1 | Amperometric biosensing principle | `glucose_biosensor_current()`, `lactate_biosensor_current()` |
| 2.2 | Potentiometric biosensing (ISE/ISFET) | `nernst_potential()`, `isfet_compute_ph()`, `nikolskii_eisenman_potential()` |
| 2.3 | Impedimetric biosensing (EIS) | `eis_randles_impedance()`, `eis_fit_charge_transfer_resistance()` |
| 2.4 | Optical biosensing (absorbance, fluorescence, SPR) | `spr_resonance_angle()`, `fluorescence_intensity()`, `fiber_evanescent_absorbance()` |
| 2.5 | Piezoelectric mass sensing (QCM) | QCM parameters in types, Sauerbrey sensitivity calc |
| 2.6 | Enzyme-based biorecognition | Full Michaelis-Menten + inhibition models |
| 2.7 | Immunoassay-based biorecognition | `antibody_antigen_complex()`, `binding_association_curve()` |
| 2.8 | DNA-based biorecognition | `dna_hybridization_efficiency()`, `dna_probe_melting_temperature()` |
| 2.9 | Signal generation in biosensors | `glucose_biosensor_current()`, `fluorescence_intensity()` |
| 2.10 | Calibration principles | All calibration functions in `biosensor_calibration.h` |

## L3: Mathematical Structures (COMPLETE)

| # | Knowledge Item | Implementation |
|---|---------------|----------------|
| 3.1 | Michaelis-Menten equation | `michaelis_menten_rate()` - rational function form |
| 3.2 | Hill equation | `hill_fractional_saturation()` - sigmoid with cooperativity |
| 3.3 | Langmuir isotherm | `langmuir_coverage()` - hyperbolic binding |
| 3.4 | Nernst equation | `nernst_potential()` - logarithmic concentration dependence |
| 3.5 | Butler-Volmer equation | `butler_volmer_current()` - exponential overpotential |
| 3.6 | Cottrell equation | `cottrell_current()` - t^(-1/2) diffusion decay |
| 3.7 | Randles-Ševčík equation | `randles_sevcik_peak_current()` - √v dependence |
| 3.8 | Beer-Lambert law (vector form) | Two-component matrix solution in `two_component_absorbance()` |
| 3.9 | FRET distance relation | `fret_efficiency()`, `fret_distance()` - 1/r⁶ dependence |
| 3.10 | Fickian diffusion | `fickian_concentration()` - erfc profile |
| 3.11 | Damköhler number | `damkohler_number()` - dimensionless analysis |
| 3.12 | Debye-Hückel theory | `debye_huckel_activity_coeff()` |
| 3.13 | Fresnel multi-layer optics | `spr_reflectivity()` - transfer matrix method |
| 3.14 | 4PL/5PL logistic models | `fourpl_fit()`, `fivepl_fit()` |

## L4: Fundamental Laws (COMPLETE)

| # | Law | C Implementation | Lean Theorem |
|---|-----|-----------------|--------------|
| 4.1 | Nernst equation | `nernst_potential()` | `nernst_potential` definition |
| 4.2 | Butler-Volmer equation | `butler_volmer_current()` | — |
| 4.3 | Cottrell equation | `cottrell_current()` | — |
| 4.4 | Randles-Ševčík equation | `randles_sevcik_peak_current()` | — |
| 4.5 | Beer-Lambert law | `beer_lambert_absorbance()` | `beer_lambert_additivity`, `beer_lambert_zero_conc` |
| 4.6 | Michaelis-Menten equation | `michaelis_menten_rate()` | `mm_rate_zero_at_zero_substrate`, `mm_rate_nonneg`, `mm_rate_bounded_by_vmax` |
| 4.7 | Hill equation | `hill_fractional_saturation()` | `classifyHill`, `classify_hill_exhaustive` |
| 4.8 | Langmuir isotherm | `langmuir_coverage()` | `langmuir_zero_conc`, `langmuir_bounded_mille` |
| 4.9 | Fick's laws of diffusion | `fickian_concentration()` | — |
| 4.10 | FRET efficiency | `fret_efficiency()` | `fret_efficiency_at_r0`, `fret_efficiency_le_one`, `fret_efficiency_positive` |
| 4.11 | Sauerbrey equation (QCM) | Mass sensitivity in `qcm_params_init()` | — |
| 4.12 | Nikolskii-Eisenman equation | `nikolskii_eisenman_potential()` | — |
| 4.13 | Wallace rule (DNA Tm) | `dna_probe_melting_temperature()` | `wallaceTm`, `gc_rich_higher_tm`, `wallace_tm_positive` |

## L5: Algorithms/Methods (COMPLETE)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 5.1 | LOD estimation (3σ method) | `lod_three_sigma()` |
| 5.2 | LOQ estimation (10σ method) | `loq_ten_sigma()` |
| 5.3 | LOD estimation (ISO 11843) | `lod_iso_11843()` |
| 5.4 | ALS baseline correction | `als_baseline_correction()` |
| 5.5 | Moving-minimum baseline | `moving_minimum_baseline()` |
| 5.6 | Exponential drift correction | `exponential_drift_correction()` |
| 5.7 | Savitzky-Golay smoothing | `savitzky_golay_smooth()` |
| 5.8 | Median filter (spike removal) | `median_filter_spike_removal()` |
| 5.9 | SNR computation | `compute_snr()` |
| 5.10 | Linear regression (OLS) | `linear_calibration_fit()` |
| 5.11 | Weighted least squares (WLS) | `weighted_linear_calibration_fit()` |
| 5.12 | 4PL Levenberg-Marquardt fit | `fourpl_fit()` |
| 5.13 | 5PL fit (grid search + 4PL) | `fivepl_fit()` |
| 5.14 | Log-logistic linearization (4PL init) | `fourpl_linearized_initial_guess()` |
| 5.15 | Standard addition method | `standard_addition_method()` |
| 5.16 | Grubbs outlier test | `grubbs_outlier_test()` |
| 5.17 | Dixon Q-test | `dixon_q_test()` |
| 5.18 | Coefficient of variation (CV%) | `coefficient_of_variation()` |
| 5.19 | Lineweaver-Burk linearization | `lineweaver_burk_estimate()` |
| 5.20 | Eadie-Hofstee linearization | `eadie_hofstee_estimate()` |
| 5.21 | Hanes-Woolf linearization | `hanes_woolf_estimate()` |
| 5.22 | Michaelis-Menten NLS (Gauss-Newton) | `michaelis_menten_nls_fit()` |
| 5.23 | Hill plot analysis | `hill_plot_estimate()` |
| 5.24 | Mandel lack-of-fit test | `mandel_fitting_test()` |
| 5.25 | Durbin-Watson test | `durbin_watson_test()` |
| 5.26 | CV peak detection | `detect_voltammetry_peaks()` |
| 5.27 | Coulometric integration | `coulometric_charge()` |
| 5.28 | Lagrange interpolation for peak position | (embedded in `detect_voltammetry_peaks()`) |
| 5.29 | Abramowitz-Stegun erf approximation | (embedded in `fickian_concentration()`) |
| 5.30 | Thomas algorithm (tridiagonal solver) | (embedded in `als_baseline_correction()`) |

## L6: Canonical Problems (COMPLETE)

| # | Problem | Implementation |
|---|---------|---------------|
| 6.1 | Glucose biosensor amperometric calibration | `glucose_biosensor_current()`, example_glucose_biosensor.c |
| 6.2 | pH-ISFET calibration and measurement | `isfet_compute_ph()`, `nernst_potential()` |
| 6.3 | ELISA immunoassay calibration (4PL) | `fourpl_fit()`, `elisa_fourpl_predict()`, example_elisa_calibration.c |
| 6.4 | Cyclic voltammetry peak detection | `detect_voltammetry_peaks()`, `randles_sevcik_peak_current()` |
| 6.5 | Electrochemical impedance spectroscopy | `eis_randles_impedance()`, `eis_fit_charge_transfer_resistance()` |
| 6.6 | Ion-selective electrode (ISE) interference | `nikolskii_eisenman_potential()`, `debye_huckel_activity_coeff()` |
| 6.7 | SPR immunoassay binding kinetics | `spr_resonance_angle()`, `binding_association_curve()` |
| 6.8 | Enzyme inhibition analysis | `competitive_inhibition_rate()`, `uncompetitive_inhibition_rate()`, `noncompetitive_inhibition_rate()` |
| 6.9 | Two-point sensor recalibration | `two_point_calibration_adjust()` |
| 6.10 | Standard addition for complex matrices | `standard_addition_method()` |
| 6.11 | Sensor aging and lifetime prediction | `estimate_sensor_remaining_life()` |

## L7: Applications (COMPLETE — 5 applications)

| # | Application | Implementation |
|---|-------------|---------------|
| 7.1 | Point-of-care glucose monitoring | `glucose_biosensor_current()`, glucose clinical classification |
| 7.2 | Clinical ELISA / immunoassay data reduction | `elisa_fourpl_predict()`, `elisa_weighted_mean()`, example_elisa_calibration.c |
| 7.3 | Point-of-care lactate for sepsis/sports | `lactate_biosensor_current()` |
| 7.4 | DNA microarray gene expression | `dna_hybridization_efficiency()`, `dna_probe_melting_temperature()` |
| 7.5 | Electronic nose for gas/odor classification | `SensorArray`, `sensor_array_normalize()`, example_sensor_array.c |
| 7.6 | Heavy metal detection in water | `standard_addition_method()`, `nikolskii_eisenman_potential()` |
| 7.7 | Fiber-optic biosensor for remote sensing | `fiber_evanescent_absorbance()`, `evanescent_power_fraction()` |
| 7.8 | Drug discovery SPR screening | `spr_adlayer_thickness()`, `binding_association_curve()` |
| 7.9 | Blood gas analyzer QC | `two_point_calibration_adjust()`, LOD/LOQ calculations |
| 7.10 | Environmental odor monitoring | `sensor_array_euclidean_distance()`, e-nose pattern matching |

## L8: Advanced Topics (PARTIAL+ — 3 topics)

| # | Topic | Implementation |
|---|-------|---------------|
| 8.1 | Nanomaterial-based chemiresistor arrays | `SensorArray`, `ChemiresistorElement`, MOS sensor modeling |
| 8.2 | Calibration transfer (PDS) | `pds_calibration_transfer()` |
| 8.3 | Microfluidic mass transport (Damköhler) | `damkohler_number()` |
| 8.4 | Aptamer-based biosensors | (modeled via Langmuir + DNA hybridization) |
| 8.5 | Multi-analyte sensor fusion | `two_component_absorbance()`, sensor array distance |

## L9: Research Frontiers (PARTIAL — documented)

| # | Topic | Status |
|---|-------|--------|
| 9.1 | CRISPR-based detection (SHERLOCK/DETECTR) | Documented in gap-report.md; not yet implemented |
| 9.2 | Single-molecule digital ELISA (Simoa) | Documented; Poisson statistics needed |
| 9.3 | Wearable sweat sensors (glucose, lactate, Na+) | Referenced in L7 applications |
| 9.4 | Nanopore sequencing biosensors | `BioreceptorType.nanopore` enumerated |
| 9.5 | Organ-on-a-chip integrated biosensing | Documented; microfluidic models partial |
