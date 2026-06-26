# Knowledge Graph - mini-instrumentation-amp

## L1: Definitions (Complete)
- `ia_topology_t` - 7 IA topologies
- `ia_coupling_t` - 5 coupling modes
- `ia_gain_mode_t` - 5 gain programming methods
- `ia_error_type_t` - 10 error types
- `ia_spec_t` - Complete IA specification struct
- `bridge_sensor_t` - Wheatstone bridge model
- `bridge_config_t` - 6 bridge configurations
- `sensor_type_t` - 14 sensor types
- `ia_state_t` - Operating state
- `ia_calibration_t` - Calibration coefficients
- `ia_error_code_t` - 10 error codes
- `thermocouple_type_t` - 8 thermocouple types
- `ia_three_opamp_config_t` - 3-op-amp config
- `ia_two_opamp_config_t` - 2-op-amp config
- `ia_current_mode_config_t` - Current-mode config
- `ia_flying_cap_config_t` - Flying-capacitor config
- `ia_indirect_current_config_t` - Indirect CFB config

## L2: Core Concepts (Complete)
- 3-op-amp IA topology and gain equation
- 2-op-amp IA topology and CMRR zero
- Current-mode IA (CCII-based)
- Flying-capacitor auto-zero topology
- Indirect current feedback topology
- DDA and CCII topologies (documented)
- CMRR resistor mismatch limitation
- R_G selection for target gain
- Output swing and headroom validation
- Topology selection based on sensor/supply
- E96/E24 standard resistor value selection
- Auto-range gain control with hysteresis

## L3: Mathematical Structures (Complete)
- Johnson-Nyquist thermal noise model
- Shot noise (Schottky) model
- 1/f (flicker) noise model
- Total RTI noise with all sources
- Noise figure (NF) computation
- Cascaded noise figure (Friis formula)
- RSS error budget synthesis
- Gain error budget (R_G tolerance + A_OL)
- SNR, MDS, dynamic range
- CMRR vs frequency roll-off model
- CMRR from impedance imbalance
- Combined CMRR from multiple sources
- Sensitivity analysis (gain to R_G, tempco)
- Offset drift, bias current drift
- Optimal source resistance for NF

## L4: Fundamental Laws (Complete)
- 3-op-amp gain theorem
- CMRR resistor mismatch (Pallás-Areny & Webster, 1991)
- Johnson-Nyquist theorem (1928)
- Schottky shot noise theorem (1918)
- 1/f noise empirical law (Hooge, 1969)
- Friis noise figure cascade formula (1944)
- Wheatstone bridge transfer function
- Noise figure definition
- Gain error propagation law
- CMRR impedance imbalance law

## L5: Algorithms/Methods (Complete)
- Two-point calibration (gain+offset)
- Multi-point linear regression calibration
- Quadratic polynomial calibration (Cramer's rule)
- Auto-zero offset measurement procedure
- Chopper stabilization ripple model
- Digital gain trim factor computation
- E96/E24 standard value selection
- Auto-range gain selection algorithm
- Ratiometric correction algorithm
- Moving average for cal readings
- Calibration standard deviation
- Self-test diagnostic procedure
- Temperature compensation (polynomial drift)
- Calibration validation check

## L6: Canonical Problems (Complete)
- Strain gauge bridge readout
- Thermocouple signal conditioning with CJC
- RTD (Pt100) resistance measurement
- Load cell weight measurement
- ECG biopotential front-end design
- EEG low-noise front-end
- EMG muscle signal conditioning
- Current shunt measurement (high-side/low-side)
- Hall effect sensor interface
- Ratiometric bridge measurement

## L7: Applications (Complete)
- Medical instrumentation: ECG front-end per IEC 60601-2-25
- Industrial weighing: load cell signal conditioning
- Temperature measurement: thermocouple + RTD
- Automotive: current shunt sensing
- Industrial: pressure bridge sensor
- Scientific: EEG/EMG biopotential acquisition

## L8: Advanced Topics (Partial, 2/5)
- [C] Chopper-stabilized IA modeling
- [C] Current-mode IA for low-voltage operation
- [D] MEMS sensor integration with IA
- [D] Wireless sensor node analog front-end
- [D] Multi-channel simultaneous sampling IAs

## L9: Research Frontiers (Partial, documented)
- Nanowatt instrumentation amplifiers for energy harvesting
- Flexible/stretchable electronics IA design
- Sub-threshold operation IA for biomedical implants
- Quantum sensor readout with cryogenic IAs
- AI-assisted auto-calibration and drift compensation