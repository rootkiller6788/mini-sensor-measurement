# mini-mems-sensor — MEMS Sensors

## Module Status: COMPLETE &#x2705;

| Level | Name | Status |
|-------|------|--------|
| L1 | Definitions | **Complete** — 13 core definitions (12 sensor types, 8 transduction mech, structs) |
| L2 | Core Concepts | **Complete** — 12 concepts (SMD dynamics, damping mechanisms, capacitive readout) |
| L3 | Math Structures | **Complete** — 11 structures (ODE, PSD/Welch, Allan var, quaternion, rotation matrix) |
| L4 | Fundamental Laws | **Complete** — 10 laws (Newton, Hooke, Coriolis, thermo-mech noise, Stoney, Johnson) |
| L5 | Algorithms/Methods | **Complete** — 14 algorithms (6-pos calib, rate table, Madgwick, Mahony, Kalman, Allan fit) |
| L6 | Canonical Problems | **Complete** — 7 problems (AHRS, roll/pitch, tilt compass, DR, full calib pipeline) |
| L7 | Applications | **Complete** — 6 apps (iPhone orientation, Quadrotor, Airbag, ESC, Pedometer, ISO 10816) |
| L8 | Advanced Topics | **Complete** — 4 topics (GPS-denied ZUPT, stationary detect, MEMS mirror, noise propagation) |
| L9 | Research Frontiers | **Partial** — Standard Quantum Limit (documented + Lean theorem) |

## Core Definitions (L1)

- **MemsSensorType**: 12 MEMS sensor types (accelerometer, gyroscope, magnetometer, pressure, microphone, resonator, temperature, humidity, gas, flow, force, tilt)
- **MemsTransduceType**: 8 transduction mechanisms (capacitive, piezoresistive, piezoelectric, thermal, optical, tunneling, resonant, magnetic)
- **MemsSpec**: Full sensor datasheet spec (FSR, sensitivity, SF error, bias, noise density, bandwidth, Q, cross-axis, g-sensitivity, ODR, etc.)
- **MemsVector3**: 3-axis vector with complete math library (add, sub, scale, dot, cross, normalize, magnitude)
- **MemsQuaternion**: Orientation quaternion with multiply, conjugate, normalize, rotate
- **MemsMatrix33**: 3×3 rotation/alignment matrix
- **MemsSeismicMass**: Proof mass mechanical parameters (mass, spring k, damping c, resonant freq, Q, ζ)
- **MemsSpringMassDamper**: Full state for 2nd-order ODE simulation (displacement, velocity, acceleration)

## Core Theorems (L4)

- **Newton's 2nd Law**: F = ma → a = (F_ext - c·ẋ - k·x) / m  (drives `mems_smd_step()`)
- **Hooke's Law**: F = kx → x_static = m·a/k  (static MEMS accelerometer sensitivity)
- **Coriolis Force**: F_c = 2·m·v_drive·Ω_input  (MEMS gyroscope principle)
- **Thermo-mechanical Noise Limit**: NEA = sqrt(4·k_B·T·ω_n/(m·Q))  (Gabrielson 1993)
- **Nyquist-Johnson Noise**: v_n = sqrt(4·k_B·T·R·Δf)  (Johnson 1928, Nyquist 1928)
- **Shot Noise**: i_n = sqrt(2·q_e·I_DC)  (Schottky 1918)
- **Stoney's Formula**: σ_f = (E_s·t_s²)/(6·(1-ν_s)·t_f·R)  (Stoney 1909)
- **Resonance Gain**: |H(jω_n)| = 1/(2ζ)  (2nd-order system response)
- **Barometric Formula**: h = (T₀/L)·(1 - (P/P₀)^(R·L/g))  (ISA standard)
- **Zener's TED**: Q⁻¹ = (E·α²·T₀/(ρ·C_p))·(ω·τ/(1+ω²·τ²))  (thermoelastic damping)

## Core Algorithms (L5)

- **6-Position Accelerometer Calibration**: Least-squares bias + scale factor from ±g per axis
- **Rate Table Gyro Calibration**: Linear regression per axis from multiple rotation rates
- **Magnetometer Ellipsoid Fitting**: Hard-iron center + soft-iron scaling from min/max
- **Least-Squares Axis Alignment**: S = (V_meas·V_ref^T)·(V_ref·V_ref^T)⁻¹ (Gebre-Egziabher 2006)
- **Polynomial Temp Compensation**: Up to 4th-order Gaussian elimination (IEEE 2700-2017)
- **Complementary Filter** (SO(3)): q = α·q_gyro + (1-α)·q_accel (Euston 2008)
- **Mahony Filter**: Nonlinear complementary with PI correction, proven stability (Mahony 2008)
- **Madgwick AHRS**: Gradient descent with 6-axis IMU + 9-axis MARG (Madgwick 2010)
- **1D Kalman Filter**: Predict-update cycle for angle estimation
- **Allan Variance Analysis**: Overlapping estimator with 5-slope noise identification (IEEE 1139)
- **Welch PSD Estimation**: Averaged periodogram with Hann window, radix-2 FFT
- **Box-Muller Gaussian Generator**: For noise simulation
- **1/f Flicker Noise Generator**: Kasdin (1995) multi-pole approximation

## Canonical Problems (L6)

- **AHRS Orientation Estimation**: Full gyro+accel+mag fusion (Mahony/Madgwick/Complementary)
- **Roll/Pitch from Accelerometer**: Static tilt computation for smartphones
- **Tilt-Compensated Compass**: Magnetometer yaw with accel-based tilt removal
- **IMU Dead Reckoning**: Acceleration double integration with gravity removal
- **IMU Calibration Pipeline**: End-to-end 6-position + rate table + ellipsoid fit
- **Tilt Sensing Demo**: Phone orientation scenarios with temperature error budget
- **Sensor Fusion Comparison**: Side-by-side comp/Mahony/Madgwick on synthetic data

## Course Mapping

| School | Course | Key MEMS Topics |
|--------|--------|----------------|
| MIT | 6.630 EM Waves & Sensors | Capacitive sensing, noise limits |
| Stanford | EE247 MEMS & Sensors | SMD, damping, pull-in, TED, Coriolis gyro |
| Berkeley | EE117 EM / EE123 DSP / EE245 MEMS | MEMS design, noise analysis, Allan var |
| Illinois | ECE 310 DSP / ECE 451 EM | Sensor signal processing, PSD |
| Michigan | EECS 351 DSP / EECS 411 / ME 552 | Resonators, RF MEMS, thin-film stress |
| Georgia Tech | ECE 4270 DSP / ECE 6350 EM / ECE 6450 | Sensor fusion, gyro dynamics |
| TU Munich | Signal Processing / High-Frequency | MEMS oscillators, structured light |
| ETH | 227-0427 DSP / 227-0455 EM / 151-0172 MEMS | Kalman, Coriolis, Stoney, TED |
| Tsinghua | 传感器原理 / 惯性导航 / 汽车电子 | Calibration, AHRS, ESC, Airbag |

## File Structure

```
mini-mems-sensor/
├── Makefile                 # make / make test / make examples
├── README.md                # This file
├── include/                 (6 headers)
│   ├── mems_basics.h        — L1-L2: Types, vector, quaternion, spec
│   ├── mems_dynamics.h      — L2-L3: SMD, frequency response, damping, Coriolis
│   ├── mems_noise.h         — L3-L4: Thermo-mech noise, Allan var, PSD, Welch
│   ├── mems_calibration.h   — L5: 6-pos, rate table, mag ellipsoid, temp comp
│   ├── mems_fusion.h        — L5-L6: Comp filter, Mahony, Madgwick, Kalman
│   └── mems_applications.h  — L7-L8: Phone IMU, drone, airbag, ESC, pedometer
├── src/                     (6 C + 1 Lean)
│   ├── mems_basics.c        — Core types, vector math, unit conversions, spec init
│   ├── mems_dynamics.c      — SMD solver, harmonic response, damping models
│   ├── mems_noise.c         — Noise limits, Allan analysis, Welch PSD, generators
│   ├── mems_calibration.c   — LS calibration, alignment matrix, temp comp
│   ├── mems_fusion.c        — Quaternion ops, 4 fusion filters, dead reckoning
│   ├── mems_applications.c  — Phone/drone/airbag/ESC/pedometer/INS/mirror
│   └── mems_sensor.lean     — 25 formal theorems (L4 laws in Nat, Lean 4)
├── tests/                   (4 test suites)
│   ├── test_basics.c        — 10 test groups: types, vectors, quat, specs
│   ├── test_dynamics.c      — 13 test groups: SMD, harmonics, damping, Coriolis
│   ├── test_noise.c         — 12 test groups: NEA, Johnson, shot, Allan, Welch
│   └── test_calibration.c   — 9 test groups: 6-pos, gyro, mag, temp comp
├── examples/                (3 end-to-end examples)
│   ├── example_imu_calib.c      — Full IMU calibration pipeline
│   ├── example_accel_tilt.c     — Tilt sensing with phone orientation demo
│   └── example_sensor_fusion.c  — Sensor fusion comparison + drift analysis
├── demos/
├── benches/
└── docs/
    ├── knowledge-graph.md
    ├── coverage-report.md
    ├── gap-report.md
    ├── course-alignment.md
    └── course-tree.md
```

## Build & Test

```bash
make          # Compile all source files
make test     # Build and run all tests
make examples # Build and run all examples
make clean    # Remove build artifacts
```

## References

- IEEE Std 1293-2018 — Inertial Sensor Terminology
- IEEE Std 2700-2017 — MEMS Performance Parameters
- IEEE Std 1139-2008 — Allan Variance (Frequency Stability)
- IEEE Std 952-2020 — Inertial Sensor Noise
- IEEE Std 1554-2005 — Inertial Sensor Test Equipment
- IEEE Std 1431-2004 — Coriolis Vibratory Gyroscopes
- ISO 10816-1:1995 — Vibration Severity
- ISO 12097-2 — Airbag Testing
- Bao (2005) — Analysis and Design Principles of MEMS Devices
- Senturia (2001) — Microsystem Design
- Acar & Shkel (2009) — MEMS Vibratory Gyroscopes
- Titterton & Weston (2004) — Strapdown Inertial Navigation Technology
- Groves (2013) — Principles of GNSS, Inertial, and Multisensor Navigation
- Markley & Crassidis (2014) — Fundamentals of Spacecraft Attitude
- Gabrielson (1993) — Mechanical-Thermal Noise in MEMS
- Mahony et al. (2008) — Nonlinear Complementary Filters on SO(3)
- Madgwick (2010) — An Efficient Orientation Filter for IMU and MARG
- Foxlin (2005) — Pedestrian Tracking with Shoe-Mounted Inertial Sensors
- Kasdin (1995) — Discrete Simulation of Colored Noise
- Zener (1937) — Internal Friction in Solids
