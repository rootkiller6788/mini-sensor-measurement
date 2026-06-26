# Course Tree — mini-sensor-principle

## Prerequisites Dependency Tree

```
mini-sensor-principle
├── mini-signal-system-theory (0)
│   └── L3: Fourier/Laplace transforms (frequency response)
│   └── L3: Convolution (filtering)
│   └── L4: Nyquist sampling theorem (anti-alias filter)
│
├── mini-circuit-analysis (1)
│   └── L2: Kirchhoff's laws (bridge circuit analysis)
│   └── L2: Voltage dividers (bridge output)
│   └── L3: Thevenin equivalent (bridge output impedance)
│
├── mini-analog-electronics (2)
│   └── L2: Operational amplifiers (instrumentation amp, gain)
│   └── L5: Filter design (anti-alias, noise filtering)
│   └── L2: Current sources (RTD/thermistor excitation)
│
├── mini-digital-electronics (3)
│   └── L2: ADC quantization noise
│   └── L2: ENOB calculation
│
├── mini-digital-signal-process (6)
│   └── L5: Digital filters (moving avg, exponential, median)
│   └── L5: FFT-based noise PSD analysis
│   └── L3: Allan variance computation
│
└── mini-electromagnetic-wave (7)
    └── L4: Johnson-Nyquist (thermal noise) fluctuation-dissipation
    └── L4: Shot noise (quantum nature of charge)
    └── L2: Seebeck effect (thermocouple)
```

## Internal Dependency Graph

```
sensor_defs.h ─────────────────────────────────────────────┐
    │ (L1: Types, structs, enums)                          │
    ├── sensor_transfer.h/c                                │
    │   ├── poly_transfer (L2: polynomial model)           │
    │   ├── thermistor (L2+L4: Steinhart-Hart)             │
    │   ├── rtd_cvd (L2+L4: Callendar-Van Dusen)           │
    │   ├── tc_model (L2+L4: Seebeck + ITS-90)             │
    │   ├── sensor_1st_order (L3: ODE solutions)           │
    │   ├── sensor_2nd_order (L3: damped response)         │
    │   ├── sensor_lut (L2: interpolation table)           │
    │   └── linearity/hysteresis (L3: fit analysis)        │
    │                                                       │
    ├── sensor_noise.h/c                                   │
    │   ├── johnson_noise (L4: Nyquist 1928)               │
    │   ├── shot_noise (L4: Schottky 1918)                 │
    │   ├── flicker_noise (L3: 1/f PSD integration)        │
    │   ├── allan_variance (L3: overlapping estimator)     │
    │   ├── nep/d_star (L3: NEP, D*)                       │
    │   └── noise_factor_cascade (L4: Friis 1944)          │
    │                                                       │
    ├── sensor_bridge.h/c                                  │
    │   ├── bridge_output (L2+L3: exact + approx)          │
    │   ├── bridge_nonlinearity (L3: error analysis)       │
    │   ├── lead_wire_error (L2: 2-wire effect)            │
    │   ├── strain calculation (L2: GF → ε)                │
    │   └── amplifier requirements (L3: CMRR)              │
    │                                                       │
    ├── sensor_calibration.h/c                             │
    │   ├── linear LS (L5: normal equations)               │
    │   ├── polynomial LS (L5: Vandermonde + Gauss elim)   │
    │   ├── weighted LS (L5: w=1/σ²)                       │
    │   ├── temp compensation (L5: α, β)                   │
    │   └── cross-sensitivity (L5: matrix inversion)       │
    │                                                       │
    ├── sensor_signal_cond.h/c                             │
    │   ├── amplifier gains (L2: op-amp configs)           │
    │   ├── anti-alias filter (L5: order selection)        │
    │   ├── 4-20mA current loop (L2: industrial standard)  │
    │   ├── digital filters (L5: MA/IIR/median)            │
    │   └── excitation (L2: V/I/AC sources)                │
    │                                                       │
    └── sensor_formal.lean                                 │
        ├── Johnson-Nyquist (L4)                           │
        ├── Wheatstone Bridge Balance (L4)                 │
        ├── First-Order Response (L4)                      │
        ├── Steinhart-Hart Monotonicity (L4)               │
        ├── Noise Additivity (L4)                          │
        ├── Sensitivity (L4)                               │
        └── RTD Monotonicity (L4)                          │
```

## Knowledge Level Progression

```
L1: sensor_defs ──────────────────────┐
    (types, structs, enums)            │
                                       │
L2: sensor_transfer ──────────────────┤
    sensor_bridge ────────────────────┤
    sensor_signal_cond ───────────────┤
    (core concepts, transfer functions)│
                                       │
L3: sensor_transfer (dynamic) ────────┤
    sensor_noise (PSD, Allan) ────────┤
    sensor_bridge (analysis) ─────────┤
    (mathematical structures)          │
                                       │
L4: sensor_noise (Johnson, shot) ─────┤
    sensor_transfer (Seebeck, CVD) ───┤
    sensor_formal.lean ───────────────┤
    (fundamental laws + verification)  │
                                       │
L5: sensor_calibration ───────────────┤
    sensor_signal_cond (filtering) ───┤
    (algorithms & methods)             │
                                       │
L6: examples/*.c ─────────────────────┘
    (canonical problems)
```

## Research Frontiers (L9) Growth Path

```
Current L7 (Applications)
  ├── Industrial Temperature → L9: Smart factory sensors (IIoT)
  ├── Structural Strain Gauging → L9: Self-powered wireless strain nodes
  └── MEMS Inertial → L9: Chip-scale atomic sensors

Current L8 (Advanced)
  ├── Noise Cascade → L9: Quantum-limited amplifier chains
  └── Allan Decomposition → L9: AI-driven noise source identification
```
