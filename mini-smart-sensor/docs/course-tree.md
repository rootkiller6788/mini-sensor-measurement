# Course Dependency Tree — mini-smart-sensor

Prerequisite knowledge tree for this module. Shows what foundational
knowledge is required before studying smart sensor theory.

## Prerequisite Tree

```
                    ┌─────────────────────────────┐
                    │  mini-smart-sensor           │
                    │  (This Module)               │
                    └─────────────┬───────────────┘
                                  │
        ┌─────────────┬───────────┼───────────┬─────────────┐
        │             │           │           │             │
        ▼             ▼           ▼           ▼             ▼
┌──────────────┐ ┌──────────┐ ┌─────────┐ ┌─────────┐ ┌──────────┐
│ sensor       │ │ signal   │ │ analog  │ │ digital │ │ mcu-     │
│ principle    │ │ cond.    │ │ electro │ │ electro │ │ embedded │
│ (module 8)   │ │ (module 8)│ │ (mod 2) │ │ (mod 3) │ │ (mod 4)  │
└──────┬───────┘ └────┬─────┘ └────┬────┘ └────┬────┘ └────┬─────┘
       │               │            │           │            │
       ▼               ▼            ▼           ▼            ▼
┌─────────────────────────────────────────────────────────────────┐
│                  mini-signal-system-theory (module 0)           │
│  Fourier analysis, Laplace/Z transform, filtering theory       │
└─────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                  mini-circuit-analysis (module 1)               │
│  DC/AC circuits, Ohm's law, Kirchhoff's laws, Thevenin/Norton  │
└─────────────────────────────────────────────────────────────────┘
```

## Topic Dependencies

### L1: Definitions
```
Prerequisites: None (self-contained terminology)
Dependencies within module: None
```

### L2: Core Concepts
```
Prerequisites: L1 (definitions)
External dependencies: 
  - mini-circuit-analysis: Ohm's law, voltage dividers
  - mini-sensor-principle: transduction principles
  - mini-analog-electronics: amplifier basics
```

### L3: Mathematical Structures
```
Prerequisites: L1, L2
External dependencies:
  - mini-signal-system-theory: Fourier/Laplace, frequency response
  - mini-circuit-analysis: Wheatstone bridge
  - mini-digital-electronics: ADC theory
```

### L4: Fundamental Laws
```
Prerequisites: L3
External dependencies:
  - mini-signal-system-theory: SNR, sampling theorem
  - mini-circuit-analysis: Ohm's law, bridge balance
  - mini-communication-principle: Shannon-Hartley, information theory
```

### L5: Algorithms/Methods
```
Prerequisites: L1-L4
External dependencies:
  - mini-signal-system-theory: FIR/IIR filter theory
  - mini-digital-signal-process: FFT, spectral analysis
  - mini-control-automation: Kalman filter, state-space
```

### L6: Canonical Problems
```
Prerequisites: L1-L5
External dependencies: All of the above
```

### L7: Applications
```
Prerequisites: L1-L6
External dependencies:
  - mini-iot-edge-computing: IoT architectures
  - mini-wireless-mobile-comm: wireless interfaces
  - mini-industrial-fieldbus: industrial protocols (CAN, Modbus, HART)
```

### L8: Advanced Topics
```
Prerequisites: L1-L7
External dependencies:
  - mini-power-electronics: energy harvesting
  - mini-radar-remote-sensing: signal processing for sensors
```

### L9: Research Frontiers
```
Prerequisites: L1-L8
External dependencies: Current research literature
```

## Learning Pathway

### Recommended Study Order

1. **mini-circuit-analysis** (module 1) — DC/AC fundamentals
2. **mini-signal-system-theory** (module 0) — Fourier, filtering
3. **mini-analog-electronics** (module 2) — Amplifiers, op-amps
4. **mini-digital-electronics** (module 3) — ADC/DAC
5. **mini-sensor-principle** (module 8) — Transduction physics
6. **mini-signal-conditioning** (module 8) — Bridge circuits, filtering
7. **mini-smart-sensor** (this module) — Integration of all above

### Time Estimate

| Phase | Topic | Estimated Hours |
|-------|-------|----------------|
| Foundation | Modules 0-3 | 40-60 hours |
| Sensor Basics | Sensor principle + signal conditioning | 15-20 hours |
| Smart Sensor | This module | 20-30 hours |
| Advanced | L7-L9 topics | 10-15 hours |
| **Total** | | **85-125 hours** |
