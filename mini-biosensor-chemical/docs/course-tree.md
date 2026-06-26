# Course Tree — mini-biosensor-chemical

## Prerequisite Dependency Graph

```
                           ┌─────────────────────────────────┐
                           │   mini-biosensor-chemical       │
                           │   (This Module)                 │
                           └───────────────┬─────────────────┘
                                           │
              ┌────────────────────────────┼────────────────────────────┐
              │                            │                            │
    ┌─────────▼──────────┐    ┌────────────▼───────────┐    ┌──────────▼──────────┐
    │ mini-signal-       │    │ mini-circuit-           │    │ mini-analog-         │
    │ system-theory      │    │ analysis                │    │ electronics          │
    │ (Fourier, Laplace, │    │ (Impedance, Thevenin,   │    │ (Op-amps, filters,   │
    │  convolution)      │    │  Kirchhoff)             │    │  noise analysis)     │
    └─────────┬──────────┘    └────────────┬───────────┘    └──────────┬──────────┘
              │                            │                            │
              └────────────────────────────┼────────────────────────────┘
                                           │
                              ┌────────────▼───────────┐
                              │ mini-digital-          │
                              │ electronics            │
                              │ (ADC, DAC, sampling)   │
                              └────────────┬───────────┘
                                           │
                              ┌────────────▼───────────┐
                              │ mini-digital-signal-   │
                              │ process                │
                              │ (FIR/IIR filters, FFT) │
                              └────────────────────────┘
```

## Internal Dependency Graph (This Module)

```
L1: Types & Definitions
│
├──► L2: Core Concepts (electrochemical, optical, piezoelectric)
│    │
│    ├──► L3: Math Structures (MM, Hill, Langmuir, Nernst, BV, BL)
│    │    │
│    │    ├──► L4: Fundamental Laws (with proofs)
│    │    │    │
│    │    │    ├──► L5: Algorithms (calibration, signal processing, kinetics)
│    │    │    │    │
│    │    │    │    ├──► L6: Canonical Problems (glucose, ELISA, CV, ISE)
│    │    │    │    │    │
│    │    │    │    │    ├──► L7: Applications (PoC, e-nose, DNA microarray)
│    │    │    │    │    │    │
│    │    │    │    │    │    ├──► L8: Advanced Topics (PDS, transfer)
│    │    │    │    │    │    │    │
│    │    │    │    │    │    │    └──► L9: Research Frontiers (documented)
│    │    │    │    │    │    │
│    │    │    │    │    │    └──► docs/knowledge-graph.md
│    │    │    │    │    │
│    │    │    │    │    └──► examples/ (3 canonical examples)
│    │    │    │    │
│    │    │    │    └──► tests/ (assert-based test suite)
│    │    │    │
│    │    │    └──► src/biosensor_lean.lean (formal proofs)
│    │    │
│    │    └──► src/biosensor_kinetics.c, biosensor_electrochemical.c, biosensor_optical.c
│    │
│    └──► include/ + src/ for each transduction modality
│
└──► include/biosensor_types.h → src/biosensor_types.c
```

## Knowledge Paths

### Path 1: Electrochemical Biosensors
L1 Types → L2 Electrochemical → L3 Nernst/BV/Cottrell → L4 Fundamental Laws → L5 Calibration → L6 Glucose Biosensor → L7 PoC Lactate

### Path 2: Optical Biosensors
L1 Types → L2 Optical → L3 Beer-Lambert/FRET/Fresnel → L4 Laws → L5 4PL Fit → L6 ELISA → L7 DNA Microarray

### Path 3: Biorecognition Kinetics
L1 Types → L3 MM/Hill/Langmuir → L4 Enzyme Kinetics → L5 Parameter Estimation → L6 Inhibition Analysis → L7 Immunoassay

### Path 4: Sensor Array / E-Nose
L1 Types → L2 Chemiresistor → L5 Signal Processing → L6 CV/Voltammetry → L7 Pattern Classification → L8 Calibration Transfer

## Cross-Module Dependencies

| This Module Uses | From Module | Knowledge Item |
|-----------------|-------------|----------------|
| OLS regression | mini-signal-system-theory | Linear algebra, matrix inversion |
| Convolution/SG filter | mini-digital-signal-process | FIR filter design |
| SNR, noise analysis | mini-analog-electronics | Thermal noise, shot noise |
| Op-amp transimpedance | mini-analog-electronics | Current-to-voltage conversion |
| ADC quantization | mini-digital-electronics | Bit depth, ENOB |
| Impedance concepts | mini-circuit-analysis | Complex impedance, phasors |
