# mini-sensor-measurement

**Sensor Measurement & Instrumentation** — from-scratch, zero-dependency C implementations.

Covers the full sensor signal chain: transduction principles, Wheatstone bridge circuits,
signal conditioning, ADC/DAC data conversion, instrumentation amplifiers, MEMS sensors,
measurement error analysis (GUM), biosensor electrochemistry, and smart sensor systems.

Based on Fraden (*Handbook of Modern Sensors*), Pallas-Areny & Webster (*Sensors and Signal
Conditioning*), Doebelin (*Measurement Systems*), and IEEE/JCGM measurement standards.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-adc-dac-converter](mini-adc-dac-converter/) | Nyquist sampling, SAR ADC, ΔΣ modulator, quantization theory, dynamic metrics (SNR/SFDR/ENOB) | MIT 6.301, Stanford EE315, Berkeley EE140 |
| [mini-biosensor-chemical](mini-biosensor-chemical/) | Electrochemical transduction, Nernst equation, amperometry, potentiometry, glucose biosensor | MIT 6.555, UC Berkeley BioE 121 |
| [mini-instrumentation-amp](mini-instrumentation-amp/) | IA topologies (3-op-amp, 2-op-amp), gain structures, CMRR, bridge sensor interface | Stanford EE247, Berkeley EE105, Michigan EECS 411 |
| [mini-measurement-error](mini-measurement-error/) | Error types (bias/random/hysteresis), GUM uncertainty propagation, calibration regression, statistics | MIT 6.630, ETH 227-0455, Berkeley EE117 |
| [mini-mems-sensor](mini-mems-sensor/) | MEMS dynamics (spring-mass-damper), Allan variance, IMU sensor fusion (Mahony/Madgwick/Kalman) | MIT 6.630, Stanford EE247, MIT 16.485 |
| [mini-sensor-principle](mini-sensor-principle/) | Sensor definitions, Wheatstone bridge analysis, transfer functions, Johnson/Nyquist noise | MIT 6.630, Stanford EE247, Berkeley EE117 |
| [mini-signal-conditioning](mini-signal-conditioning/) | Amplification, filtering, linearization, cold-junction compensation, 4–20 mA current loop, isolation | MIT 6.003, Stanford EE247, Berkeley EE105 |
| [mini-smart-sensor](mini-smart-sensor/) | Energy harvesting, Allan variance analysis, edge ML feature extraction (RMS, crest factor, zero-crossing) | Stanford EE267, ETH 227-0436, MIT 6.630 |

## Design Philosophy

1. **From signal to digit** — Every module traces the physical-to-digital path: sensor element → signal conditioning → ADC → digital processing.
2. **Standards-based** — IEEE Std 1241 (ADC), IEEE Std 1139/952 (noise/Allan), JCGM GUM (uncertainty), ISO 5725 (accuracy).
3. **Zero dependencies** — Pure C99/C11 with `libc` and `libm` only; no external libraries, no vendor SDKs.
4. **Nine-layer knowledge** — L1 Definitions → L2 Concepts → L3 Math → L4 Laws → L5 Algorithms → L6 Canonical Problems → L7 Applications → L8 Advanced → L9 Frontiers.

## Building

```sh
# Build a specific sub-module
cd mini-adc-dac-converter && make

# Build all sub-modules
for d in mini-*/; do (cd "$d" && make); done

# Run tests
for d in mini-*/; do (cd "$d" && make test); done
```

## Project Structure

```
mini-sensor-measurement/
├── mini-adc-dac-converter/        # Nyquist sampling, SAR ADC, ΔΣ modulator, quantization, ADC metrics
├── mini-biosensor-chemical/       # Electrochemical transduction, Nernst, Butler-Volmer, glucose biosensor
├── mini-instrumentation-amp/      # 3-op-amp/2-op-amp topologies, gain-stage analysis, bridge interfacing
├── mini-measurement-error/        # Error taxonomy, Type-A/Type-B uncertainty, calibration, statistics, filtering
├── mini-mems-sensor/              # MEMS dynamics, noise processes, 6-position calibration, IMU sensor fusion
├── mini-sensor-principle/         # Sensor definitions, Wheatstone bridge, transfer functions, noise models
├── mini-signal-conditioning/      # Amplification stages, filter design, linearization, CJC, 4–20 mA loop
├── mini-smart-sensor/             # Energy harvesting, Allan variance fitting, edge ML feature extraction
├── .gitignore                     # Build artifact and IDE exclusions
├── README.md                      # This file (English)
└── README-CN.md                   # Chinese version
```

## License

MIT — See individual sub-module directories for details.
