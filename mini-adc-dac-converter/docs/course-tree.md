# Course Tree — mini-adc-dac-converter

## Prerequisites

```
mini-adc-dac-converter
  ├── mini-signal-system-theory     (L3: Fourier/Laplace/Z transforms)
  │   └── FFT for spectral analysis of ADC data
  ├── mini-circuit-analysis          (L2: RC settling, Thévenin/Norton)
  │   └── R-2R ladder, CDAC charge redistribution
  ├── mini-analog-electronics        (L2: op-amp, comparator)
  │   └── Comparator design, sample-and-hold
  ├── mini-digital-electronics       (L2: registers, binary arithmetic)
  │   └── SAR logic, thermometer decoder
  ├── mini-digital-signal-process    (L3: sampling, multirate DSP)
  │   └── Decimation, interpolation, window functions
  └── mini-semiconductor-devices     (L8: MOS capacitor matching)
      └── Pelgrom's law, kT/C noise physics
```

## Dependent Modules

```
mini-adc-dac-converter
  ├──→ mini-sensor-measurement       (uses ADC theory for sensor interfacing)
  ├──→ mini-communication-principle  (uses ADC/DAC in transceiver design)
  ├──→ mini-audio-video-eng          (uses ΣΔ for audio, pipeline for video)
  └──→ mini-radar-remote-sensing     (uses high-speed ADC for IF sampling)
```

## Concept Dependencies

```
Sampling Theorem (L4)
  ├── Quantization Theory (L4)
  │   ├── DNL/INL (L5)
  │   │   └── SAR ADC Design (L6)
  │   ├── Oversampling (L4)
  │   │   └── ΣΔ Modulation (L5)
  │   │       └── ΣΔ ADC Design (L6)
  │   └── Dithering (L5)
  │       └── Audio ADC (L7)
  ├── Anti-alias Filtering (L5)
  │   └── Decimation (L5)
  └── Reconstruction (L5)
      └── DAC Design (L6)
```

## L9 Research Frontiers Dependencies

```
Compressive Sensing ADC
  └── Requires: L3 Fourier analysis, L4 Sampling theorem, L8 Random matrix theory

AI-Enhanced Calibration
  └── Requires: L5 SAR calibration, L8 Neural network approximation theory

Quantum DAC
  └── Requires: L4 Josephson effect, L9 Quantum metrology
```
