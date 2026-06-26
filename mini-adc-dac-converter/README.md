# mini-adc-dac-converter

ADC/DAC Converter Theory and Implementation Library.

**Codebase**: 6513 lines (include/ + src/), 6 headers, 7 source files.

---

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Complete (3 applications: sensor data, audio ADC, instrumentation)
- **L8**: Complete (4 advanced topics: time-interleaved, dithering, sub-radix-2 SAR, MASH)
- **L9**: Partial (documented, not implemented)

**Self-Check Score**: 17/18 (threshold: ≥16)

---

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | ✅ Complete | 20+ struct typedefs: ADC/DAC specs, quantization, sampling, ΣΔ, SAR |
| **L2** | Core Concepts | ✅ Complete | Flash, SAR, pipeline, ΣΔ, R-2R, string, current-steering architectures |
| **L3** | Math Structures | ✅ Complete | FFT, Chebyshev model, Bennett noise, Pelgrom mismatch, Poisson summation |
| **L4** | Fundamental Laws | ✅ Complete | SNR=6.02N+1.76, q²/12, kT/C, Lee's criterion, Nyquist, OSR gain |
| **L5** | Algorithms | ✅ Complete | 21 algorithms: SAR search, Lloyd-Max, CIC, FFT, dithering, calibration |
| **L6** | Canonical Problems | ✅ Complete | 10 problems: ADC SNR sweep, SAR sim, ΣΔ design, glitch, interpolation |
| **L7** | Applications | ✅ Complete | Sensor acquisition, audio ΣΔ ADC, instrumentation-grade ADC |
| **L8** | Advanced Topics | ✅ Complete | TPDF dithering, sub-radix-2 SAR, MASH 1-1, time-interleaved |
| **L9** | Research Frontiers | ⚡ Partial | AI calibration, compressive sensing, quantum ADC (documented) |

---

## Core Definitions (L1)

- **ADC Specification**: resolution (bits), sample rate (f_s), reference voltage, LSB, aperture jitter, conversion time, coding scheme
- **DAC Specification**: resolution, update rate, settling time, glitch energy, output impedance
- **Performance Metrics**: SNR, SINAD, ENOB, THD, SFDR, noise floor, DNL, INL
- **Quantization**: step size, uniform/non-uniform, mid-tread/mid-rise, μ-law, A-law
- **Sampling**: Nyquist frequency, oversampling ratio, anti-alias filter, coherent sampling
- **ΣΔ Modulation**: modulator order, OSR, NTF, STF, quantizer levels
- **SAR ADC**: CDAC, comparator offset, bit-cycling, redundancy, mismatch

## Core Theorems (L4)

| Theorem | Formula | Implementation |
|---------|---------|----------------|
| Ideal ADC SNR | SNR = 6.02·N + 1.76 dB | `adc_ideal_snr_db()` |
| Quantization noise power | σ² = q²/12 | `quantization_noise_power()` |
| Nyquist-Shannon | f_s > 2·BW | `nyquist_check()` |
| Oversampling SNR gain | ΔSNR = 10·log₁₀(OSR) | `oversampling_snr_gain()` |
| ΔΣ in-band noise | IBN ≈ σ_q²·π^{2L} / ((2L+1)·OSR^{2L+1}) | `sdm_inband_noise_power()` |
| kT/C thermal noise | v_n = √(kT/C) | `ktc_noise_rms()` |
| Lee's stability | max\|NTF\| < 2.0 | `sdm_check_lee_criterion()` |
| Jitter SNR limit | SNR = -20·log₁₀(2π·f·t_j) | `adc_jitter_snr_limit()` |
| Pelgrom's law | σ(ΔC/C) = A_C/√(W·L) | `capacitor_with_mismatch()` |
| Sinc interpolation | x(t) = Σx[n]·sinc((t-nT)/T) | `sinc_interpolate()` |

## Core Algorithms (L5)

- **SAR binary search**: N successive bit trials, MSB first → `sar_bit_trial()`
- **Sub-radix-2 SAR**: Redundant search with error correction → `sar_redundant_conversion()`
- **Lloyd-Max quantizer**: Iterative centroid/nearest-neighbor → `lloyd_max_quantizer_design()`
- **CIC decimation**: Cascaded integrator-comb filter → `cic_decimate()`
- **FIR rate conversion**: Polyphase decimation/interpolation → `fir_decimate()`, `fir_interpolate()`
- **Noise shaping**: 1st/2nd order error feedback → `noise_shape_first_order()`, `noise_shape_second_order()`
- **ΔΣ MASH 1-1**: Cascaded noise shaping → `sdm_mash_1_1_process()`
- **Radix-2 FFT**: Cooley-Tukey DIT → `fft_real_to_complex()` (internal)
- **Dynamic testing**: SINAD/SNR/THD/SFDR/ENOB → `adc_full_dynamic_test()`
- **Window functions**: Hann, Hamming, Blackman, Blackman-Harris → `adc_apply_window()`
- **TPDF dithering**: Non-subtractive triangular PDF → `dither_nonsubtractive_tpdf()`
- **μ-law / A-law**: Telephony companding per ITU-T G.711 → `mu_law_compress()`, `a_law_compress()`

## Classic Problems (L6)

1. **SNR vs Resolution** — 6 to 24 bits sweep with jitter floor
2. **SAR ADC 10-bit Simulation** — kT/C noise, offset, mismatch
3. **ΔΣ Audio ADC Design** — 1st/2nd/3rd order, OSR 32-256
4. **DAC Glitch Energy** — MSB transition timing skew
5. **R-2R Linearity** — Resistor mismatch impact
6. **Oversampling Benefit** — OSR 2-256, ΔENOB
7. **Band-Pass Sampling** — 70-90 MHz IF band
8. **Sinc Reconstruction** — Interpolation of sampled data
9. **IEEE 1241 Test Suite** — Full dynamic characterization
10. **Two-Tone IMD** — IM2/IM3 measurement

## Nine-School Course Mapping

| School | Course | Alignment |
|--------|--------|-----------|
| MIT | 6.301, 6.775 | SAR, ΣΔ, CMOS data converters |
| Stanford | EE315, EE315B | Full data converter theory |
| Berkeley | EE247, EE140 | ΣΔ modulators, analog IC |
| Illinois | ECE 483 | Comparator, CDAC design |
| Michigan | EECS 511 | Converter fundamentals |
| Georgia Tech | ECE 6416 | High-speed ADC |
| TU Munich | — | RF sampling, bandpass |
| ETH | 227-0455 | — |
| Tsinghua | 信号与系统, 模拟集成电路 | Sampling, quantization, data converters |

## Quick Start

```bash
make          # Build library and run all tests
make examples # Build example programs
make clean    # Remove build artifacts
```

## File Structure

```
mini-adc-dac-converter/
├── Makefile                    # Build system
├── README.md                   # This file
├── include/                    # Headers (6 files, 2361 lines)
│   ├── adc_dac_core.h         # Core ADC/DAC types, SNR, calibration
│   ├── adc_metrics.h          # SINAD, THD, SFDR, ENOB, IMD
│   ├── quantization.h         # Uniform/non-uniform, Lloyd-Max, μ-law
│   ├── sampling.h             # Nyquist, aliasing, decimation, interpolation
│   ├── sigma_delta.h          # ΔΣ modulation, NTF, MASH
│   └── sar_adc.h              # CDAC, binary search, kT/C, mismatch
├── src/                        # Implementations (7 files, 4152 lines)
│   ├── adc_core.c             # ADC ideal models, DNL/INL, calibration
│   ├── dac_core.c             # R-2R, binary-weighted, string, segmented
│   ├── adc_metrics.c          # FFT-based dynamic testing (IEEE 1241)
│   ├── quantization.c         # Lloyd-Max, dithering, noise shaping
│   ├── sampling.c             # CIC/FIR decimation, sinc interpolation
│   ├── sigma_delta.c          # ΔΣ modulator, MASH, NTF design
│   └── sar_adc.c              # CDAC model, SAR search, mismatch
├── tests/                      # assert-based tests (7 files)
├── examples/                   # End-to-end examples (4 files)
├── demos/                      # Visualization/plot scripts
├── benches/                    # Performance benchmarks
└── docs/                       # Knowledge documentation (5 files)
    ├── knowledge-graph.md
    ├── coverage-report.md
    ├── gap-report.md
    ├── course-alignment.md
    └── course-tree.md
```

## References

- Maloberti, "Data Converters" (2007)
- Schreier & Temes, "Understanding Delta-Sigma Data Converters" (2005)
- Kester, "Data Conversion Handbook" (2005)
- IEEE Std 1241-2010: ADC Test Methods
- IEEE Std 1658-2011: DAC Test Methods
- Oppenheim & Schafer, "Discrete-Time Signal Processing" (2010)
- McCreary & Gray, "All-MOS Charge Redistribution ADC" JSSC (1975)
- Hogenauer, "An Economical Class of Digital Filters" IEEE ASSP (1981)
- Pelgrom et al., "Matching Properties of MOS Transistors" JSSC (1989)
