# Knowledge Graph — mini-adc-dac-converter

## L1: Definitions ✅ Complete
- ADC resolution (bits), sampling rate (f_s), reference voltage (V_ref)
- LSB = V_ref / 2^N, quantization step
- DAC settling time, glitch energy
- ENOB (Effective Number of Bits)
- SNR, SINAD, SFDR, THD
- DNL (Differential Non-Linearity), INL (Integral Non-Linearity)
- Quantization error: e = Q(x) - x
- Anti-aliasing filter, Nyquist frequency
- Aperture jitter, sample-and-hold
- DAC monotonicity
- Oversampling ratio (OSR)
- ΣΔ modulator order, quantizer levels
- CDAC (capacitive DAC), unit capacitor
- SAR register, binary search, bit-cycling
- Mid-tread vs mid-rise quantizer
- μ-law, A-law companding
- Thermometer code, Gray code, offset binary, two's complement

## L2: Core Concepts ✅ Complete
- Nyquist-Shannon sampling theorem
- Uniform quantization process
- Binary-weighted CDAC, R-2R ladder DAC
- Sigma-Delta noise shaping principle
- Pipeline ADC architecture (coarse + fine)
- SAR successive approximation algorithm
- Flash ADC (full-parallel)
- Dual-slope integrating ADC
- Oversampling and decimation
- Time-interleaved ADC
- Dithering (subtractive and non-subtractive)
- Charge redistribution in SAR ADC
- Current-steering DAC
- String DAC (Kelvin divider) — inherent monotonicity
- Segmented/hybrid DAC architecture
- Coherent vs non-coherent sampling for FFT testing
- Window functions (Hann, Hamming, Blackman, Blackman-Harris)

## L3: Mathematical Structures ✅ Complete
- Fourier transform of sampled signal (Poisson summation)
- Z-transform for ΣΔ loop filter analysis
- Quantization error as additive white noise (Bennett model)
- Statistical model of DNL/INL
- Probability density of quantization error (uniform [-q/2, q/2])
- Transfer characteristic polynomial model (Chebyshev)
- CDAC charge equation and Thévenin equivalent
- Comparator metastability exponential model
- Capacitor mismatch — Pelgrom's law
- Box-Muller transform for noise generation
- Cooley-Tukey radix-2 FFT

## L4: Fundamental Laws ✅ Complete
- **Ideal ADC SNR**: SNR[dB] = 6.02N + 1.76 (uniform quantization, full-scale sine)
- **Quantization noise power**: σ_q² = q² / 12 (derived via uniform distribution)
- **Nyquist-Shannon Theorem**: f_s > 2·BW for perfect reconstruction
- **Oversampling SNR gain**: ΔSNR = 10·log₁₀(OSR) [dB]
- **ΔΣ in-band noise**: IBN ≈ σ_q² · π^{2L} / ((2L+1) · OSR^{2L+1})
- **kT/C noise**: v_n_rms = √(kT/C)
- **Lee's criterion**: max|NTF(e^{jω})| < 2.0 for single-bit ΔΣ stability
- **Sinc interpolation (Whittaker-Shannon)**: x(t) = Σ x[n]·sinc((t-nT)/T)
- **Jitter SNR limit**: SNR_j = -20·log₁₀(2π·f_in·t_jitter_rms)
- **Bennett's conditions** for quantization noise whitening
- **Pelgrom's law**: σ(ΔC/C) = A_C / √(W·L)

## L5: Algorithms/Methods ✅ Complete
- SAR binary search algorithm (N bit trials)
- Sub-radix-2 redundant SAR search
- Lloyd-Max iterative optimal quantizer design
- CIC decimation filter (Hogenauer 1981)
- FIR decimation and interpolation (polyphase)
- Rational rate conversion L/M
- First-order noise shaping: NTF(z) = 1 - z⁻¹
- Second-order noise shaping: NTF(z) = (1 - z⁻¹)²
- ΔΣ MASH 1-1 cascade
- ΔΣ Butterworth NTF design
- Radix-2 FFT (Cooley-Tukey, decimation-in-time)
- Three-parameter sine fit for SINAD
- Interpolated FFT for frequency estimation
- Two-point ADC calibration (offset + gain)
- SAR capacitor mismatch calibration
- CIC compensation filter design
- DNL/INL computation from code transition voltages
- SFDR computation from spectrum (excl. harmonics)
- μ-law / A-law companding (ITU-T G.711)
- TPDF dithering for moment decorrelation
- Window function application and coherent gain normalization

## L6: Canonical Problems ✅ Complete
- **ADC resolution vs SNR sweep** (6 to 24 bits, jitter limit)
- **SAR ADC 10-bit simulation** (kT/C noise, comparator offset, mismatch)
- **ΣΔ audio ADC design** (1st/2nd/3rd order, OSR 32-256)
- **DAC glitch energy** at MSB transition (timing skew model)
- **R-2R ladder linearity** (resistor mismatch analysis)
- **Oversampling SNR improvement** (OSR 2-256, ΔENOB)
- **Band-pass sampling** (70-90 MHz IF band, valid fs ranges)
- **Sinc interpolation** reconstruction of sampled sinusoid
- **Full IEEE 1241 dynamic test suite** (SNR, SINAD, THD, SFDR, ENOB)
- **Two-tone IMD measurement** (IM2, IM3 products)

## L7: Applications ✅ Partial+ (3 applications)
- Sensor data acquisition (ramp input through SAR ADC with noise)
- Audio ADC (ΣΔ modulator for 20 kHz bandwidth, CD-quality target)
- Instrumentation grade ADC (16-bit ENOB target analysis)

## L8: Advanced Topics ✅ Partial+ (4 topics with implementation)
- Time-interleaved ADC architecture (documented in adc_dac_core.h)
- Stochastic ADC / dithering (TPDF non-subtractive implemented)
- Sub-radix-2 redundant SAR (implemented with error correction)
- ΣΔ MASH cascade architectures (MASH 1-1 implemented, ext to 2-1)

## L9: Research Frontiers ⚡ Partial (documented)
- AI-enhanced ADC calibration (noted in gap report)
- Compressive sensing ADC (noted in course tree)
- Quantum ADC concepts (documented in research frontiers section)
