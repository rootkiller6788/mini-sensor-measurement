# Course Alignment — mini-sensor-principle

## Nine-School Curriculum Mapping

| School | Course | Covered Topics |
|--------|--------|---------------|
| **MIT** | 6.003 Signal Processing | Dynamic sensor models (1st/2nd order), transfer functions, frequency response, Bode plots |
| | 6.630 EM Waves | Sensor noise physics (Johnson, shot), electromagnetic sensing principles |
| **Stanford** | EE102A Signal Processing | Digital filtering for sensors, moving average, exponential filter, median filter |
| | EE359 Wireless | Noise figure cascade (Friis), Allan variance for oscillator/sensor characterization |
| **Berkeley** | EE16A/B Circuits | Wheatstone bridge, instrumentation amplifier, signal conditioning |
| | EE105 Analog | Op-amp gain configurations, anti-alias filter design, current loop |
| | EE123 DSP | Allan variance, digital filtering, sensor linearization |
| **Illinois** | ECE 310 DSP | Sensor noise PSD, 1/f noise, noise equivalent bandwidth |
| | ECE 459 Comm | Noise figure, noise temperature, cascade analysis |
| **Michigan** | EECS 351 DSP | Moving average, exponential filter, frequency-domain analysis |
| | EECS 411 Microwave | Noise figure, thermal noise, detector sensitivity |
| **Georgia Tech** | ECE 4270 DSP | Sensor calibration (LS, polynomial), temperature compensation |
| | ECE 6601 Comm | Noise models (Johnson, shot, flicker), NEP, D* |
| **TU Munich** | Signal Processing | Filter design for sensors, sampling theorem application |
| | High-Frequency Engineering | Thermal noise, shot noise, noise figure |
| **ETH** | 227-0427 Signal Processing | Allan variance, noise characterization, sensor dynamics |
| | 227-0436 Comm | Noise in measurement systems, detectivity limits |
| **清华** | 信号与系统 | 1st/2nd-order dynamic response, frequency response |
| | 通信原理 | Noise figure, noise temperature, SNR analysis |
| | 数字信号处理 | Digital filtering, moving average, exponential smoothing |

## Reference Textbooks Mapped

| Textbook | Chapters | Implementation |
|----------|----------|---------------|
| Fraden, "Handbook of Modern Sensors" | Ch. 1-2 (Definitions), Ch. 3 (Transfer Functions) | sensor_defs, sensor_transfer |
| Doebelin, "Measurement Systems" | Ch. 3-4 (Dynamic Models), Ch. 5 (Bridge) | sensor_transfer, sensor_bridge |
| Pallas-Areny & Webster, "Sensors and Signal Conditioning" | Ch. 2-4 (Bridge, Amplifiers) | sensor_bridge, sensor_signal_cond |
| Bevington & Robinson, "Data Reduction and Error Analysis" | Ch. 6-8 (Least-Squares) | sensor_calibration |
| NIST Monograph 175 | Thermocouple Reference Functions | sensor_transfer (tc_model) |
| IEC 60751 | RTD Standard | sensor_transfer (rtd_cvd) |
| IEEE Std 952-1997 | Allan Variance for Gyros | sensor_noise (allan_variance) |
| ISO/IEC Guide 98-3 (GUM) | Uncertainty Budget | sensor_defs (uncertainty_budget) |
