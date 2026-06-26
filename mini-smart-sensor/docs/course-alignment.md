# Course Alignment — mini-smart-sensor

Nine-school curriculum mapping for Smart Sensor topics.

## MIT

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| 6.003 Signal Processing | Signal conditioning, filtering | L3/L5: Wheatstone bridge, digital filters (SMA, median, IIR) |
| 6.003 Signal Processing | Noise analysis, SNR | L3/L4: ADC quantization noise, SNR, MDS |
| 6.450 Digital Communications | Estimation theory | L5: Kalman filter, least-squares regression |
| 6.630 EM Waves | Sensor electromagnetic principles | L1: Hall effect, magnetoresistive sensors |

## Stanford

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| EE102A Signal Processing | Digital filter design | L5: IIR lowpass/highpass, moving average, median filter |
| EE102A Signal Processing | Parameter estimation | L5: Linear/polynomial regression, least squares |
| EE359 Wireless Communications | IoT sensor networks | L7: Multi-channel sensor node, wireless interfaces |
| EE247 Optical Devices | Photodetectors | L1: Photoelectric, optical transducers |

## UC Berkeley

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| EE16A/B Designing Information Devices | Sensor circuits, Wheatstone bridge | L3: Bridge analysis, resistive sensors |
| EE105 Analog Circuits | Instrumentation amplifiers | L2/L3: Analog front-end, amplifier gain calculation |
| EE123 Digital Signal Processing | Numerical methods | L5: Regression, Gaussian elimination, Horner's method |
| EE123 Digital Signal Processing | Kalman filtering | L5: 1D Kalman, IMU Kalman filter |

## Illinois

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| ECE 310 DSP | Digital filtering | L5: SMA, median, IIR filters |
| ECE 310 DSP | Statistical signal processing | L5: Allan variance, Welford algorithm |
| ECE 459 Communications | Estimation, detection | L5: Kalman filter, outlier rejection |

## Michigan

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| EECS 351 DSP | Sensor data processing | L5: Digital filters, feature extraction |
| EECS 351 DSP | Statistical analysis | L8: RMS, crest factor, skewness, ZCR |
| EECS 411 Microwave Circuits | Automotive radar/sensors | L7: TPMS, vibration monitoring |
| EECS 455 Communications | Stochastic estimation | L5: Kalman filter |

## Georgia Tech

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| ECE 4270 DSP | Adaptive filtering | L5: Kalman filter, complementary filter |
| ECE 6601 Communications | IoT sensor networks | L7: Multi-channel node, energy harvesting |
| ECE 6350 EM | Sensing principles | L1: Transducer types, Hall effect |

## TU Munich

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| Signal Processing | Filter theory | L5: IIR design, frequency response |
| Signal Processing | Estimation theory | L5: Kalman, least squares |
| Communications | Sensor networks | L7: IoT node, wireless interfaces |
| High-Frequency Engineering | RF sensors | L1: Wireless interfaces, radar sensors |

## ETH Zürich

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| 227-0427 Signal Processing | Statistical signal processing | L5/L8: Kalman, Allan variance, outlier detection |
| 227-0427 Signal Processing | Parameter estimation | L5: Regression, calibration |
| 227-0436 Communications | Stochastic estimation, fusion | L5: Weighted average, Kalman, complementary filter |
| 227-0455 EM/High-Frequency | RF energy harvesting | L8: Energy harvesting power budget |

## Tsinghua (清华)

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| 信号与系统 (Signals & Systems) | Filtering, SNR, sampling | L3/L4: ADC metrics, digital filters |
| 通信原理 (Communication Principles) | Estimation, detection | L5: Kalman filter, sensor fusion |
| 传感器技术 (Sensor Technology) | Transducer principles | L1: All transducer types, TEDS |
| 传感器技术 (Sensor Technology) | Signal conditioning | L2/L3: Bridge circuits, analog front-end |
| 数字信号处理 (DSP) | Digital filtering | L5: Moving average, median, IIR |
| 数字信号处理 (DSP) | Statistical signal processing | L5/L8: Welford, Allan variance, feature extraction |

## Cross-Cutting Themes

| Theme | Schools | Module Implementation |
|-------|---------|---------------------|
| Signal conditioning | MIT, Berkeley, Tsinghua | `sensor_conditioning.c` |
| Digital filtering | Stanford, Illinois, ETH, Tsinghua | `sensor_conditioning.c` (filters) |
| Estimation/Kalman | MIT, Stanford, Berkeley, ETH | `sensor_fusion.c` |
| Sensor fusion | All 9 schools | `sensor_fusion.c` |
| Calibration/Regression | Stanford, ETH, Tsinghua | `sensor_calibration.c` |
| IoT/Applications | Stanford, Michigan, Georgia Tech | `sensor_applications.c` |
| Energy harvesting | ETH, Georgia Tech | `sensor_advanced.c` |
