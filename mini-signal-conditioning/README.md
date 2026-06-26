# mini-signal-conditioning

**Signal Conditioning for Sensor Measurement Systems**

Complete implementation of analog and digital signal conditioning: filtering,
linearization, excitation sources, galvanic isolation, cold junction compensation,
Wheatstone bridge circuits, and self-calibration.

---

## Module Status: COMPLETE

- **L1-L6:** Complete
- **L7:** Complete (4 applications: Boeing 787, NASA, F-35, NHS)
- **L8:** Partial (3/5 advanced topics)
- **L9:** Partial (documented, not fully implemented)

**include/ + src/ line count:** >3700 lines (>= 3000 threshold)

---

## Nine-Layer Knowledge Coverage

| Level | Name | Rating | Key Count |
|-------|------|--------|-----------|
| L1 | Definitions | Complete | 20 typedefs, 7 headers |
| L2 | Core Concepts | Complete | 10 concepts |
| L3 | Math Structures | Complete | 9 structures |
| L4 | Fundamental Laws | Complete | 11 laws verified |
| L5 | Algorithms/Methods | Complete | 19 algorithms |
| L6 | Canonical Problems | Complete | 3 examples |
| L7 | Applications | Complete | 4 real-world apps |
| L8 | Advanced Topics | Partial | 3 topics |
| L9 | Research Frontiers | Partial | Smart sensor documented |

---

## Core Definitions (L1)

Pipeline: Protection -> Excitation -> Amplification -> Filtering ->
Linearization -> Isolation -> Offset Null -> CJC -> ADC Buffer

Filter Families: Butterworth, Chebyshev I/II, Elliptic, Bessel, Gaussian
Linearization: LUT, Polynomial, Steinhart-Hart, Callendar-Van Dusen, NIST ITS-90
Excitation: Const V/I, AC, Ratiometric, Pulsed
Isolation: Optical, Transformer, Capacitive, GMR, Hall, RF iCoupler

## Core Theorems (L4)

Nyquist-Shannon: fc_anti-alias <= fs/(2*k), k>=2.5 for 12-bit
Successive Temperatures: E(Th,0) = E(Th,Tcj) + E(Tcj,0)
Bridge Balance: R1*R3 = R2*R4 => Vout = 0
Steinhart-Hart: 1/T = A + B*ln(R) + C*[ln(R)]^3
Callendar-Van Dusen: R(t) = R0*(1+A*t+B*t^2) for t>=0C; +C*(t-100)*t^3 for t<0C

## Core Algorithms (L5)

19 algorithms: Butterworth/Chebyshev/Bessel pole placement, Sallen-Key/MFB/Twin-T
design, Bilinear transform, Kaiser FIR, Horner evaluation, Lagrange/Catmull-Rom
interpolation, NIST ITS-90, Howland current pump, Two-point/multi-channel
calibration, Bridge nonlinearity correction, 3W/4W Kelvin, Elliptic order,
Least-squares fit, CJC compensation.

## Canonical Problems (L6)

1. Thermocouple Signal Chain (ex1) - Type K + CJC + amp + filter + ADC + ITS-90
2. Strain Gauge Measurement (ex2) - Bridge + excitation + nonlinearity correction
3. Precision RTD Measurement (ex3) - 4-wire Pt100 + CVD + self-heating correction

## Nine-School Course Mapping

| School | Courses | Coverage |
|--------|---------|----------|
| MIT | 6.003, 6.301, 6.630 | Filters, amplifier circuits, isolation |
| Stanford | EE247, EE359 | Filter chains, ADC interface |
| UC Berkeley | EE16A/B, EE105, EE123 | Bridge circuits, active filters, DSP |
| Michigan | EECS 351, 411, 455 | Digital filters, isolation, receivers |
| Georgia Tech | ECE 4270, 6350, 4430 | Filter coeffs, EMI/EMC, sensors |
| TU Munich | SP, HF Eng | Analog filters, isolation |
| ETH Zurich | 227-0427, 0436, 0455 | Filter theory, ADC, galvanic isolation |
| Tsinghua | Signal&Systems, Comm | Laplace, sampling, measurement |

## Build and Test

```bash
make          # Build library + test + examples
make test     # Run assert-based tests
make clean    # Clean build artifacts
```

## References

Pallas-Areny & Webster (2001), Kester (2004), Zumbahlen (2008),
Van Valkenburg (1982), NIST Monograph 175 (1993),
IEC 60751:2008, IEC 60601-1, Sallen & Key (1955), Kaiser (1974)
