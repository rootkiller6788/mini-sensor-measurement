#!/usr/bin/env python3
"""
Demo: ADC SNR Curve visualization.
Generates the classic 6.02N+1.76 SNR vs resolution plot.
"""
import math
import sys

def main():
    print("ADC Ideal SNR Curve")
    print("=" * 50)
    print(f"{'Bits':>6s}  {'SNR[dB]':>10s}  {'ENOB':>8s}")
    print("-" * 30)
    for n in range(6, 25):
        snr = 6.02 * n + 1.76
        enob = (snr - 1.76) / 6.02
        print(f"{n:6d}  {snr:10.2f}  {enob:8.2f}")
    print()
    print("Formula: SNR[dB] = 6.02·N + 1.76")
    return 0

if __name__ == "__main__":
    sys.exit(main())
