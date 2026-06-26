/**
 * @file test_dynamics.c
 * @brief Tests for L2-L3: Spring-mass-damper dynamics, frequency response,
 *        damping models, Coriolis force, Stoney's formula
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include "mems_dynamics.h"

#define PASS(msg) printf("  PASS: %s\n", msg)

static void test_smd_step(void) {
    printf("test_smd_step:\n");
    MemsSpringMassDamper sys = {
        .mass = 1e-9, .spring_k = 100.0, .damping_c = 1e-8,
        .displacement = 0.0, .velocity = 0.0, .acceleration = 0.0,
        .force_external = 1e-6, .time = 0.0
    };
    /* Step once: ẍ = (1e-6 - 0 - 0)/1e-9 = 1000 m/s² */
    mems_smd_step(&sys, 1e-6);
    assert(fabs(sys.acceleration - 1000.0) < 10.0);
    PASS("ẍ = F/m = 1000 m/s²");
    assert(sys.velocity > 0.0);
    PASS("Velocity increased");
    assert(sys.displacement > 0.0);
    PASS("Displacement increased");
    assert(fabs(sys.time - 1e-6) < 1e-12);
    PASS("Time advanced by dt");

    /* Null pointer safety */
    mems_smd_step(NULL, 0.001);
    PASS("NULL pointer handled");
}

static void test_harmonic_response(void) {
    printf("test_harmonic_response:\n");
    /* m=1, k=100, c=1, F0=10, ω=5 rad/s → natural ω_n = sqrt(100) = 10
     * X = 10 / sqrt((100-25)² + (1*5)²) = 10 / sqrt(75²+25) ≈ 10/75.17 ≈ 0.133 */
    double amp = mems_harmonic_amplitude(1.0, 100.0, 1.0, 10.0, 5.0);
    assert(fabs(amp - 0.133) < 0.01);
    PASS("Harmonic amplitude ≈ 0.133");

    double phase = mems_harmonic_phase(1.0, 100.0, 1.0, 5.0);
    /* φ = atan2(5, 75) ≈ 0.0666 rad */
    assert(fabs(phase - 0.0666) < 0.001);
    PASS("Harmonic phase ≈ 0.0666 rad");

    /* At resonance ω = ω_n = 10, amplitude should be large */
    double amp_res = mems_harmonic_amplitude(1.0, 100.0, 1.0, 10.0, 10.0);
    /* X = 10 / sqrt(0 + (10)²) = 10/10 = 1.0 */
    assert(fabs(amp_res - 1.0) < 0.01);
    PASS("Resonance amplitude = F0/(c·ω_n)");
}

static void test_frequency_response(void) {
    printf("test_frequency_response:\n");
    /* r = 0 (DC): |H| = 1/1 = 1 */
    double mag0 = mems_freq_response_mag(0.0, 0.1);
    assert(fabs(mag0 - 1.0) < 1e-10);
    PASS("|H(0)| = 1");

    /* r = 1 (resonance): |H| = 1/(2ζ) = 1/0.2 = 5 for ζ=0.1 */
    double mag1 = mems_freq_response_mag(1.0, 0.1);
    assert(fabs(mag1 - 5.0) < 0.01);
    PASS("|H(1)| = 1/(2ζ) = 5");

    /* r → large: |H| ≈ 1/r² → small */
    double mag_large = mems_freq_response_mag(10.0, 0.1);
    assert(mag_large < 0.02);
    PASS("|H(10)| << 1");

    /* Phase at resonance → -π/2 as ζ → 0+.
     * For ζ near 0, phase = -atan2(2ζ, 0+) = -π/2 */
    double phase = mems_freq_response_phase(1.0, 1e-10);
    assert(fabs(phase + M_PI/2.0) < 1e-8);
    PASS("∠H(1) → -90° for ζ→0");
}

static void test_natural_freq(void) {
    printf("test_natural_freq:\n");
    double wn = mems_natural_freq_rad(1.0, 100.0);
    assert(fabs(wn - 10.0) < 1e-10);
    PASS("ω_n = sqrt(100/1) = 10 rad/s");

    double wd = mems_damped_freq_rad(10.0, 0.1);
    assert(fabs(wd - 10.0 * sqrt(0.99)) < 1e-10);
    PASS("ω_d = ω_n·sqrt(1-ζ²)");

    /* Overdamped → 0 */
    assert(mems_damped_freq_rad(10.0, 1.5) == 0.0);
    PASS("Critically/overdamped → 0");
}

static void test_static_displacement(void) {
    printf("test_static_displacement:\n");
    /* x = m·a/k = 1e-9 * 9.81 / 100 = 9.81e-11 m */
    double x = mems_static_displacement(1e-9, 100.0, 9.81);
    assert(fabs(x - 9.81e-11) < 1e-13);
    PASS("x = m·a/k = 98.1 pm");
}

static void test_squeeze_film_damping(void) {
    printf("test_squeeze_film_damping:\n");
    /* Air: μ=1.85e-5 Pa·s, A=1e-8 m², gap=2e-6 m, w=100e-6 m, L=100e-6 m
     * β ≈ 1 - 0.21 + 0.09 = 0.88 (for w/L=1)
     * c ≈ 1.85e-5 * 1e-8 * 0.88 / (2e-6)³ ≈ 2.0e4 N·s/m */
    double c = mems_squeeze_film_damping(1.85e-5, 1e-8, 2e-6, 100e-6, 100e-6);
    assert(c > 1e4 && c < 3e4);
    PASS("c_sf ≈ 2×10⁴ N·s/m for μN gap");

    assert(mems_squeeze_film_damping(1.85e-5, 1e-8, 0.0, 100e-6, 100e-6) == 0.0);
    PASS("Zero gap → 0");
}

static void test_couette_damping(void) {
    printf("test_couette_damping:\n");
    double c = mems_couette_damping(1.85e-5, 1e-8, 2e-6);
    /* c = 1.85e-5 * 1e-8 / 2e-6 = 9.25e-8 */
    assert(fabs(c - 9.25e-8) < 1e-10);
    PASS("c_couette = μA/g");
}

static void test_thermoelastic_damping(void) {
    printf("test_thermoelastic_damping:\n");
    /* Si beam: E=170e9, α=2.6e-6, T=300, ρ=2330, Cp=700,
     * ω=2π·1000, b=10e-6, χ=9e-5
     * Q⁻¹ should be small but positive */
    double q_inv = mems_thermoelastic_q_inv(170e9, 2.6e-6, 300.0, 2330.0,
                                             700.0, 2.0 * M_PI * 1000.0,
                                             10e-6, 9e-5);
    assert(q_inv >= 0.0);
    PASS("Q⁻¹_TED ≥ 0");
    assert(q_inv < 1.0);
    PASS("Q⁻¹_TED < 1 (reasonable Q)");
}

static void test_coriolis_force(void) {
    printf("test_coriolis_force:\n");
    /* m=1e-9 kg, v=0.01 m/s, Ω=1 rad/s → F_c = 2*1e-9*0.01*1 = 2e-11 N */
    double Fc = mems_coriolis_force(1e-9, 0.01, 1.0);
    assert(fabs(Fc - 2e-11) < 1e-14);
    PASS("F_c = 2·m·v·Ω = 2e-11 N");

    /* Zero rate → zero force */
    assert(mems_coriolis_force(1e-9, 0.01, 0.0) == 0.0);
    PASS("Zero Ω → F_c = 0");
}

static void test_gyro_scale_factor(void) {
    printf("test_gyro_scale_factor:\n");
    /* m=1e-9, v=0.01, Q=1000, k=100 → SF = 2*1e-9*0.01*1000/100 = 2e-10 */
    double sf = mems_gyro_scale_factor(1e-9, 0.01, 1000.0, 100.0);
    assert(fabs(sf - 2e-10) < 1e-13);
    PASS("SF = 2·m·v·Q/k");

    assert(mems_gyro_scale_factor(1e-9, 0.01, 1000.0, 0.0) == 0.0);
    PASS("Zero k_sense → 0");
}

static void test_stoney_formula(void) {
    printf("test_stoney_formula:\n");
    /* Si substrate: E=170e9 Pa, t_s=500e-6 m, ν=0.22, t_f=1e-6 m, R=10 m
     * σ = 170e9*(500e-6)²/(6*(1-0.22)*1e-6*10) ≈ 170e9*2.5e-7/(6*0.78*1e-5)
     *   = 42500/(4.68e-5) ≈ 908 MPa (tensile) */
    double stress = mems_stoney_stress(170e9, 500e-6, 0.22, 1e-6, 10.0);
    assert(stress > 1e7 && stress < 1e10);
    PASS("σ_f in MPa range for thin film");

    /* Inverse: curvature from known stress */
    double R = mems_stoney_curvature(170e9, 500e-6, 0.22, 1e-6, stress);
    assert(fabs(R - 10.0) < 1e-6);
    PASS("Inverse Stoney recovers R");

    /* Invalid inputs */
    assert(mems_stoney_stress(170e9, 0.0, 0.22, 1e-6, 10.0) == 0.0);
    PASS("Zero t_s → 0");
}

static void test_capacitive_readout(void) {
    printf("test_capacitive_readout:\n");
    /* C1=210fF, C2=190fF, d0=2μm → x = 2μm * 20/400 = 0.1μm */
    double x = mems_capacitive_displacement(210e-15, 190e-15, 2e-6);
    assert(fabs(x - 0.1e-6) < 1e-9);
    PASS("x = d₀·(C1-C2)/(C1+C2) = 0.1 μm");

    /* Equal caps → zero displacement */
    double x0 = mems_capacitive_displacement(200e-15, 200e-15, 2e-6);
    assert(fabs(x0) < 1e-12);
    PASS("C1=C2 → x=0");
}

static void test_pullin_voltage(void) {
    printf("test_pullin_voltage:\n");
    /* k=100 N/m, d0=2μm, ε=8.854e-12, A=1e-8 m²
     * V_pi = sqrt(8*100*(2e-6)³/(27*8.854e-12*1e-8))
     *      = sqrt(8*100*8e-18/(27*8.854e-20))
     *      = sqrt(6.4e-15/(2.39e-18)) = sqrt(2678) ≈ 51.8 V */
    double Vpi = mems_pullin_voltage(100.0, 2e-6, 8.854e-12, 1e-8);
    assert(Vpi > 10.0 && Vpi < 200.0);
    PASS("V_pi in reasonable range");

    assert(isinf(mems_pullin_voltage(0.0, 2e-6, 8.854e-12, 1e-8)));
    PASS("Zero k → INF");
}

int main(void) {
    printf("=== test_dynamics ===\n\n");
    test_smd_step();
    test_harmonic_response();
    test_frequency_response();
    test_natural_freq();
    test_static_displacement();
    test_squeeze_film_damping();
    test_couette_damping();
    test_thermoelastic_damping();
    test_coriolis_force();
    test_gyro_scale_factor();
    test_stoney_formula();
    test_capacitive_readout();
    test_pullin_voltage();
    printf("\n=== All dynamics tests PASSED ===\n");
    return 0;
}
