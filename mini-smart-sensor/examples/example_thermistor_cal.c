/**
 * @file    example_thermistor_cal.c
 * @brief   L6: Thermistor calibration and temperature measurement
 *
 * Demonstrates the complete thermistor signal chain:
 *   1. Steinhart-Hart 3-point calibration
 *   2. Beta parameter simplified calibration
 *   3. Resistance-to-temperature conversion
 *   4. Temperature compensation of readings
 *   5. Multi-point calibration table lookup
 *
 * This is a canonical problem in sensor calibration:
 *   Given calibration data (R, T) pairs, determine the sensor model
 *   parameters and use them to convert resistance readings to temperature.
 *
 * Reference:
 *   Steinhart & Hart (1968), Deep Sea Research 15(4):497-503
 *   NTC Thermistor Datasheet: Murata NCP series
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "smart_sensor.h"
#include "sensor_calibration.h"
#include "sensor_conditioning.h"

int main(void) {
    printf("============================================================\n");
    printf("  L6: Thermistor Calibration and Temperature Measurement\n");
    printf("  Steinhart-Hart Model + Beta Parameter Fitting\n");
    printf("============================================================\n\n");

    /* ============================================================
     * Part 1: Steinhart-Hart 3-point calibration
     *
     * Calibration data for a 10k NTC thermistor:
     *   At   0 degC (273.15 K): R = 27,280 Ohm
     *   At  25 degC (298.15 K): R = 10,000 Ohm
     *   At  50 degC (323.15 K): R =  3,588 Ohm
     * ============================================================ */

    printf("--- Part 1: Steinhart-Hart 3-Point Calibration ---\n\n");

    ss_steinhart_hart_t sh_model;
    int ret = ss_cal_steinhart_hart_3pt(
        27280.0, 273.15,   /* Point 1: 0 degC */
        10000.0, 298.15,   /* Point 2: 25 degC */
         3588.0, 323.15,   /* Point 3: 50 degC */
        &sh_model);

    if (ret != 0) {
        printf("ERROR: Calibration failed (degenerate points?)\n");
        return 1;
    }

    printf("Calibrated Steinhart-Hart coefficients:\n");
    printf("  A = %e\n", sh_model.coeff_A);
    printf("  B = %e\n", sh_model.coeff_B);
    printf("  C = %e\n\n", sh_model.coeff_C);

    /* Verify at calibration points */
    printf("Verification at calibration points:\n");
    printf("  T(27280 Ohm) = %.2f K  (expected 273.15 K)\n",
           ss_cal_steinhart_hart_temp(&sh_model, 27280.0));
    printf("  T(10000 Ohm) = %.2f K  (expected 298.15 K)\n",
           ss_cal_steinhart_hart_temp(&sh_model, 10000.0));
    printf("  T( 3588 Ohm) = %.2f K  (expected 323.15 K)\n\n",
           ss_cal_steinhart_hart_temp(&sh_model, 3588.0));

    /* Test at intermediate points */
    printf("Temperature vs. Resistance (Steinhart-Hart):\n");
    printf("  %-12s %-12s %-12s\n", "R (Ohm)", "T (K)", "T (degC)");
    printf("  %-12s %-12s %-12s\n", "--------", "--------", "--------");

    {
        double test_R[] = {50000.0, 20000.0, 15000.0, 5000.0, 2000.0};
        int i;
        for (i = 0; i < 5; i++) {
            double tK = ss_cal_steinhart_hart_temp(&sh_model, test_R[i]);
            printf("  %-12.0f %-12.2f %-12.2f\n", test_R[i], tK, tK - 273.15);
        }
    }

    /* ============================================================
     * Part 2: Beta parameter simplified model
     *
     * Using only R0 at T0=25 degC and Beta value.
     * Beta = ln(R1/R0) / (1/T1 - 1/T0)
     * ============================================================ */

    printf("\n--- Part 2: Beta Parameter Simplified Model ---\n\n");

    ss_steinhart_hart_t beta_model;
    ret = ss_cal_thermistor_beta(10000.0, 298.15, 3588.0, 323.15, &beta_model);
    if (ret == 0) {
        printf("Beta model parameters:\n");
        printf("  R0 (at 25 degC) = %.0f Ohm\n", beta_model.r_at_25c);
        printf("  Beta            = %.1f K\n\n", beta_model.beta_value);

        printf("Beta model vs Steinhart-Hart comparison:\n");
        printf("  %-12s %-12s %-12s %-12s\n",
               "R (Ohm)", "T_SH (degC)", "T_Beta (degC)", "Delta (degC)");
        printf("  %-12s %-12s %-12s %-12s\n",
               "--------", "----------", "-----------", "-----------");

        {
            double test_R[] = {30000.0, 20000.0, 10000.0, 5000.0, 3000.0};
            int i;
            for (i = 0; i < 5; i++) {
                double t_SH = ss_cal_steinhart_hart_temp(&sh_model, test_R[i]) - 273.15;
                double t_B  = ss_cal_steinhart_hart_temp(&beta_model, test_R[i]) - 273.15;
                printf("  %-12.0f %-12.2f %-12.2f %-12.4f\n",
                       test_R[i], t_SH, t_B, fabs(t_SH - t_B));
            }
        }
    }

    /* ============================================================
     * Part 3: Polynomial calibration for nonlinear sensors
     *
     * For a pressure sensor with known nonlinearity, use polynomial
     * regression to create a calibration model.
     * ============================================================ */

    printf("\n--- Part 3: Polynomial Calibration for Pressure Sensor ---\n\n");

    /* Synthetic calibration data: pressure sensor with slight nonlinearity */
    /* True pressure (kPa):    0,   20,  40,   60,   80,  100 */
    /* Sensor output (mV):    0,  198, 404,  618,  840, 1070 */

    double p_true[]  = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0};
    double v_sensor[] = {0.0, 198.0, 404.0, 618.0, 840.0, 1070.0};
    ss_polynomial_model_t poly;

    ret = ss_cal_polynomial_regression(p_true, v_sensor, 6, 2, &poly);
    if (ret == 0) {
        printf("2nd-order polynomial fit:\n");
        printf("  p(mV) = %.6f * mV^2 + %.6f * mV + %.6f\n",
               poly.coeffs[0], poly.coeffs[1], poly.coeffs[2]);
        printf("  R-squared = %.6f\n", poly.r_squared);
        printf("  Max residual = %.4f mV\n\n", poly.max_residual);

        /* Apply calibration */
        printf("Calibration verification:\n");
        printf("  %-12s %-12s %-12s %-12s\n",
               "True (kPa)", "Sensor (mV)", "Calc (kPa)", "Error (kPa)");
        printf("  %-12s %-12s %-12s %-12s\n",
               "----------", "----------", "----------", "----------");

        {
            int i;
            for (i = 0; i < 6; i++) {
                double p_calc = ss_cal_polynomial_apply(&poly, v_sensor[i]);
                printf("  %-12.1f %-12.1f %-12.2f %-12.4f\n",
                       p_true[i], v_sensor[i], p_calc, fabs(p_calc - p_true[i]));
            }
        }

        /* Inverse: find mV for target pressure 50 kPa */
        double mv_target = ss_cal_polynomial_inverse(&poly, 50.0, 500.0, 1e-6, 50);
        printf("\nInverse: For target pressure 50 kPa, sensor should read %.1f mV\n",
               mv_target);
    }

    /* ============================================================
     * Part 4: Temperature compensation in action
     * ============================================================ */

    printf("\n--- Part 4: Temperature Compensation ---\n\n");

    /* Simulate a pressure sensor at different temperatures */
    /* Zero drift: 0.1 kPa/degC, Span drift: 200 ppm/degC (0.0002/degC) */

    printf("Uncompensated vs Compensated readings:\n");
    printf("  %-10s %-12s %-12s %-12s\n",
           "T (degC)", "Raw (kPa)", "Corr (kPa)", "Error (kPa)");
    printf("  %-10s %-12s %-12s %-12s\n",
           "--------", "----------", "----------", "----------");

    {
        double temps[] = {-20.0, 0.0, 25.0, 50.0, 85.0};
        double true_pressure = 100.0; /* Constant true value */
        int i;
        for (i = 0; i < 5; i++) {
            /* Simulate raw reading affected by temperature */
            double dt = temps[i] - 25.0;
            double zero_shift = 0.1 * dt;
            double span_factor = 1.0 + 0.0002 * dt;
            double raw = true_pressure * span_factor + zero_shift;

            double compensated = ss_cond_temp_compensate(
                raw, temps[i], 25.0, 0.1, 0.0002);

            printf("  %-10.1f %-12.3f %-12.3f %-12.4f\n",
                   temps[i], raw, compensated, fabs(compensated - true_pressure));
        }
    }

    printf("\n============================================================\n");
    printf("  Thermistor Calibration Example Complete\n");
    printf("============================================================\n");
    return 0;
}
