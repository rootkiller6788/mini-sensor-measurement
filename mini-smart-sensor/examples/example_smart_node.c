/**
 * @file    example_smart_node.c
 * @brief   L7: Complete Smart Sensor Node — Multi-Channel IoT Monitor
 *
 * Demonstrates a complete smart sensor node with:
 *   - 3 sensor channels (temperature, humidity, vibration)
 *   - Analog front-end configuration
 *   - Digital signal processing chain
 *   - Calibration application
 *   - Measurement buffer management
 *   - Alarm/threshold detection
 *   - Energy budget calculation
 *
 * This represents a typical IoT sensor node for industrial
 * condition monitoring (ISO 13374) or environmental monitoring
 * applications.
 *
 * Reference:
 *   IEEE 1451.0-2007 — Smart Transducer Interface Standard
 *   ISO 13374 — Condition Monitoring and Diagnostics of Machines
 *   Frank, R. (2013), "Understanding Smart Sensors", 3rd ed., Artech House
 */

#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <math.h>
#include <string.h>
#include "smart_sensor.h"
#include "sensor_conditioning.h"
#include "sensor_calibration.h"
#include "sensor_fusion.h"
#include "sensor_advanced.h"
#include <stdlib.h>

int main(void) {
    printf("============================================================\n");
    printf("  L7: Smart Sensor Node — Multi-Channel IoT Monitor\n");
    printf("  Industrial Condition Monitoring Application\n");
    printf("============================================================\n\n");

    /* ============================================================
     * Step 1: Create sensor node with 3 channels
     * ============================================================ */

    printf("--- Step 1: Initialize Sensor Node ---\n");

    ss_node_t node;
    int ret __attribute__((unused)) = ss_node_init(&node, 3, 0xA1B2C3D4E5F60001ULL);
    printf("Node initialized: ID = 0x%016llX, %d channels\n",
           (unsigned long long)node.node_id, node.n_channels);

    /* ============================================================
     * Step 2: Configure Channel 0 — Temperature (Thermocouple Type K)
     * ============================================================ */

    printf("\n--- Step 2: Configure Channel 0 — Thermocouple Type K ---\n");

    ss_sensor_spec_t spec_temp;
    ss_sensor_spec_init(&spec_temp,
                        SS_TRANSDUCER_THERMOELECTRIC,
                        SS_MEASURAND_TEMPERATURE,
                        0.041,      /* mV/degC (Type K Seebeck coefficient) */
                        1372.0,     /* Full scale: 0 to 1372 degC */
                        54.8);      /* 54.8 mV at 1372 degC */

    ss_accuracy_spec_t acc_temp;
    ss_accuracy_spec_init(&acc_temp, 0.75, 0.1);
    /* Type K: ±0.75% accuracy */

    ss_analog_frontend_t afe_temp;
    ss_analog_frontend_init(&afe_temp);
    afe_temp.amplifier_gain = 200.0;  /* Amplify mV to V range */
    afe_temp.amplifier_topology = 1;  /* Differential amplifier */
    afe_temp.adc_resolution_bits = 16; /* High resolution for precision */
    afe_temp.excitation_type = 3;      /* No excitation (self-powered TC) */

    ss_digital_processing_t dp_temp;
    ss_digital_processing_init(&dp_temp);
    dp_temp.lowpass_alpha = 0.05;     /* Heavy filtering for slow temp changes */

    ss_sensor_config_t cfg_temp;
    ss_sensor_config_init(&cfg_temp, &spec_temp, &acc_temp,
                          &afe_temp, &dp_temp, SS_INTERFACE_SPI);
    cfg_temp.sample_period_ms = 1000; /* 1 Hz sampling */

    /* Customize TEDS for this channel */
    strncpy(cfg_temp.teds.manufacturer_name, "Omega Engineering",
            sizeof(cfg_temp.teds.manufacturer_name) - 1);
    strncpy(cfg_temp.teds.model_number, "KMQXL-062G-6",
            sizeof(cfg_temp.teds.model_number) - 1);
    cfg_temp.teds.cal_interval_days = 365.0;

    int ch0 = ss_node_add_channel(&node, &cfg_temp);
    printf("Added channel %d: Type K Thermocouple, 0-1372 degC\n", ch0);

    /* ============================================================
     * Step 3: Configure Channel 1 — Humidity (Capacitive)
     * ============================================================ */

    printf("\n--- Step 3: Configure Channel 1 — Capacitive Humidity ---\n");

    ss_sensor_spec_t spec_hum;
    ss_sensor_spec_init(&spec_hum,
                        SS_TRANSDUCER_CAPACITIVE,
                        SS_MEASURAND_HUMIDITY,
                        32.0,       /* pF / %RH (Sensirion SHT4x typical) */
                        100.0,      /* 0-100% RH */
                        3200.0);    /* 3200 pF at 100% RH */

    ss_accuracy_spec_t acc_hum;
    ss_accuracy_spec_init(&acc_hum, 2.0, 0.05); /* ±2% RH accuracy */

    ss_analog_frontend_t afe_hum;
    ss_analog_frontend_init(&afe_hum);
    afe_hum.amplifier_topology = 3;  /* Transimpedance amplifier */
    afe_hum.adc_resolution_bits = 14;

    ss_sensor_config_t cfg_hum;
    ss_sensor_config_init(&cfg_hum, &spec_hum, &acc_hum,
                          &afe_hum, NULL, SS_INTERFACE_I2C);
    cfg_hum.interface_address = 0x44; /* SHT4x default I2C address */
    cfg_hum.sample_period_ms = 2000;  /* 0.5 Hz */

    int ch1 = ss_node_add_channel(&node, &cfg_hum);
    printf("Added channel %d: Capacitive Humidity, 0-100%% RH, I2C addr 0x44\n", ch1);

    /* ============================================================
     * Step 4: Configure Channel 2 — Vibration (Piezoelectric Accelerometer)
     * ============================================================ */

    printf("\n--- Step 4: Configure Channel 2 — Piezoelectric Accelerometer ---\n");

    ss_sensor_spec_t spec_vib;
    ss_sensor_spec_init(&spec_vib,
                        SS_TRANSDUCER_PIEZOELECTRIC,
                        SS_MEASURAND_ACCELERATION,
                        100.0,       /* mV/g (IEPE standard) */
                        50.0,        /* ±50 g */
                        5000.0);     /* 5000 mV at 50 g */

    ss_accuracy_spec_t acc_vib;
    ss_accuracy_spec_init(&acc_vib, 1.0, 0.001); /* ±1%, good repeatability */

    ss_analog_frontend_t afe_vib;
    ss_analog_frontend_init(&afe_vib);
    afe_vib.excitation_type = 2;     /* AC coupled (IEPE needs 4mA current) */
    afe_vib.excitation_current = 4.0;/* IEPE standard bias current */
    afe_vib.amplifier_gain = 10.0;
    afe_vib.filter_cutoff_freq = 10000.0; /* High bandwidth for vibration */
    afe_vib.adc_resolution_bits = 24;     /* High dynamic range */
    afe_vib.adc_sample_rate_hz = 25600.0; /* 25.6 kSPS for 10 kHz BW */

    ss_digital_processing_t dp_vib;
    ss_digital_processing_init(&dp_vib);
    dp_vib.moving_avg_window = 8;
    dp_vib.highpass_alpha = 0.001;     /* Remove DC, keep AC vibration */
    dp_vib.upper_threshold = 20.0;     /* Alarm at 20g */
    dp_vib.rate_of_change_limit = 100.0; /* g/s, for impact detection */

    ss_sensor_config_t cfg_vib;
    ss_sensor_config_init(&cfg_vib, &spec_vib, &acc_vib,
                          &afe_vib, &dp_vib, SS_INTERFACE_ANALOG_VOLTAGE);
    cfg_vib.sample_period_ms = 1; /* 1000 Hz effective (with decimation) */

    int ch2 = ss_node_add_channel(&node, &cfg_vib);
    printf("Added channel %d: IEPE Accelerometer, ±50g, 25.6 kSPS\n", ch2);

    /* ============================================================
     * Step 5: Simulate measurement acquisition and processing
     * ============================================================ */

    printf("\n--- Step 5: Simulated Measurements ---\n\n");

    /* Simulate 10 measurement cycles */
    {
        int cycle;
        for (cycle = 0; cycle < 10; cycle++) {
            /* --- Channel 0: Temperature --- */
            {
                /* Simulate: raw ADC reading → voltage → EMF → temperature */
                double true_temp = 150.0 + cycle * 10.0; /* 150-240 degC */
                double emf_mv = true_temp * 0.041;       /* 0.041 mV/degC */
                double raw_adc = emf_mv * afe_temp.amplifier_gain
                                 * 65536.0 / afe_temp.adc_reference_voltage
                                 + 50.0 * ((double)rand() / RAND_MAX - 0.5);

                /* Signal conditioning: convert ADC to voltage */
                double voltage = raw_adc * afe_temp.adc_reference_voltage
                                 / 65536.0 / afe_temp.amplifier_gain;

                /* Apply thermocouple calibration */
                ss_thermocouple_model_t tc_model;
                ss_cal_thermocouple_init(&tc_model, 0, 25.0);
                double temp_degc = ss_cal_thermocouple_to_temp(&tc_model,
                                                               voltage * 1000.0);

                /* Update sensor state */
                ss_sensor_state_update(&node.states[ch0], raw_adc,
                                       voltage, temp_degc,
                                       (uint64_t)(cycle * 1000000));

                if (cycle == 0 || cycle == 9) {
                    printf("  [Ch0-TEMP] Cycle %2d: %.1f degC (raw ADC: %.0f)\n",
                           cycle, temp_degc, raw_adc);
                }
            }

            /* --- Channel 1: Humidity --- */
            {
                double true_rh = 55.0 + 3.0 * sin(cycle * 0.5);
                double raw_adc = true_rh * 32.0 * 10.0  /* pF * gain */
                                 * 16384.0 / 3.3
                                 + 20.0 * ((double)rand() / RAND_MAX - 0.5);

                /* Simple linear calibration */
                ss_linear_model_t hum_model;
                ss_cal_linear_regression(
                    (double[]){0.0, 100.0}, (double[]){0.0, 32000.0}, 2, &hum_model);
                double rh = ss_cal_linear_apply(&hum_model, raw_adc * 3.3 / 16384.0);

                ss_sensor_state_update(&node.states[ch1], raw_adc,
                                       raw_adc, rh,
                                       (uint64_t)(cycle * 1000000));

                if (cycle == 0 || cycle == 9) {
                    printf("  [Ch1-HUM]  Cycle %2d: %.1f %%RH (raw ADC: %.0f)\n",
                           cycle, rh, raw_adc);
                }
            }

            /* --- Channel 2: Vibration --- */
            {
                double true_accel = 0.5 * sin(2.0 * M_PI * 50.0 * cycle * 0.001);
                double raw_adc = true_accel * 100.0 * afe_vib.amplifier_gain
                                 * 16777216.0 / afe_vib.adc_reference_voltage
                                 + 100.0 * ((double)rand() / RAND_MAX - 0.5);

                double accel_g = raw_adc * afe_vib.adc_reference_voltage
                                 / 16777216.0 / afe_vib.amplifier_gain / 100.0;

                /* Apply highpass filter to remove DC bias */
                static double prev_out = 0.0, prev_in = 0.0;
                double filtered = ss_cond_iir_highpass(&prev_out, &prev_in,
                                                       accel_g, 0.001);

                /* Check threshold */
                static int alarm_state = 0;
                int alarm = ss_cond_threshold_hysteresis(
                    &alarm_state, fabs(filtered), 19.0, 20.0);

                ss_sensor_state_update(&node.states[ch2], raw_adc,
                                       filtered, filtered,
                                       (uint64_t)(cycle * 1000));

                if (cycle == 0 || cycle == 9) {
                    printf("  [Ch2-VIB]  Cycle %2d: %.4f g (filtered), Alarm: %s\n",
                           cycle, filtered, alarm ? "YES" : "no");
                }
            }
        }
    }

    /* ============================================================
     * Step 6: Statistical summary from measurement buffer
     * ============================================================ */

    printf("\n--- Step 6: Statistical Summary ---\n\n");

    {
        int ch;
        const char *names[] = {"Temperature", "Humidity", "Vibration"};
        const char *units[] = {"degC", "%RH", "g"};

        for (ch = 0; ch < 3; ch++) {
            printf("Channel %d (%s):\n", ch, names[ch]);
            printf("  Samples:       %llu\n",
                   (unsigned long long)node.states[ch].sample_count);
            printf("  Current value: %.4f %s\n",
                   node.states[ch].current_engineering_units, units[ch]);
            printf("  Moving mean:   %.4f %s\n",
                   node.states[ch].moving_average, units[ch]);
            printf("\n");
        }
    }

    /* ============================================================
     * Step 7: Sensor Fusion — Combine redundant measurements
     * ============================================================ */

    printf("--- Step 7: Temperature Multi-Sensor Fusion ---\n\n");

    /* Simulate 3 temperature sensors in close proximity */
    double t_readings[] = {152.3, 151.1, 153.8};
    double t_variances[] = {0.25, 0.64, 0.09};

    printf("Redundant temperature sensors:\n");
    {
        int i;
        for (i = 0; i < 3; i++) {
            printf("  Sensor %d: %.1f degC (sigma=%.2f)\n",
                   i + 1, t_readings[i], sqrt(t_variances[i]));
        }
    }

    double fused_t, fused_var;
    fused_t = ss_fusion_weighted_average(t_readings, t_variances, 3, &fused_var);

    printf("\nFused temperature: %.2f degC (sigma=%.4f degC)\n",
           fused_t, sqrt(fused_var));

    /* ============================================================
     * Step 8: Energy Budget Calculation
     * ============================================================ */

    printf("\n--- Step 8: Energy Budget for Battery-Free Operation ---\n\n");

    /* Estimate power consumption */
    double p_active_mw = 50.0;   /* 50 mW active (MCU + radio + sensors) */
    double p_sleep_uw = 10.0;    /* 10 uW deep sleep */
    double duty = 0.01;          /* 1% duty cycle (wake 10ms every 1s) */
    double wakeup_uj = 100.0;    /* 100 uJ per wakeup */
    double wakeup_hz = 1.0;      /* Wake once per second */

    double p_avg_uw = ss_advanced_avg_power_uw(p_active_mw, p_sleep_uw,
                                                duty, wakeup_uj, wakeup_hz);
    printf("Power budget:\n");
    printf("  Active power:   %.0f mW\n", p_active_mw);
    printf("  Sleep power:    %.0f uW\n", p_sleep_uw);
    printf("  Duty cycle:     %.1f%%\n", duty * 100.0);
    printf("  Avg power:      %.1f uW\n\n", p_avg_uw);

    /* Check if indoor PV can sustain this */
    double pv_efficiency = 0.08;
    double pv_area = 15.0;        /* 15 cm^2 */
    double indoor_light = 200.0;  /* 200 uW/cm^2 indoor */
    double pv_power = ss_advanced_pv_power_uw(pv_efficiency, pv_area, indoor_light);

    printf("Energy harvesting (indoor PV):\n");
    printf("  Cell efficiency: %.0f%%\n", pv_efficiency * 100.0);
    printf("  Cell area:       %.0f cm^2\n", pv_area);
    printf("  Indoor light:    %.0f uW/cm^2\n", indoor_light);
    printf("  PV output:       %.1f uW\n", pv_power);
    printf("  Energy balance:  %+.1f uW (%s)\n",
           pv_power - p_avg_uw,
           pv_power > p_avg_uw ? "SUSTAINABLE" : "DEFICIT");

    printf("\n============================================================\n");
    printf("  Smart Sensor Node Example Complete\n");
    printf("============================================================\n");
    return 0;
}
