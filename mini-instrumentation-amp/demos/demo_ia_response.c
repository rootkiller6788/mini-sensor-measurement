/**
 * demo_ia_response.c - IA Response Visualization Demo
 *
 * Outputs frequency response and gain sweep data in CSV format
 * for plotting and visualization.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include "inst_amp_analysis.h"

int main(void) {
    printf("=== IA Frequency Response Demo ===\n");
    printf("# CMRR vs Frequency for different topologies\n");
    printf("# freq_hz,cmrr_3opamp_db,cmrr_2opamp_db,cmrr_currmode_db\n");

    for (int f = 10; f <= 100000; f *= 2) {
        double c3 = ia_topology_cmrr_at_freq(IA_TOPO_THREE_OPAMP, 90.0, f, 1.0, 100.0);
        double c2 = ia_topology_cmrr_at_freq(IA_TOPO_TWO_OPAMP, 80.0, f, 1.0, 100.0);
        double cm = ia_topology_cmrr_at_freq(IA_TOPO_CURRENT_MODE, 70.0, f, 1.0, 100.0);
        printf("%d,%.1f,%.1f,%.1f\n", f, c3, c2, cm);
    }

    printf("\n# Noise vs Source Resistance\n");
    printf("# rs_ohm,en_total_nv_rhz,en_thermal_nv_rhz\n");
    for (int rs = 100; rs <= 100000; rs *= 2) {
        double en_total = noise_total_rti_nv_per_rhz(9.0, 0.1, rs, 500.0, 25000.0, 100.0, 27.0);
        double en_therm = noise_thermal_nv_per_rhz(rs, 27.0);
        printf("%d,%.2f,%.2f\n", rs, en_total, en_therm);
    }

    printf("\n# Gain Sweep (3-op-amp with R_G variation)\n");
    printf("# rg_ohm,gain,bw_khz\n");
    ia_three_opamp_config_t cfg = {25000.0, 10000.0, 10000.0, 1000.0, 0.1};
    for (double rg = 50000.0; rg >= 50.0; rg /= 2.0) {
        cfg.rg_ohm = rg;
        double G = ia_three_opamp_gain(&cfg);
        double bw = ia_topology_bandwidth(IA_TOPO_THREE_OPAMP, G, 1.0, 10.0) / 1000.0;
        printf("%.1f,%.1f,%.1f\n", rg, G, bw);
    }

    printf("\n=== Demo Complete ===\n");
    return 0;
}
