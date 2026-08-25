#!/usr/bin/env python3
"""
plot_energy_check_new_ligands.py

Plots and numerically summarizes the EM/NVT/NPT stability check for the
4 new-ligand test sequences (CDCA, GLCA, LCA3S, and the LCA control), from
the .xvg files produced by run_energy_check.sh (standard EM/NVT/NPT naming,
no _qfix suffix -- these systems were parameterized correctly from the
start, so there's no prior buggy version to distinguish from):

    EM/em_potential.xvg   -- potential energy during minimization
    NVT/nvt_temp.xvg      -- temperature during NVT
    NPT/npt_density.xvg   -- density during NPT
    NPT/npt_volume.xvg    -- volume during NPT

For each sequence, saves a 2x2 panel PNG to
    {WORKDIR}/{seq_id}_EM_EQ_energy.png
and prints a stability verdict per quantity (temperature target 300 K per
nvt.mdp; density/volume checked for plateau via second-half noise and
drift rather than an absolute target).

Usage:
    conda activate biosensors
    bash run_energy_check.sh seq_ids_new_ligands.txt
    python plot_energy_check_new_ligands.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt

BASE_DIR = "/scratch/alpine/ivta1597/LCA_boltz_models"
TARGET_TEMP_K = 300.0

SEQUENCES = [
    ("binders", "cdca_001_binder"),
    ("binders", "glca_001_binder"),
    ("binders", "lca3s_001_binder"),
    ("binders", "lca_001_binder"),
]


def read_xvg(path):
    """Skip GROMACS .xvg header lines (# and @); return (x, y) arrays."""
    x, y = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            x.append(float(parts[0]))
            y.append(float(parts[1]))
    return np.array(x), np.array(y)


def second_half_stats(x, y):
    """Mean/std/total-drift-fraction over the second half of the series
    (post-plateau): drift_frac is the fitted trend's total change across
    the window, as a fraction of the mean -- catches a slow monotonic
    trend that a low std alone could miss."""
    n = len(y)
    half = y[n // 2:]
    xhalf = x[n // 2:]
    mean, std = half.mean(), half.std()
    span = xhalf[-1] - xhalf[0] if len(xhalf) > 1 else 0.0
    if span > 0:
        slope = np.polyfit(xhalf, half, 1)[0]
        drift_frac = abs(slope * span / mean) if mean else float("inf")
    else:
        drift_frac = 0.0
    return mean, std, drift_frac


def verdict_temp(mean, std):
    if abs(mean - TARGET_TEMP_K) > 5 or std > 5:
        return "UNSTABLE", f"mean={mean:.1f} K (target {TARGET_TEMP_K:.0f} K), std={std:.2f} K"
    return "stable", f"mean={mean:.1f} K, std={std:.2f} K"


def verdict_plateau(label, mean, std, drift_frac, unit):
    # relative std (noise) and total fitted drift (trend) over the second
    # half of the series -- either one alone can miss instability the
    # other catches.
    rel_std = std / abs(mean) if mean else float("inf")
    if rel_std > 0.02 or drift_frac > 0.01:
        return "UNSTABLE", (f"{label} mean={mean:.2f} {unit}, std={std:.2f} "
                             f"({rel_std*100:.1f}% of mean), drift={drift_frac*100:.1f}% of mean")
    return "stable", (f"{label} mean={mean:.2f} {unit}, std={std:.2f} "
                       f"({rel_std*100:.2f}% of mean), drift={drift_frac*100:.2f}% of mean")


def main():
    summary = []
    for type_dir, seq_id in SEQUENCES:
        workdir = os.path.join(BASE_DIR, type_dir, seq_id)
        em_f = os.path.join(workdir, "EM", "em_potential.xvg")
        nvt_f = os.path.join(workdir, "NVT", "nvt_temp.xvg")
        dens_f = os.path.join(workdir, "NPT", "npt_density.xvg")
        vol_f = os.path.join(workdir, "NPT", "npt_volume.xvg")

        if not all(os.path.isfile(f) for f in [em_f, nvt_f, dens_f, vol_f]):
            print(f"SKIP {seq_id}: missing one or more .xvg files "
                  f"(run `bash run_energy_check.sh seq_ids_new_ligands.txt` first)")
            continue

        em_x, em_y = read_xvg(em_f)
        nvt_x, nvt_y = read_xvg(nvt_f)
        dens_x, dens_y = read_xvg(dens_f)
        vol_x, vol_y = read_xvg(vol_f)

        fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)

        axes[0, 0].plot(em_x, em_y, color="#648FFF")
        axes[0, 0].set_title("EM: Potential energy")
        axes[0, 0].set_xlabel("Step")
        axes[0, 0].set_ylabel("Potential (kJ/mol)")
        axes[0, 0].grid(True, alpha=0.4)

        axes[0, 1].plot(nvt_x, nvt_y, color="#DC267F")
        axes[0, 1].axhline(TARGET_TEMP_K, color="black", linestyle="--", linewidth=0.8)
        axes[0, 1].set_title("NVT: Temperature")
        axes[0, 1].set_xlabel("Time (ps)")
        axes[0, 1].set_ylabel("Temperature (K)")
        axes[0, 1].grid(True, alpha=0.4)

        axes[1, 0].plot(dens_x, dens_y, color="#FE6100")
        axes[1, 0].set_title("NPT: Density")
        axes[1, 0].set_xlabel("Time (ps)")
        axes[1, 0].set_ylabel("Density (kg/m^3)")
        axes[1, 0].grid(True, alpha=0.4)

        axes[1, 1].plot(vol_x, vol_y, color="#FFB000")
        axes[1, 1].set_title("NPT: Volume")
        axes[1, 1].set_xlabel("Time (ps)")
        axes[1, 1].set_ylabel("Volume (nm^3)")
        axes[1, 1].grid(True, alpha=0.4)

        fig.suptitle(f"{seq_id} -- EM/NVT/NPT stability check")

        out_png = os.path.join(workdir, f"{seq_id}_EM_EQ_energy.png")
        fig.savefig(out_png, dpi=300)
        plt.close(fig)

        temp_mean, temp_std, _ = second_half_stats(nvt_x, nvt_y)
        dens_mean, dens_std, dens_drift = second_half_stats(dens_x, dens_y)
        vol_mean, vol_std, vol_drift = second_half_stats(vol_x, vol_y)

        t_status, t_msg = verdict_temp(temp_mean, temp_std)
        d_status, d_msg = verdict_plateau("density", dens_mean, dens_std, dens_drift, "kg/m^3")
        v_status, v_msg = verdict_plateau("volume", vol_mean, vol_std, vol_drift, "nm^3")

        overall = "UNSTABLE" if "UNSTABLE" in (t_status, d_status, v_status) else "stable"
        print(f"\n=== {seq_id} ===  -> {out_png}")
        print(f"  Temperature : {t_status:8s} {t_msg}")
        print(f"  Density     : {d_status:8s} {d_msg}")
        print(f"  Volume      : {v_status:8s} {v_msg}")
        print(f"  OVERALL     : {overall}")

        summary.append((seq_id, overall))

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for seq_id, overall in summary:
        print(f"  {overall:10s} {seq_id}")


if __name__ == "__main__":
    main()
