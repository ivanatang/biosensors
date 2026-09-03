#!/usr/bin/env python3
"""
plot_energy_check.py

Plots and numerically summarizes the EM/NVT/NPT stability check from the
.xvg files produced by run_energy_check.sh:

    EM{suffix}/em_potential{suffix}.xvg   -- potential energy during minimization
    NVT{suffix}/nvt_temp{suffix}.xvg      -- temperature during NVT
    NPT{suffix}/npt_density{suffix}.xvg   -- density during NPT
    NPT{suffix}/npt_volume{suffix}.xvg    -- volume during NPT

For each sequence, saves a 2x2 panel PNG to
    {WORKDIR}/{seq_id}_EM_EQ_energy{suffix}.png
and prints a stability verdict per quantity (temperature target 300 K per
nvt.mdp; density/volume checked for plateau via second-half noise and drift
rather than an absolute target).

Reads sequences from a seq_ids.txt-style file (same format the bash scripts
use), same convention as run_energy_check.sh -- pass a smaller file, or a
process-substitution filter, to scope to specific sequences rather than
creating a separate tracked seq-list file.

Usage:
    conda activate biosensors
    python plot_energy_check.py --seq-list seq_ids.txt --suffix _qfix \
        --filter bind_022_binder bind_019_binder bind_020_binder \
                 nonb_006_nb nonb_008_nb nonb_009_nb
    python plot_energy_check.py --filter cdca_001 glca_001 lca3s_001 lca_001
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

BASE_DIR = "/scratch/alpine/ivta1597/LCA_boltz_models"
TARGET_TEMP_K = 300.0

DIR_TYPE = {
    "Binder": "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry": "neg_fail_gate",
}


def read_seq_list(path, name_filter):
    """Parses a seq_ids.txt-style file into (type_dir, seq_id) pairs.

    Args:
        path (str): Path to the file. Blank lines and lines starting with
            "#" are skipped.
        name_filter (list[str] | None): If given, only keep seq_ids
            containing any of these substrings.

    Returns:
        list[tuple]: (type_dir, seq_id) pairs, with type_dir mapped via
        DIR_TYPE (or left as-is if not a known type).
    """
    sequences = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            seq_id, seq_type = parts[0], parts[1]
            if name_filter and not any(sub in seq_id for sub in name_filter):
                continue
            sequences.append((DIR_TYPE.get(seq_type, seq_type), seq_id))
    return sequences


def read_xvg(path):
    """Reads a GROMACS .xvg file, skipping header lines (# and @).

    Args:
        path (str): Path to the .xvg file.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (x, y) columns 1 and 2.
    """
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
    """Computes mean/std/drift over the second half of a series (post-plateau).

    drift_frac is the fitted trend's total change across the window, as a
    fraction of the mean, catching a slow monotonic trend that a low std
    alone could miss.

    Args:
        x: X values (e.g. time).
        y: Y values (e.g. temperature, density).

    Returns:
        tuple[float, float, float]: (mean, std, drift_frac) over the
        second half of `y`.
    """
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
    """Classifies NVT temperature stability against TARGET_TEMP_K.

    Args:
        mean (float): Second-half mean temperature (K).
        std (float): Second-half temperature std (K).

    Returns:
        tuple[str, str]: ("stable" | "UNSTABLE", message). Unstable if the
        mean is off target by >5 K or std exceeds 5 K.
    """
    if abs(mean - TARGET_TEMP_K) > 5 or std > 5:
        return "UNSTABLE", f"mean={mean:.1f} K (target {TARGET_TEMP_K:.0f} K), std={std:.2f} K"
    return "stable", f"mean={mean:.1f} K, std={std:.2f} K"


def verdict_plateau(label, mean, std, drift_frac, unit):
    """Classifies NPT density/volume stability from noise and drift.

    Checks relative std (noise) and total fitted drift (trend) over the
    second half of the series; either one alone can miss instability the
    other catches.

    Args:
        label (str): Quantity name, used in the message (e.g. "density").
        mean (float): Second-half mean.
        std (float): Second-half std.
        drift_frac (float): Second-half drift fraction (see second_half_stats).
        unit (str): Unit string for the message.

    Returns:
        tuple[str, str]: ("stable" | "UNSTABLE", message). Unstable if
        relative std > 2% or drift > 1% of the mean.
    """
    rel_std = std / abs(mean) if mean else float("inf")
    if rel_std > 0.02 or drift_frac > 0.01:
        return "UNSTABLE", (f"{label} mean={mean:.2f} {unit}, std={std:.2f} "
                             f"({rel_std*100:.1f}% of mean), drift={drift_frac*100:.1f}% of mean")
    return "stable", (f"{label} mean={mean:.2f} {unit}, std={std:.2f} "
                       f"({rel_std*100:.2f}% of mean), drift={drift_frac*100:.2f}% of mean")


def parse_args():
    """Parses CLI arguments for the EM/NVT/NPT stability check."""
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seq-list", default="seq_ids.txt",
                    help="seq_ids.txt-style file to read sequences from (default: seq_ids.txt)")
    p.add_argument("--suffix", default="",
                    help='directory/filename suffix, e.g. "_qfix" (default: none, standard EM/NVT/NPT)')
    p.add_argument("--filter", nargs="+", default=None,
                    help="only process seq_ids containing any of these substrings "
                         "(default: all rows in --seq-list)")
    return p.parse_args()


def main():
    """Plots and prints an EM/NVT/NPT stability verdict for each sequence."""
    args = parse_args()
    sequences = read_seq_list(args.seq_list, args.filter)
    suffix = args.suffix

    summary = []
    for type_dir, seq_id in sequences:
        workdir = os.path.join(BASE_DIR, type_dir, seq_id)
        em_f = os.path.join(workdir, f"EM{suffix}", f"em_potential{suffix}.xvg")
        nvt_f = os.path.join(workdir, f"NVT{suffix}", f"nvt_temp{suffix}.xvg")
        dens_f = os.path.join(workdir, f"NPT{suffix}", f"npt_density{suffix}.xvg")
        vol_f = os.path.join(workdir, f"NPT{suffix}", f"npt_volume{suffix}.xvg")

        if not all(os.path.isfile(f) for f in [em_f, nvt_f, dens_f, vol_f]):
            print(f"SKIP {seq_id}: missing one or more .xvg files "
                  f"(run `bash run_energy_check.sh {args.seq_list} {suffix}` first)")
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

        fig.suptitle(f"{seq_id}{suffix} -- EM/NVT/NPT stability check")

        out_png = os.path.join(workdir, f"{seq_id}_EM_EQ_energy{suffix}.png")
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

        summary.append((seq_id, overall, t_msg, d_msg, v_msg))

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for seq_id, overall, t_msg, d_msg, v_msg in summary:
        print(f"  {overall:10s} {seq_id}")


if __name__ == "__main__":
    main()
