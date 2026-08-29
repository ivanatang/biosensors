"""
extract_gate_latch_rmsd_feats.py

Computes Ca RMSD-to-own-Boltz-reference (protein_<ID>.pdb), independently
for each region, for every sequence in seq_ids.txt. Regions covered: gate
(84-90), latch (114-118), loop Lb7a5 (148-155), the C-terminal recoil helix
(154-166), and the whole protein.

This is a different axis from the existing gate-latch pairwise distance
(which conflates gate and latch motion into one relative coordinate) and
from per-region RMSF (which measures spread around a trajectory's own mean
structure, blind to where that mean sits). Two separate superpositions are
used here:
  - A "core" alignment (protein Ca atoms excluding gate, latch, Lb7a5, and
    the recoil helix) so none of those four regions' own motion can
    contaminate the alignment used to score them. Gate/Latch/Lb7a5/Recoil
    Ca RMSD are all scored against this one alignment.
  - A separate whole-protein alignment (every common Ca atom, no exclusions)
    -- the conventional global-drift metric, answering a different question
    than the four region-specific scores above and so deliberately not
    reusing their alignment.
Both test whether a region's (or the whole structure's) average position
drifts toward or away from its modeled starting pose over the trajectory.

Outputs:
    <run_dir>/gate_rmsd_to_ref.xvg           per-frame gate Ca RMSD to reference
    <run_dir>/latch_rmsd_to_ref.xvg          per-frame latch Ca RMSD to reference
    <run_dir>/lb7a5_rmsd_to_ref.xvg          per-frame Lb7a5 Ca RMSD to reference
    <run_dir>/recoil_rmsd_to_ref.xvg         per-frame recoil helix Ca RMSD to reference
    <run_dir>/whole_protein_rmsd_to_ref.xvg  per-frame whole-protein Ca RMSD to reference
    <REPO_DIR>/analysis/gate_latch_rmsd_to_ref_summary{TAG}.csv   per-sequence summary,
        including early/late window means (default 100 ns) and a
        full-trajectory linear regression slope (A/ns) for every region -
        the recommended ML features are the wide-window late mean and/or
        the regression slope, since those were the most robust
        binder/nonbinder discriminators found for gate/latch in this
        analysis.

Runs on Alpine (not locally), since inputs are read from the PetaLibrary
archive rather than scratch. After running, pull the resulting summary CSV
back to a local checkout (e.g. via dtn.rc.colorado.edu) for
rmsd_to_ref_significance.py / plot_gate_latch_rmsd_to_ref.ipynb.

Usage:
    conda activate biosensors
    python extract_gate_latch_rmsd_feats.py [seq_ids_orig.txt] [--tag _500ns]

    # or via SLURM:
    sbatch submit_extract_rmsd_to_ref.sh [seq_ids_orig.txt] [--tag _500ns]
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import mdtraj as md
from scipy.stats import linregress

# ── Configurable paths / window tag ───────────────────────────────────────────
# Trajectory inputs are read from the PetaLibrary archive, not scratch --
# scratch auto-deletes after 90 days and older runs' xtc/gro/medoid files are
# already gone (mirrors water_analysis/R_score_calc.py and
# LIG_contacts/contact_type_analysis.py's PetaLibrary convention).
BASE     = "/pl/active/shirts_archive/IvanaTang/biosensors"
REPO_DIR = "/projects/ivta1597/biosensors"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
TAG    = "_500ns"
NM_TO_ANG = 10.0
# ─────────────────────────────────────────────────────────────────────────────

TYPE_SUBDIR = {
    "Binder":         "binders",
    "False Positive": "nonbinders",
    "Low Confidence": "neg_low_pkt",
    "Fail Geometry":  "neg_fail_gate",
}

GATE   = (84, 90)
LATCH  = (114, 118)
LB7A5  = (148, 155)
RECOIL = (154, 166)

# Regions scored against the core alignment (below), which excludes all of
# these residues from the atoms used to compute that alignment.
REGIONS = {
    "Gate":   GATE,
    "Latch":  LATCH,
    "Lb7a5":  LB7A5,
    "Recoil": RECOIL,
}

REGION_XVG = {
    "Gate":   "gate_rmsd_to_ref.xvg",
    "Latch":  "latch_rmsd_to_ref.xvg",
    "Lb7a5":  "lb7a5_rmsd_to_ref.xvg",
    "Recoil": "recoil_rmsd_to_ref.xvg",
    "Whole":  "whole_protein_rmsd_to_ref.xvg",
}

EXCLUDED_RESIDS = set()
for _lo, _hi in (GATE, LATCH, LB7A5, RECOIL):
    EXCLUDED_RESIDS.update(range(_lo, _hi + 1))


def seq_run_dir(seq_id, group_label, end_ns):
    """Directory containing a sequence's per-frame trajectory outputs.
    Mirrors extract_rmsf_feats.py's rmsf_run_dir: newer pipeline runs write
    under runrel/{end_ns}ns/, older ones write directly into runrel/."""
    run_dir = os.path.join(BASE, TYPE_SUBDIR[group_label], seq_id, RUNREL)
    nested_dir = os.path.join(run_dir, f"{int(end_ns)}ns")
    if os.path.exists(os.path.join(nested_dir, "medoid_PL.pdb")):
        return nested_dir
    return run_dir


def find_reference_pdb(seq_id, group_label):
    """Per-sequence Boltz-predicted protein-only structure, protein_<ID>.pdb.
    Naming isn't fully uniform across naming conventions (e.g. protein_019.pdb,
    protein_pair3061.pdb, protein_seq10_binder.pdb), so glob rather than
    parse the ID out of seq_id."""
    seq_dir = os.path.join(BASE, TYPE_SUBDIR[group_label], seq_id)
    candidates = sorted(
        p for p in glob.glob(os.path.join(seq_dir, "protein_*.pdb"))
        if "fixed_H" not in os.path.basename(p)
    )
    return candidates[0] if candidates else None


def ca_index_map(topology):
    """resSeq -> atom index, for protein Ca atoms only."""
    return {
        res.resSeq: atom.index
        for res in topology.residues
        if res.is_protein
        for atom in res.atoms
        if atom.name == "CA"
    }


def region_indices(ca_map, lo, hi):
    return [ca_map[r] for r in range(lo, hi + 1) if r in ca_map]


def _rmsd_to_ref(traj_super, traj_idx, ref, ref_idx):
    """Per-frame RMSD (Angstrom) between traj_super's atoms at traj_idx and
    ref's atoms at ref_idx, given a trajectory already superposed onto ref
    using whatever alignment atom set the caller chose."""
    diff = traj_super.xyz[:, traj_idx, :] - ref.xyz[0, ref_idx, :]
    return np.sqrt(np.mean(np.sum(diff ** 2, axis=2), axis=1)) * NM_TO_ANG


def compute_all_rmsd(seq_id, group_label, end_ns):
    ref_pdb = find_reference_pdb(seq_id, group_label)
    run_dir = seq_run_dir(seq_id, group_label, end_ns)
    top_pdb = os.path.join(run_dir, "medoid_PL.pdb")
    xtc     = os.path.join(run_dir, f"PL_only_40_{int(end_ns)}ns.xtc")

    if ref_pdb is None or not os.path.exists(top_pdb) or not os.path.exists(xtc):
        return None

    ref  = md.load(ref_pdb)
    traj = md.load(xtc, top=top_pdb)

    ref_ca  = ca_index_map(ref.topology)
    traj_ca = ca_index_map(traj.topology)

    # ── Core alignment: excludes gate/latch/Lb7a5/recoil so their own
    # motion can't contaminate the alignment used to score them ──
    core_resids = sorted((set(ref_ca) & set(traj_ca)) - EXCLUDED_RESIDS)
    if not core_resids:
        return None
    core_ref_idx  = [ref_ca[r]  for r in core_resids]
    core_traj_idx = [traj_ca[r] for r in core_resids]
    traj_core = traj.superpose(ref, atom_indices=core_traj_idx, ref_atom_indices=core_ref_idx)

    region_rmsd = {}
    for name, (lo, hi) in REGIONS.items():
        ref_idx  = region_indices(ref_ca, lo, hi)
        traj_idx = region_indices(traj_ca, lo, hi)
        if not ref_idx or not traj_idx:
            return None
        region_rmsd[name] = _rmsd_to_ref(traj_core, traj_idx, ref, ref_idx)

    # ── Whole-protein alignment: a SEPARATE superposition on every common Ca
    # atom (no exclusions) -- the conventional global-drift metric, not the
    # core alignment above. mdtraj's Trajectory has no public .copy(), but
    # re-superposing traj_core (== traj; superpose() mutates in place and
    # returns self) is safe here regardless: Kabsch superposition always
    # finds the globally optimal fit for whatever atom_indices are given,
    # independent of the trajectory's current orientation, so re-running it
    # with a different atom selection produces the same result as superposing
    # fresh from the original coordinates. This is only correct because the
    # four core-aligned region_rmsd arrays above were already computed (as
    # independent numpy arrays) before this second, whole-protein-atom
    # superposition mutates traj_core's coordinates again. ──
    whole_resids = sorted(set(ref_ca) & set(traj_ca))
    if not whole_resids:
        return None
    whole_ref_idx  = [ref_ca[r]  for r in whole_resids]
    whole_traj_idx = [traj_ca[r] for r in whole_resids]
    traj_whole = traj_core.superpose(ref, atom_indices=whole_traj_idx, ref_atom_indices=whole_ref_idx)
    region_rmsd["Whole"] = _rmsd_to_ref(traj_whole, whole_traj_idx, ref, whole_ref_idx)

    time_ns = traj_core.time / 1000.0  # mdtraj reports ps
    return time_ns, region_rmsd


def write_xvg(path, time_ns, values, ylabel):
    with open(path, "w") as f:
        f.write('@    xaxis label "Time (ns)"\n')
        f.write(f'@    yaxis label "{ylabel}"\n')
        for t, v in zip(time_ns, values):
            f.write(f"{t:.4f}\t{v:.5f}\n")


def load_seq_ids(seq_list_path):
    """Yields (seq_id, group_label) tuples from a seq_ids.txt-style file."""
    with open(seq_list_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            seq_id      = parts[0].strip()
            group_label = parts[1].strip() if len(parts) > 1 else ""
            yield seq_id, group_label


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list", nargs="?", default="seq_ids.txt")
    parser.add_argument("--tag", default=TAG,
                        help=f"Suffix appended to output CSV filename (default: {TAG})")
    parser.add_argument("--wide-window-ns", type=float, default=100.0,
                        help="Window size (ns) for early/late mean comparison (default: 100.0), "
                             "used for the recommended ML feature.")
    parser.add_argument("--end-ns", type=float, default=500.0,
                        help="End of analysis window in ns; selects which "
                             "windowed subdirectory / PL_only xtc to read "
                             "(default: 500)")
    parser.add_argument("--suffix", default="",
                        help="Run-directory suffix, e.g. '_qfix' for the "
                             "bond-order/charge-fix systems (default: '', the "
                             "standard production directory). Applies to every "
                             "sequence in seq_list -- run _qfix sequences via a "
                             "separate seq_list from standard ones.")
    args = parser.parse_args()
    tag = args.tag

    global RUNREL
    RUNREL = f"{RUNREL}{args.suffix}"
    wide_window = args.wide_window_ns
    end_ns = args.end_ns

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    all_systems = list(load_seq_ids(args.seq_list))

    rows, missing = [], []
    for seq_id, group_label in all_systems:
        if group_label not in TYPE_SUBDIR:
            missing.append(seq_id)
            continue
        result = compute_all_rmsd(seq_id, group_label, end_ns)
        if result is None:
            missing.append(seq_id)
            continue
        time_ns, region_rmsd = result

        run_dir = seq_run_dir(seq_id, group_label, end_ns)
        for name, xvg_name in REGION_XVG.items():
            write_xvg(os.path.join(run_dir, xvg_name), time_ns, region_rmsd[name],
                      f"{name} Ca RMSD to Boltz reference (A)")

        t_end = time_ns[-1]
        early_wide_mask = time_ns <= (time_ns[0] + wide_window)
        late_wide_mask  = time_ns >= (t_end - wide_window)

        row = {
            "Sequence": seq_id,
            "Group": group_label,
            "N frames": len(time_ns),
        }

        for name in REGION_XVG:
            rmsd = region_rmsd[name]
            row[f"{name} RMSD mean (A)"] = round(float(rmsd.mean()), 4)
            row[f"{name} RMSD SD (A)"]   = round(float(rmsd.std(ddof=1)), 4)

            early_wide = float(rmsd[early_wide_mask].mean())
            late_wide  = float(rmsd[late_wide_mask].mean())
            slope      = float(linregress(time_ns, rmsd).slope)
            row[f"{name} RMSD early{int(wide_window)} mean (A)"] = round(early_wide, 4)
            row[f"{name} RMSD late{int(wide_window)} mean (A)"]  = round(late_wide, 4)
            row[f"{name} drift{int(wide_window)} (A)"]          = round(late_wide - early_wide, 4)
            row[f"{name} slope (A/ns)"]                          = round(slope, 6)

        rows.append(row)

    if missing:
        print(f"Skipped {len(missing)} sequences (missing reference/topology/trajectory): "
              f"{missing[:5]}{'...' if len(missing) > 5 else ''}")

    if not rows:
        print("\nNo sequences processed - nothing to write.")
        sys.exit(1)

    df = pd.DataFrame(rows).set_index("Sequence")
    out_dir = os.path.join(REPO_DIR, "analysis")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"gate_latch_rmsd_to_ref_summary{tag}.csv")
    df.to_csv(csv_path)
    print(f"Saved ({len(df)} rows) -> {csv_path}")


if __name__ == "__main__":
    main()
