"""
select_binding_pocket.py

Identifies the LCA binding site pocket from mdpocket exploration output for
each sequence in a seq_ids.txt file, by finding alpha spheres in the frequency
isosurface PDB within a distance cutoff of the ligand centroid.

Outputs a selected_pocket.pdb per sequence for use in mdpocket characterization:
    mdpocket --trajectory_file protein_only.xtc \
             --trajectory_format xtc \
             --selected_pocket selected_pocket.pdb \
             -f protein_only_ref.pdb \
             -o mdpocket_<seq_id>

Usage:
    python select_binding_pocket.py [seq_ids.txt] [--cutoff 8.0] [--ligand_resname LIG]
"""

import os
import sys
import argparse
import numpy as np

# ── Configurable paths ────────────────────────────────────────────────────────
BASE   = "/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
# ─────────────────────────────────────────────────────────────────────────────


def get_dir_type(seq_type):
    mapping = {
        "Binder":         "binders",
        "False Positive": "nonbinders",
        "Low Confidence": "neg_low_pkt",
        "Fail Geometry":  "neg_fail_gate",
    }
    return mapping.get(seq_type, seq_type)


def parse_pdb_atoms(pdb_path, record_types=("ATOM", "HETATM")):
    """Return (n,3) coordinate array and list of raw lines."""
    coords = []
    lines  = []
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() in record_types:
                try:
                    coords.append([float(line[30:38]),
                                   float(line[38:46]),
                                   float(line[46:54])])
                    lines.append(line)
                except ValueError:
                    continue
    return np.array(coords) if coords else np.empty((0, 3)), lines


def get_ligand_coords(pdb_path, ligand_resname):
    """Extract coordinates of all ligand atoms matching resname."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() in ("ATOM", "HETATM"):
                if line[17:20].strip() == ligand_resname:
                    try:
                        coords.append([float(line[30:38]),
                                       float(line[38:46]),
                                       float(line[46:54])])
                    except ValueError:
                        continue
    if not coords:
        raise ValueError(
            f"No atoms with resname '{ligand_resname}' found in {pdb_path}.\n"
            f"  Check with: grep HETATM {pdb_path} | awk '{{print $4}}' | sort -u"
        )
    return np.array(coords)


def process_sequence(seq_id, run_dir, cutoff, ligand_resname):
    """Run pocket selection for a single sequence. Returns 'ok', 'skip', or 'fail'."""

    pl_pdb   = os.path.join(run_dir, "medoid_PL.pdb")
    freq_pdb = os.path.join(run_dir, f"mdpocket_{seq_id}_freq_iso_0_5.pdb")
    out_pdb  = os.path.join(run_dir, "selected_pocket.pdb")

    # ── Validate inputs ───────────────────────────────────────────────────────
    if not os.path.isdir(run_dir):
        print(f"  SKIP: run directory not found: {run_dir}")
        return "skip"

    for path, label in [(pl_pdb, "medoid_PL.pdb"), (freq_pdb, "freq_iso PDB")]:
        if not os.path.exists(path):
            print(f"  SKIP: {label} not found: {path}")
            return "skip"

    # ── Skip if already done ──────────────────────────────────────────────────
    if os.path.exists(out_pdb):
        print(f"  SKIP: selected_pocket.pdb already exists")
        return "skip"

    # ── Ligand centroid ───────────────────────────────────────────────────────
    try:
        lig_coords = get_ligand_coords(pl_pdb, ligand_resname)
    except ValueError as e:
        print(f"  FAIL: {e}")
        return "fail"

    centroid = lig_coords.mean(axis=0)
    print(f"  Ligand atoms    : {len(lig_coords)}")
    print(f"  Ligand centroid : ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f}) Å")

    # ── Load alpha spheres ────────────────────────────────────────────────────
    pocket_coords, pocket_lines = parse_pdb_atoms(freq_pdb)
    if len(pocket_coords) == 0:
        print(f"  FAIL: no atoms found in {freq_pdb}")
        return "fail"

    print(f"  Alpha spheres   : {len(pocket_coords)} total")

    # ── Distance filter ───────────────────────────────────────────────────────
    dists = np.linalg.norm(pocket_coords - centroid, axis=1)
    mask  = dists <= cutoff
    n_kept = mask.sum()

    print(f"  Cutoff          : {cutoff} Å  →  {n_kept} spheres kept")

    if n_kept == 0:
        nearest = dists.min()
        print(f"  WARN: nearest sphere is {nearest:.2f} Å away — "
              f"try --cutoff {int(np.ceil(nearest + 2))}")
        return "fail"

    # ── Write selected pocket ─────────────────────────────────────────────────
    with open(out_pdb, "w") as f:
        for line in np.array(pocket_lines)[mask]:
            f.write(line)
        f.write("END\n")

    print(f"  → {out_pdb}")
    return "ok"


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq_list",        nargs="?", default="seq_ids.txt",
                        help="Tab-separated seq_ids file (default: seq_ids.txt)")
    parser.add_argument("--cutoff",        type=float, default=8.0,
                        help="Distance cutoff in Å from ligand centroid (default: 8.0)")
    parser.add_argument("--ligand_resname",default="LIG",
                        help="Residue name of the ligand in medoid_PL.pdb (default: LIG)")
    args = parser.parse_args()

    if not os.path.exists(args.seq_list):
        print(f"ERROR: seq list not found: {args.seq_list}")
        sys.exit(1)

    counts = {"ok": 0, "skip": 0, "fail": 0}

    with open(args.seq_list) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            seq_id      = parts[0].strip()
            seq_type    = parts[1].strip() if len(parts) > 1 else ""
            custom_path = parts[2].strip() if len(parts) > 2 else ""

            print(f"\n{seq_id}  [{seq_type}]")

            if custom_path:
                run_dir = os.path.join(custom_path, RUNREL)
            else:
                run_dir = os.path.join(BASE, get_dir_type(seq_type), seq_id, RUNREL)

            result = process_sequence(seq_id, run_dir, args.cutoff, args.ligand_resname)
            counts[result] += 1

    print(f"\n{'='*35}")
    print(f" Processed : {counts['ok']}")
    print(f" Skipped   : {counts['skip']}")
    print(f" Failed    : {counts['fail']}")
    print(f"{'='*35}")


if __name__ == "__main__":
    main()
