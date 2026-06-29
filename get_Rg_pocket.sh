#!/bin/bash
# get_Rg_pocket.sh — Compute Rg for pocket residues only (serial, one process per sequence)
# Pocket: 27 consensus positions from the PYR1/LCA binding pocket (Leonard et al. 2024)
# Usage: bash get_Rg_pocket.sh [seq_ids.txt]

set -uo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
SEQ_IDS="${1:-seq_ids.txt}"

# 27 consensus pocket positions (from config.yaml: medoid.pocket_resids)
POCKET_RESIDS="59 60 61 62 79 81 83 87 88 89 91 92 94 108 109 110 115 117 120 122 141 158 159 160 163 164 167"

declare -A TYPE_SUBDIR=(
    ["Binder"]="binders"
    ["False Positive"]="nonbinders"
    ["Low Confidence"]="neg_low_pkt"
    ["Fail Geometry"]="neg_fail_gate"
)

# ── Counters ──────────────────────────────────────────────────────────────────
skip=0; run=0; err=0

# ── Main loop ─────────────────────────────────────────────────────────────────
while IFS=$'\t' read -r folder_name label custom_base; do
    [[ -z "$folder_name" || "$folder_name" == \#* ]] && continue

    subdir="${TYPE_SUBDIR[$label]:-}"
    if [[ -z "$subdir" ]]; then
        echo "WARNING: unknown label '$label' for $folder_name — skipping"
        continue
    fi

    if [[ -n "$custom_base" ]]; then
        rundir="${custom_base}/${RUNREL}"
    else
        rundir="${BASE}/${subdir}/${folder_name}/${RUNREL}"
    fi

    if [[ ! -d "$rundir" ]]; then
        echo "WARNING: directory not found — $rundir"
        continue
    fi

    if [[ -f "${rundir}/Rg_pocket.xvg" ]]; then
        echo "SKIP: $folder_name (Rg_pocket.xvg already exists)"
        skip=$(( skip + 1 ))
        continue
    fi

    echo "── $folder_name ──"

    missing=""
    for f in medoid_PL.pdb PL_only_40_500ns.xtc; do
        [[ ! -f "${rundir}/${f}" ]] && missing+=" $f"
    done
    if [[ -n "$missing" ]]; then
        echo "  ERROR: missing inputs:$missing"
        err=$(( err + 1 ))
        continue
    fi

    (
        cd "$rundir"

        # Build a pocket-residue index (all atoms of the 27 pocket positions)
        "$GMX" select \
            -s medoid_PL.pdb \
            -on pocket_res.ndx \
            -select "protein and resid $POCKET_RESIDS" \
            2>> gyrate_pocket.log

        if [[ ! -f "pocket_res.ndx" ]]; then
            echo "  ERROR: failed to create pocket_res.ndx (check gyrate_pocket.log)"
            exit 1
        fi

        # Group 0 is the single selection written by gmx select
        echo "0" | "$GMX" gyrate \
            -s medoid_PL.pdb \
            -f PL_only_40_500ns.xtc \
            -n pocket_res.ndx \
            -o Rg_pocket.xvg \
            2>> gyrate_pocket.log || true

        if [[ -f "Rg_pocket.xvg" ]]; then
            echo "  → Rg_pocket.xvg done"
        else
            echo "  → ERROR: Rg_pocket.xvg not created (check gyrate_pocket.log)"
            exit 1
        fi
    )

    if [[ -f "${rundir}/Rg_pocket.xvg" ]]; then
        run=$(( run + 1 ))
    else
        err=$(( err + 1 ))
    fi

done < "$SEQ_IDS"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════"
printf "%-12s  %8s  %8s  %8s\n" ""            "Skipped" "Generated" "Errors"
printf "%-12s  %8d  %8d  %8d\n" "Rg_pocket"  "$skip"    "$run"      "$err"
echo "════════════════════════════════════════"
