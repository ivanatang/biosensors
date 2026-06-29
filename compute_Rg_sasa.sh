#!/bin/bash
# compute_Rg_sasa.sh — Compute Rg and SASA for a specified protein region (serial)
# Usage: bash compute_Rg_sasa.sh [--region pocket|whole] [seq_ids.txt]
#
#   pocket: 27 consensus pocket positions (Leonard et al. 2024)
#           outputs Rg_pocket.xvg, sasa_pocket.xvg
#   whole:  Protein + ligand (Protein_LIG group in index.ndx)
#           outputs Rg_PL.xvg, sasa_PL.xvg

set -uo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
POCKET_RESIDS="59 60 61 62 79 81 83 87 88 89 91 92 94 108 109 110 115 117 120 122 141 158 159 160 163 164 167"
REGION="whole"
SEQ_IDS=""

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --region)   REGION="$2"; shift 2 ;;
        --region=*) REGION="${1#*=}"; shift ;;
        *)          SEQ_IDS="$1"; shift ;;
    esac
done
SEQ_IDS="${SEQ_IDS:-seq_ids.txt}"

if [[ "$REGION" != "whole" && "$REGION" != "pocket" ]]; then
    echo "ERROR: --region must be 'whole' or 'pocket'" >&2
    exit 1
fi

# ── Region-specific config ────────────────────────────────────────────────────
if [[ "$REGION" == "pocket" ]]; then
    SUFFIX="pocket"
    NDX_SELECTION="0"          # group 0 is the single selection written by gmx select
    REQUIRED_INPUTS=("medoid_PL.pdb" "PL_only_40_500ns.xtc")
else
    SUFFIX="PL"
    NDX_SELECTION="Protein_LIG"
    REQUIRED_INPUTS=("medoid_PL.pdb" "PL_only_40_500ns.xtc" "index.ndx")
fi

RG_OUT="Rg_${SUFFIX}.xvg"
SASA_OUT="sasa_${SUFFIX}.xvg"
RG_LOG="gyrate_${SUFFIX}.log"
SASA_LOG="sasa_${SUFFIX}.log"

declare -A TYPE_SUBDIR=(
    ["Binder"]="binders"
    ["False Positive"]="nonbinders"
    ["Low Confidence"]="neg_low_pkt"
    ["Fail Geometry"]="neg_fail_gate"
)

echo "Region: $REGION  |  seq_ids: $SEQ_IDS"
echo ""

# ── Counters ──────────────────────────────────────────────────────────────────
rg_skip=0;   rg_run=0;   rg_err=0
sasa_skip=0; sasa_run=0; sasa_err=0

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

    if [[ -f "${rundir}/${RG_OUT}" && -f "${rundir}/${SASA_OUT}" ]]; then
        echo "SKIP: $folder_name (both outputs exist)"
        rg_skip=$(( rg_skip + 1 ))
        sasa_skip=$(( sasa_skip + 1 ))
        continue
    fi

    echo "── $folder_name ──"

    missing=""
    for f in "${REQUIRED_INPUTS[@]}"; do
        [[ ! -f "${rundir}/${f}" ]] && missing+=" $f"
    done
    if [[ -n "$missing" ]]; then
        echo "  ERROR: missing inputs:$missing"
        rg_err=$(( rg_err + 1 ))
        sasa_err=$(( sasa_err + 1 ))
        continue
    fi

    # Build pocket index when needed (reuse if it already exists)
    if [[ "$REGION" == "pocket" ]]; then
        ndx="${rundir}/pocket_res.ndx"
        if [[ ! -f "$ndx" ]]; then
            echo "  Building pocket residue index..."
            "$GMX" select \
                -s "${rundir}/medoid_PL.pdb" \
                -on "$ndx" \
                -select "protein and resid $POCKET_RESIDS" \
                2>> "${rundir}/${RG_LOG}"
        fi
        if [[ ! -f "$ndx" ]]; then
            echo "  ERROR: failed to create pocket_res.ndx (check ${RG_LOG})"
            rg_err=$(( rg_err + 1 ))
            sasa_err=$(( sasa_err + 1 ))
            continue
        fi
        ndx_args=("-n" "$ndx")
    else
        ndx_args=("-n" "${rundir}/index.ndx")
    fi

    # ── Rg ────────────────────────────────────────────────────────────────────
    if [[ -f "${rundir}/${RG_OUT}" ]]; then
        echo "  SKIP rg: $RG_OUT already exists"
        rg_skip=$(( rg_skip + 1 ))
    else
        echo "  RUN  rg"
        echo "$NDX_SELECTION" | "$GMX" gyrate \
            -s "${rundir}/medoid_PL.pdb" \
            -f "${rundir}/PL_only_40_500ns.xtc" \
            "${ndx_args[@]}" \
            -o "${rundir}/${RG_OUT}" \
            2>> "${rundir}/${RG_LOG}" || true
        if [[ -f "${rundir}/${RG_OUT}" ]]; then
            echo "  → rg done"
            rg_run=$(( rg_run + 1 ))
        else
            echo "  → ERROR: $RG_OUT not created (check $RG_LOG)"
            rg_err=$(( rg_err + 1 ))
        fi
    fi

    # ── SASA ──────────────────────────────────────────────────────────────────
    if [[ -f "${rundir}/${SASA_OUT}" ]]; then
        echo "  SKIP sasa: $SASA_OUT already exists"
        sasa_skip=$(( sasa_skip + 1 ))
    else
        echo "  RUN  sasa"
        echo "$NDX_SELECTION" | "$GMX" sasa \
            -s "${rundir}/medoid_PL.pdb" \
            -f "${rundir}/PL_only_40_500ns.xtc" \
            "${ndx_args[@]}" \
            -o "${rundir}/${SASA_OUT}" \
            2>> "${rundir}/${SASA_LOG}" || true
        if [[ -f "${rundir}/${SASA_OUT}" ]]; then
            echo "  → sasa done"
            sasa_run=$(( sasa_run + 1 ))
        else
            echo "  → ERROR: $SASA_OUT not created (check $SASA_LOG)"
            sasa_err=$(( sasa_err + 1 ))
        fi
    fi

done < "$SEQ_IDS"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════"
printf "%-6s  %8s  %8s  %8s\n" ""       "Skipped" "Generated" "Errors"
printf "%-6s  %8d  %8d  %8d\n" "Rg"     "$rg_skip"   "$rg_run"   "$rg_err"
printf "%-6s  %8d  %8d  %8d\n" "SASA"   "$sasa_skip" "$sasa_run" "$sasa_err"
echo "════════════════════════════════════════"
