#!/bin/bash
#SBATCH --job-name=rg_sasa
#SBATCH --output=logs/output_%A_%a.out
#SBATCH --error=logs/error_%A_%a.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-128%20

# Usage: sbatch --array=1-<N>%20 submit_compute_Rg_sasa.sh [--region pocket|whole] [seq_ids.txt]
#   N = number of lines in seq_ids.txt  (wc -l seq_ids.txt)

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
    NDX_SELECTION="0"
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

# ── Load modules ──────────────────────────────────────────────────────────────
module purge
module load gcc
module load openmpi

# ── Get this task's line from seq_ids.txt ─────────────────────────────────────
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SEQ_IDS")
if [[ -z "$line" ]]; then
    echo "No entry for array task ${SLURM_ARRAY_TASK_ID} — exiting"
    exit 0
fi

folder_name=$(echo "$line" | cut -f1)
label=$(echo "$line"       | cut -f2)
custom_base=$(echo "$line" | cut -f3)

echo "Task ${SLURM_ARRAY_TASK_ID}: $folder_name ($label)  region: $REGION"

# ── Resolve run directory ─────────────────────────────────────────────────────
subdir="${TYPE_SUBDIR[$label]:-}"
if [[ -z "$subdir" ]]; then
    echo "ERROR: unknown label '$label' — exiting"
    exit 1
fi

if [[ -n "$custom_base" ]]; then
    rundir="${custom_base}/${RUNREL}"
else
    rundir="${BASE}/${subdir}/${folder_name}/${RUNREL}"
fi

echo "rundir: $rundir"

if [[ ! -d "$rundir" ]]; then
    echo "ERROR: directory not found — $rundir"
    exit 1
fi

# ── Check required inputs ─────────────────────────────────────────────────────
missing=""
for f in "${REQUIRED_INPUTS[@]}"; do
    [[ ! -f "${rundir}/${f}" ]] && missing+=" $f"
done
if [[ -n "$missing" ]]; then
    echo "ERROR: missing inputs:$missing"
    exit 1
fi

# ── Build pocket index when needed ────────────────────────────────────────────
if [[ "$REGION" == "pocket" ]]; then
    ndx="${rundir}/pocket_res.ndx"
    if [[ ! -f "$ndx" ]]; then
        echo "Building pocket residue index..."
        "$GMX" select \
            -s "${rundir}/medoid_PL.pdb" \
            -on "$ndx" \
            -select "protein and resid $POCKET_RESIDS" \
            2>> "${rundir}/${RG_LOG}"
    fi
    if [[ ! -f "$ndx" ]]; then
        echo "ERROR: failed to create pocket_res.ndx (check ${RG_LOG})"
        exit 1
    fi
    ndx_args=("-n" "$ndx")
else
    ndx_args=("-n" "${rundir}/index.ndx")
fi

# ── Rg ────────────────────────────────────────────────────────────────────────
if [[ -f "${rundir}/${RG_OUT}" ]]; then
    echo "SKIP: $RG_OUT already exists"
else
    echo "Running gmx gyrate..."
    echo "$NDX_SELECTION" | "$GMX" gyrate \
        -s "${rundir}/medoid_PL.pdb" \
        -f "${rundir}/PL_only_40_500ns.xtc" \
        "${ndx_args[@]}" \
        -o "${rundir}/${RG_OUT}" \
        2>> "${rundir}/${RG_LOG}" || true
    if [[ -f "${rundir}/${RG_OUT}" ]]; then
        echo "  → $RG_OUT done"
    else
        echo "  → ERROR: $RG_OUT not created (check $RG_LOG)"
    fi
fi

# ── SASA ──────────────────────────────────────────────────────────────────────
if [[ -f "${rundir}/${SASA_OUT}" ]]; then
    echo "SKIP: $SASA_OUT already exists"
else
    echo "Running gmx sasa..."
    echo "$NDX_SELECTION" | "$GMX" sasa \
        -s "${rundir}/medoid_PL.pdb" \
        -f "${rundir}/PL_only_40_500ns.xtc" \
        "${ndx_args[@]}" \
        -o "${rundir}/${SASA_OUT}" \
        2>> "${rundir}/${SASA_LOG}" || true
    if [[ -f "${rundir}/${SASA_OUT}" ]]; then
        echo "  → $SASA_OUT done"
    else
        echo "  → ERROR: $SASA_OUT not created (check $SASA_LOG)"
    fi
fi

echo "Task ${SLURM_ARRAY_TASK_ID} complete"
