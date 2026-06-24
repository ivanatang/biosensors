#!/bin/bash
# gen_rg_sasa.sh — Generate missing Rg_PL.xvg and sasa_PL.xvg for all sequences
# Usage: bash gen_rg_sasa.sh [seq_ids.txt]

set -uo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
GMX="/projects/ivta1597/pkgs/gromacs-2025.3/bin/gmx"
BASE="/scratch/alpine/ivta1597/LCA_boltz_models"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
SEQ_IDS="${1:-seq_ids.txt}"

declare -A TYPE_SUBDIR=(
    ["Binder"]="binders"
    ["False Positive"]="nonbinders"
    ["Low Confidence"]="neg_low_pkt"
    ["Fail Geometry"]="neg_fail_gate"
)

# ── Counters (per-analysis) ───────────────────────────────────────────────────
rg_skip=0;   rg_run=0;   rg_err=0
sasa_skip=0; sasa_run=0; sasa_err=0

# ── Helper: run one gmx command if output is missing ─────────────────────────
# Usage: run_if_missing <outfile> <logfile> <label> <gmx_args...>
run_if_missing() {
    local outfile="$1"; shift
    local logfile="$1"; shift
    local label="$1";   shift
    local skip_var="${label}_skip"
    local run_var="${label}_run"
    local err_var="${label}_err"

    if [[ -f "$outfile" ]]; then
        echo "  SKIP $label: $(basename "$outfile") already exists"
        printf -v "$skip_var" '%d' "$(( ${!skip_var} + 1 ))"
        return
    fi

    echo "  RUN  $label"
    ( echo "Protein_LIG" | "$GMX" "$@" 2>> "$logfile" ) || true

    if [[ -f "$outfile" ]]; then
        echo "  → $label done"
        printf -v "$run_var" '%d' "$(( ${!run_var} + 1 ))"
    else
        echo "  → ERROR: $(basename "$outfile") not created (check $logfile)"
        printf -v "$err_var" '%d' "$(( ${!err_var} + 1 ))"
    fi
}

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

    # Skip entirely only if both outputs already exist
    if [[ -f "${rundir}/Rg_PL.xvg" && -f "${rundir}/sasa_PL.xvg" ]]; then
        echo "SKIP: $folder_name (both outputs exist)"
        rg_skip=$(( rg_skip + 1 ))
        sasa_skip=$(( sasa_skip + 1 ))
        continue
    fi

    echo "── $folder_name ──"

    # Check required inputs
    missing=""
    for f in medoid_PL.pdb PL_only_40_500ns.xtc index.ndx; do
        [[ ! -f "${rundir}/${f}" ]] && missing+=" $f"
    done
    if [[ -n "$missing" ]]; then
        echo "  ERROR: missing inputs:$missing"
        rg_err=$(( rg_err + 1 ))
        sasa_err=$(( sasa_err + 1 ))
        continue
    fi

    (
        cd "$rundir"

        run_if_missing "Rg_PL.xvg" "gyrate_rg.log" "rg" \
            gyrate \
            -s medoid_PL.pdb \
            -f PL_only_40_500ns.xtc \
            -n index.ndx \
            -o Rg_PL.xvg

        run_if_missing "sasa_PL.xvg" "sasa.log" "sasa" \
            sasa \
            -s medoid_PL.pdb \
            -f PL_only_40_500ns.xtc \
            -n index.ndx \
            -o sasa_PL.xvg
    )

done < "$SEQ_IDS"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════"
printf "%-6s  %8s  %8s  %8s\n" ""        "Skipped" "Generated" "Errors"
printf "%-6s  %8d  %8d  %8d\n" "Rg"      "$rg_skip"   "$rg_run"   "$rg_err"
printf "%-6s  %8d  %8d  %8d\n" "SASA"    "$sasa_skip" "$sasa_run" "$sasa_err"
echo "════════════════════════════════════════"
