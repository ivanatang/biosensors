#!/bin/bash
# submit_water_analysis.sh
# ─────────────────────────────────────────────────────────────────────────────
# Reads seq_ids.txt and submits a separate run_water_analysis.sh SLURM job
# for each sequence. Sequences with a custom path (3rd column) are skipped
# since their trajectories live in a non-standard directory — run those
# manually with: sbatch run_water_analysis.sh <seq_id> <dir_type> [start_ns] [end_ns]
#
# Usage:
#   bash submit_water_analysis.sh                          # 40–500 ns (default)
#   bash submit_water_analysis.sh seq_ids.txt              # specify seq list
#   bash submit_water_analysis.sh seq_ids.txt 40 250       # 250 ns window
#   bash submit_water_analysis.sh seq_ids.txt 40 300       # 300 ns window
#   bash submit_water_analysis.sh seq_ids.txt 40 500 core  # steroid core only
#   bash submit_water_analysis.sh seq_ids.txt 40 500 tail  # carboxylate tail only
#   bash submit_water_analysis.sh seq_ids_qfix.txt 40 500 whole _qfix
#                                                           # _qfix systems
#
# NOTE: ligand_region=core/tail only runs R_score_calc.py (steps 2-3 of
# run_water_analysis.sh are whole-ligand only and are skipped automatically).
#
# seq_ids.txt format (tab-separated):
#   seq_id              seq_type (display)    optional_custom_path
#   pair_3069_binder    Binder
#   seq14_binder        Binder                /scratch/.../water_contacts   ← skipped
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SEQ_LIST=${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}
START_NS=${2:-40}
END_NS=${3:-500}
LIGAND_REGION=${4:-whole}
SUFFIX=${5:-}

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

FAILED_LIST="${SCRIPT_DIR}/submit_water_analysis_failed_${START_NS}_${END_NS}ns_${LIGAND_REGION}.txt"
> "$FAILED_LIST"

REGION_SUFFIX=""
[ "$LIGAND_REGION" != "whole" ] && REGION_SUFFIX="_${LIGAND_REGION}"

echo "============================================================"
echo "  Water analysis submission"
echo "  Seq list  : $SEQ_LIST"
echo "  Window    : ${START_NS}–${END_NS} ns"
echo "  Region    : ${LIGAND_REGION}"
echo "  Output dir: water_contacts_${START_NS}_${END_NS}ns${REGION_SUFFIX}/"
echo "============================================================"

# ── Map display seq_type → directory name used in the file system ─────────────
get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"      ;;
        "False Positive") echo "nonbinders"   ;;
        "Low Confidence") echo "neg_low_pkt"  ;;
        "Fail Geometry")  echo "neg_fail_gate";;
        *)                echo "$1"           ;;   # fallback: use as-is
    esac
}

submitted=0
skipped=0
failed=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    # Skip empty lines and comments
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    # Skip sequences with a custom path — their trajectories are in a
    # non-standard location that run_water_analysis.sh does not handle
    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path): $seq_id"
        ((skipped++))
        continue
    fi

    dir_type=$(get_dir_type "$seq_type")

    # Absolute path to the worker + a unique job name/log per sequence:
    # sbatch was previously invoked with a bare relative "run_water_analysis.sh"
    # (only resolves if this script happens to be run from inside water_analysis/)
    # and every job shared the same generic job-name/output_water_%j.out log, so
    # a failed submission was silent and, even when it wasn't, undiagnosable.
    sbatch_out=$(sbatch --job-name="wc_${seq_id}" \
                        --output="${SCRIPT_DIR}/output_water_${seq_id}_%j.out" \
                        --error="${SCRIPT_DIR}/error_water_${seq_id}_%j.err" \
                        "${SCRIPT_DIR}/run_water_analysis.sh" \
                        "$seq_id" "$dir_type" "$START_NS" "$END_NS" "$LIGAND_REGION" "$SUFFIX" 2>&1)
    if [[ $? -ne 0 ]]; then
        echo "SUBMIT FAILED: $seq_id  -- $sbatch_out"
        printf '%s\t%s\n' "$seq_id" "$seq_type" >> "$FAILED_LIST"
        ((failed++))
        continue
    fi

    echo "Submitting: $seq_id  [$seq_type → $dir_type]  window=${START_NS}–${END_NS}ns  region=${LIGAND_REGION}  ($sbatch_out)"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Skipped   : $skipped sequences (custom path, run manually)"
echo "  Failed    : $failed sequences (sbatch itself rejected the submission)"
echo ""
echo "  To run skipped sequences manually:"
echo "  sbatch ${SCRIPT_DIR}/run_water_analysis.sh <seq_id> <dir_type> $START_NS $END_NS $LIGAND_REGION"
if (( failed > 0 )); then
    echo ""
    echo "  $failed failed submissions written to: $FAILED_LIST"
    echo "  Re-submit just those once resolved (e.g. after a job-count cap frees up):"
    echo "  bash ${SCRIPT_DIR}/submit_water_analysis.sh $FAILED_LIST $START_NS $END_NS $LIGAND_REGION"
fi
