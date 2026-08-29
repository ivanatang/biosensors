#!/bin/bash
# submit_hydration_analysis.sh
# ─────────────────────────────────────────────────────────────────────────────
# Reads seq_ids.txt and submits a separate run_hydration_analysis.sh SLURM job
# for each sequence. Sequences with a custom path (3rd column) are skipped
# since their trajectories live in a non-standard directory — run those
# manually with: sbatch run_hydration_analysis.sh <seq_id> <dir_type> [start_ns] [end_ns]
#
# Usage:
#   bash submit_hydration_analysis.sh                              # 40–500 ns, ligand region (default)
#   bash submit_hydration_analysis.sh seq_ids.txt                  # specify seq list
#   bash submit_hydration_analysis.sh seq_ids.txt pocket_residues  # pocket-residue region
#   bash submit_hydration_analysis.sh seq_ids.txt ligand 40 250    # 250 ns window
#   bash submit_hydration_analysis.sh seq_ids.txt ligand 40 500 3,4,5,6,8 10   # explicit cutoffs/stride
#
# seq_ids.txt format (tab-separated):
#   seq_id              seq_type (display)    optional_custom_path
#   pair_3069_binder    Binder
#   seq14_binder        Binder                /scratch/.../HMR/dodecahedron   ← skipped
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SEQ_LIST=${1:-/projects/ivta1597/biosensors/seq_ids_ngs_observed.txt}
REFERENCE_REGION=${2:-ligand}
START_NS=${3:-40}
END_NS=${4:-500}
CUTOFFS=${5:-3,4,5,6,8}
STRIDE=${6:-10}
SUFFIX=${7:-}

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

FAILED_LIST="${SCRIPT_DIR}/submit_hydration_analysis_failed_${START_NS}_${END_NS}ns.txt"
> "$FAILED_LIST"

echo "============================================================"
echo "  Hydration-shell analysis submission"
echo "  Seq list : $SEQ_LIST"
echo "  Region   : ${REFERENCE_REGION}"
echo "  Window   : ${START_NS}–${END_NS} ns"
echo "  Cutoffs  : ${CUTOFFS} A"
echo "  Stride   : ${STRIDE}"
echo "  Out dir  : water_density_${START_NS}_${END_NS}ns/"
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
    # non-standard location that run_hydration_analysis.sh does not handle
    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path): $seq_id"
        ((skipped++))
        continue
    fi

    dir_type=$(get_dir_type "$seq_type")

    sbatch_out=$(sbatch --job-name="hyd_${seq_id}" \
                        --output="${SCRIPT_DIR}/output_hyd_${seq_id}_%j.out" \
                        --error="${SCRIPT_DIR}/error_hyd_${seq_id}_%j.err" \
                        "${SCRIPT_DIR}/run_hydration_analysis.sh" \
                        "$seq_id" "$dir_type" "$REFERENCE_REGION" "$START_NS" "$END_NS" "$CUTOFFS" "$STRIDE" "$SUFFIX" 2>&1)
    if [[ $? -ne 0 ]]; then
        echo "SUBMIT FAILED: $seq_id  -- $sbatch_out"
        printf '%s\t%s\n' "$seq_id" "$seq_type" >> "$FAILED_LIST"
        ((failed++))
        continue
    fi

    echo "Submitting: $seq_id  [$seq_type → $dir_type]  region=${REFERENCE_REGION}  window=${START_NS}–${END_NS}ns  ($sbatch_out)"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Skipped   : $skipped sequences (run manually)"
echo "  Failed    : $failed sequences (sbatch itself rejected the submission)"
echo ""
echo "  To run skipped sequences manually:"
echo "  sbatch ${SCRIPT_DIR}/run_hydration_analysis.sh <seq_id> <dir_type> $REFERENCE_REGION $START_NS $END_NS $CUTOFFS $STRIDE"
if (( failed > 0 )); then
    echo ""
    echo "  $failed failed submissions written to: $FAILED_LIST"
    echo "  Re-submit just those: bash ${SCRIPT_DIR}/submit_hydration_analysis.sh $FAILED_LIST $REFERENCE_REGION $START_NS $END_NS $CUTOFFS $STRIDE"
fi
