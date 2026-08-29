#!/bin/bash
# Usage: bash submit_salt_bridge.sh [seq_list] [start_ns] [end_ns] [config]
#   bash submit_salt_bridge.sh                                    # seq_ids.txt, 40-500 ns, config.yaml
#   bash submit_salt_bridge.sh ../seq_ids_orig.txt 40 250          # 250 ns window, full cohort
#   bash submit_salt_bridge.sh ../seq_ids_new_ligands.txt 40 500 ../config_qfix.yaml
#                                                                  # _qfix systems (differently-named runrel)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEQ_FILE="${1:-${SCRIPT_DIR}/../seq_ids.txt}"
START_NS="${2:-40}"
END_NS="${3:-500}"
CONFIG="${4:-${SCRIPT_DIR}/../config.yaml}"

FAILED_LIST="${SCRIPT_DIR}/submit_salt_bridge_failed_${START_NS}_${END_NS}ns.txt"
> "$FAILED_LIST"

submitted=0
failed=0

# Map feat_table group labels -> seq_type keys used in config.yaml's type_subdir
map_seq_type() {
    case "$1" in
        "Binder")          echo "binder" ;;
        "False Positive")  echo "nb" ;;
        "Low Confidence")  echo "low_pkt" ;;
        "Fail Geometry")   echo "fail_gate" ;;
        *)                 echo "" ;;
    esac
}

while IFS=$'\t' read -r seq_id label; do
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    seq_type=$(map_seq_type "$label")
    if [[ -z "$seq_type" ]]; then
        echo "WARNING: unrecognized label '$label' for $seq_id — skipping"
        continue
    fi

    sbatch_out=$(sbatch --job-name="sb_${seq_id}" \
           --output="${SCRIPT_DIR}/output_sb_${seq_id}_%j.out" \
           --error="${SCRIPT_DIR}/error_sb_${seq_id}_%j.err" \
           "${SCRIPT_DIR}/run_salt_bridge.sh" \
           "$CONFIG" "$seq_id" "$seq_type" "$SCRIPT_DIR" "$START_NS" "$END_NS" 2>&1)
    if [[ $? -ne 0 ]]; then
        echo "SUBMIT FAILED: $seq_id  -- $sbatch_out"
        printf '%s\t%s\n' "$seq_id" "$label" >> "$FAILED_LIST"
        ((failed++))
        continue
    fi

    echo "Submitted: $seq_id ($label -> $seq_type)  window=${START_NS}-${END_NS}ns  ($sbatch_out)"
    ((submitted++))
done < "$SEQ_FILE"

echo ""
echo "=== Done: submitted ${submitted} jobs, ${failed} failed ==="
if (( failed > 0 )); then
    echo "  $failed failed submissions written to: $FAILED_LIST"
    echo "  Re-submit just those: bash ${SCRIPT_DIR}/submit_salt_bridge.sh $FAILED_LIST $START_NS $END_NS"
fi
