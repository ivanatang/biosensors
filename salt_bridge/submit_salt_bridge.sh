#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../config.yaml"
SEQ_FILE="${SCRIPT_DIR}/../seq_ids.txt"

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

    sbatch --job-name="sb_${seq_id}" \
           "${SCRIPT_DIR}/run_salt_bridge.sh" \
           "$CONFIG" "$seq_id" "$seq_type" "$SCRIPT_DIR"

    echo "Submitted: $seq_id ($label -> $seq_type)"
done < "$SEQ_FILE"
