#!/usr/bin/env bash
# =============================================================================
# submit_mdpocket_exploration.sh  —  submission script
#
# Usage:
#   bash submit_mdpocket_exploration.sh [seq_ids.txt]
# =============================================================================

SEQ_LIST=${1:-seq_ids.txt}

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"       ;;
        "False Positive") echo "nonbinders"    ;;
        "Low Confidence") echo "neg_low_pkt"   ;;
        "Fail Geometry")  echo "neg_fail_gate" ;;
        *)                echo "$1"            ;;
    esac
}

submitted=0; skipped=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path): $seq_id"
        ((skipped++)); continue
    fi

    dir_type=$(get_dir_type "$seq_type")
    echo "Submitting: $seq_id  [$seq_type → $dir_type]"
    sbatch pkt_vol_exp.sh "$seq_id" "$dir_type"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "============================="
echo " Submitted : $submitted"
echo " Skipped   : $skipped"
echo "============================="
