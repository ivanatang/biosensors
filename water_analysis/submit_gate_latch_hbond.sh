#!/bin/bash
# submit_gate_latch_hbond.sh
# ─────────────────────────────────────────────────────────────────────────────
# Reads seq_ids.txt and submits a separate gate_latch_hbond_gmx.sh SLURM job
# for each sequence. Sequences with a custom path (3rd column) are skipped,
# same convention as submit_gate_latch_water_bridge.sh -- run those manually:
#   sbatch gate_latch_hbond_gmx.sh <seq_id> <dir_type> [start_ns] [end_ns]
#
# Usage:
#   bash submit_gate_latch_hbond.sh                    # 40-500 ns (default)
#   bash submit_gate_latch_hbond.sh seq_ids_orig.txt 40 500
# ─────────────────────────────────────────────────────────────────────────────

RUN_SCRIPT="/projects/ivta1597/biosensors/water_analysis/gate_latch_hbond_gmx.sh"

SEQ_LIST=${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}
START_NS=${2:-40}
END_NS=${3:-500}

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

echo "============================================================"
echo "  Gate-latch backbone/side-chain H-bond (gmx hbond) submission"
echo "  Seq list  : $SEQ_LIST"
echo "  Window    : ${START_NS}-${END_NS} ns"
echo "============================================================"

get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"      ;;
        "False Positive") echo "nonbinders"   ;;
        "Low Confidence") echo "neg_low_pkt"  ;;
        "Fail Geometry")  echo "neg_fail_gate";;
        *)                echo "$1"           ;;
    esac
}

submitted=0
skipped=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path): $seq_id"
        ((skipped++))
        continue
    fi

    dir_type=$(get_dir_type "$seq_type")

    echo "Submitting: $seq_id  [$seq_type -> $dir_type]  window=${START_NS}-${END_NS}ns"
    sbatch "$RUN_SCRIPT" "$seq_id" "$dir_type" "$START_NS" "$END_NS"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Skipped   : $skipped sequences (run manually)"
