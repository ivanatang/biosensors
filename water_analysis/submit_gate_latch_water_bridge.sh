#!/bin/bash
# submit_gate_latch_water_bridge.sh
# ─────────────────────────────────────────────────────────────────────────────
# Reads seq_ids.txt and submits a separate run_gate_latch_water_bridge.sh
# SLURM job for each sequence. Sequences with a custom path (3rd column) are
# skipped since their trajectories live in a non-standard directory -- run
# those manually with:
#   sbatch run_gate_latch_water_bridge.sh <seq_id> <dir_type> [start_ns] [end_ns] [ligand_region]
#
# Usage:
#   bash submit_gate_latch_water_bridge.sh                          # 40-500 ns, core (default)
#   bash submit_gate_latch_water_bridge.sh seq_ids_orig.txt        # specify seq list
#   bash submit_gate_latch_water_bridge.sh seq_ids_orig.txt 40 500 whole  # whole ligand
#
# seq_ids.txt format (tab-separated):
#   seq_id              seq_type (display)    optional_custom_path
#   pair_3069_binder    Binder
# ─────────────────────────────────────────────────────────────────────────────

SEQ_LIST=${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}
START_NS=${2:-40}
END_NS=${3:-500}
LIGAND_REGION=${4:-core}

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

REGION_SUFFIX=""
[ "$LIGAND_REGION" != "whole" ] && REGION_SUFFIX="_${LIGAND_REGION}"

echo "============================================================"
echo "  Gate-latch-ligand water bridge submission"
echo "  Seq list  : $SEQ_LIST"
echo "  Window    : ${START_NS}-${END_NS} ns"
echo "  Region    : ${LIGAND_REGION}"
echo "  Output dir: gate_latch_water_bridge_${START_NS}_${END_NS}ns${REGION_SUFFIX}/"
echo "============================================================"

# ── Map display seq_type -> directory name used in the file system ─────────
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

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path): $seq_id"
        ((skipped++))
        continue
    fi

    dir_type=$(get_dir_type "$seq_type")

    echo "Submitting: $seq_id  [$seq_type -> $dir_type]  window=${START_NS}-${END_NS}ns  region=${LIGAND_REGION}"
    sbatch run_gate_latch_water_bridge.sh "$seq_id" "$dir_type" "$START_NS" "$END_NS" "$LIGAND_REGION"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Skipped   : $skipped sequences (run manually)"
echo ""
echo "  To run skipped sequences manually:"
echo "  sbatch run_gate_latch_water_bridge.sh <seq_id> <dir_type> $START_NS $END_NS $LIGAND_REGION"
