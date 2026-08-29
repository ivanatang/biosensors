#!/usr/bin/env bash
# =============================================================================
# submit_mdpocket_char.sh  —  submission script
#
# Loops through SEQ_LIST and submits one SLURM job per sequence to run
# mdpocket characterization using the selected_pocket.pdb from exploration.
#
# Usage:
#   bash submit_mdpocket_char.sh [seq_ids_orig.txt] [--overwrite-existing] [--suffix _qfix]
#
# --overwrite-existing is passed through to every submitted job, forcing
# mdpocket characterization to rerun instead of skipping sequences that
# already have descriptors.
# --suffix is passed through too, e.g. "_qfix" for the bond-order/charge-fix
# systems (default: "", the standard production directory).
# =============================================================================

# Absolute path, not relative -- sbatch resolves a bare filename against the
# directory this script was invoked from, not this script's own location.
WORKER="/projects/ivta1597/biosensors/pkt_vol/pkt_vol_char.sh"

OVERWRITE=false
SUFFIX=""
SEQ_LIST="/projects/ivta1597/biosensors/seq_ids_orig.txt"
next_is_suffix=false
for arg in "$@"; do
    if [[ "$next_is_suffix" == "true" ]]; then
        SUFFIX="$arg"
        next_is_suffix=false
        continue
    fi
    case "$arg" in
        --overwrite-existing) OVERWRITE=true ;;
        --suffix)              next_is_suffix=true ;;
        *)                    SEQ_LIST="$arg" ;;
    esac
done

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

# ── Map display seq_type → directory name used in the file system ─────────────
get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"       ;;
        "False Positive") echo "nonbinders"    ;;
        "Low Confidence") echo "neg_low_pkt"   ;;
        "Fail Geometry")  echo "neg_fail_gate" ;;
        *)                echo "$1"            ;;
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

    echo "Submitting: $seq_id  [$seq_type → $dir_type]"
    sbatch "$WORKER" "$seq_id" "$dir_type" "$OVERWRITE" "$SUFFIX"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "============================="
echo " Submitted : $submitted"
echo " Skipped   : $skipped"
echo "============================="
