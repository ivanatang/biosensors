#!/bin/bash
# submit_water_spatial_run.sh
# ─────────────────────────────────────────────────────────────────────────────
# Reads seq_ids.txt and submits a separate water_spatial_run.sh SLURM job for
# each sequence. Requires water_spatial_prep.sh to have already been run for
# these sequences (produces <run_dir>/water_spatial/{protein_lig.ndx,
# fit_trim.xtc,fit_trim_ref.pdb}) -- water_spatial_prep.sh is a single serial
# loop script (like pkt_vol/pkt_vol_prep.sh), not something submitted here,
# since PBC/fit trjconv calls are I/O-bound rather than the heavy step.
#
# Sequences with a custom path (3rd column) are skipped, same convention as
# submit_water_analysis.sh / submit_pkt_*.sh -- run those manually with:
#   sbatch water_spatial_run.sh <seq_id> <dir_type> [bin_nm]
#
# Usage:
#   bash submit_water_spatial_run.sh                      # seq_ids_ngs_observed.txt
#   bash submit_water_spatial_run.sh seq_ids.txt           # specify seq list
#   bash submit_water_spatial_run.sh seq_ids.txt 0.05      # explicit bin width
#   bash submit_water_spatial_run.sh seq_ids.txt 0.05 300 _qfix
#                                                           # _qfix systems (must match
#                                                           # what water_spatial_prep.sh
#                                                           # was run with)
# ─────────────────────────────────────────────────────────────────────────────

SEQ_LIST=${1:-/projects/ivta1597/biosensors/seq_ids_ngs_observed.txt}
BIN_NM=${2:-0.05}
NAB=${3:-300}
SUFFIX=${4:-}

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

echo "============================================================"
echo "  water_spatial_run.sh submission"
echo "  Seq list     : $SEQ_LIST"
echo "  Bin (nm)     : $BIN_NM"
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

    echo "Submitting: $seq_id  [$seq_type -> $dir_type]"
    sbatch water_spatial_run.sh "$seq_id" "$dir_type" "$BIN_NM" "$NAB" "$SUFFIX"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Skipped   : $skipped sequences (run manually)"
echo ""
echo "  To run skipped sequences manually:"
echo "  sbatch water_spatial_run.sh <seq_id> <dir_type> $BIN_NM"
