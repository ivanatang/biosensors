#!/bin/bash
# submit_multi_ID.sh
# ─────────────────────────────────────────────────────────────────────────────
# Reads seq_ids.txt and submits a separate run_gate_latch.sh SLURM job
# for each sequence. The full WORKDIR path is constructed from BASE_DIR,
# the directory type, the sequence ID, and RUN_SUBDIR.
#
# Sequences with a custom path (3rd column) are skipped since their
# trajectories live in a non-standard directory — run those manually with:
#   sbatch run_gate_latch.sh <full_workdir_path>
#
# Usage:
#   bash submit_multi_ID.sh                   # uses seq_ids.txt by default
#   bash submit_multi_ID.sh my_seq_list.txt   # pass a different file
#
# seq_ids.txt format (tab-separated):
#   seq_id              seq_type (display)    optional_custom_path
#   pair_3069_binder    Binder
#   seq14_binder        Binder                /scratch/.../water_contacts   ← skipped
# ─────────────────────────────────────────────────────────────────────────────

SEQ_LIST=${1:-seq_ids.txt}

# ── Path configuration ────────────────────────────────────────────────────────
BASE_DIR="/scratch/alpine/ivta1597/LCA_boltz_models"
RUN_SUBDIR="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
# ─────────────────────────────────────────────────────────────────────────────

if [ ! -f "$SEQ_LIST" ]; then
    echo "ERROR: seq list file not found: $SEQ_LIST"
    exit 1
fi

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
missing=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do

    # Skip empty lines and comments
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue

    # Skip sequences with a custom path — their trajectories are in a
    # non-standard location that run_gate_latch.sh does not handle
    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path): $seq_id"
        ((skipped++))
        continue
    fi

    dir_type=$(get_dir_type "$seq_type")

    WORKDIR="${BASE_DIR}/${dir_type}/${seq_id}/${RUN_SUBDIR}"

    # Warn if the directory does not exist, but still submit so SLURM captures
    # the error in the job log rather than silently dropping the sequence
    if [[ ! -d "$WORKDIR" ]]; then
        echo "WARNING: Directory not found (will fail at runtime): $WORKDIR"
        ((missing++))
    fi

    echo "Submitting: $seq_id  [$seq_type → $dir_type]"
    echo "           WORKDIR: $WORKDIR"
    sbatch run_50pct_analysis.sh "$WORKDIR"
    ((submitted++))

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  Submitted : $submitted jobs"
echo "  Skipped   : $skipped sequences (custom path — run manually)"
if (( missing > 0 )); then
    echo "  WARNING   : $missing submitted jobs have a missing WORKDIR — check error logs"
fi
