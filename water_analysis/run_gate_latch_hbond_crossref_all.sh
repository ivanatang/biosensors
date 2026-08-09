#!/bin/bash

#SBATCH --job-name=gl_hbond_crossref
#SBATCH --output=output_glhbondxref_%j.out
#SBATCH --error=error_glhbondxref_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# run_gate_latch_hbond_crossref_all.sh
# ------------------------------------------------------------
# Runs parse_gate_latch_hbond.py for every sequence in seq_ids.txt, one
# process each, sequentially. Each call is pure CSV/xvg parsing (no
# trajectory I/O), so unlike gate_latch_hbond_gmx.sh this does NOT need one
# SLURM job per sequence -- a single short job loops over all of them.
#
# Must run AFTER, for every sequence:
#   1. gate_latch_water_bridge.py --dump-frames  (submit_gate_latch_water_bridge.sh ... true)
#   2. gate_latch_hbond_gmx.sh                    (submit_gate_latch_hbond.sh)
# Sequences missing either input are skipped with a warning, not a hard
# failure, so one bad run doesn't block the other 193.
#
# Usage:
#   sbatch run_gate_latch_hbond_crossref_all.sh [seq_list] [hbond_start_ns] [hbond_end_ns]
#
# Example:
#   sbatch run_gate_latch_hbond_crossref_all.sh seq_ids_orig.txt 40 500
# ============================================================

set -uo pipefail   # NOT -e: one sequence's failure shouldn't kill the loop

module purge
module load anaconda
conda activate biosensors

SCRIPT_DIR="/projects/ivta1597/biosensors/water_analysis"
SEQ_LIST=${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}
START_NS=${2:-40}
END_NS=${3:-500}

get_dir_type() {
    case "$1" in
        "Binder")         echo "binders"      ;;
        "False Positive") echo "nonbinders"   ;;
        "Low Confidence") echo "neg_low_pkt"  ;;
        "Fail Geometry")  echo "neg_fail_gate";;
        *)                echo "$1"           ;;
    esac
}

ok=0; failed=0; skipped=0

while IFS=$'\t' read -r seq_id seq_type custom_path || [[ -n "$seq_id" ]]; do
    [[ -z "$seq_id" || "$seq_id" == \#* ]] && continue
    if [[ -n "$custom_path" ]]; then
        echo "SKIP (custom path): $seq_id"
        ((skipped++))
        continue
    fi

    dir_type=$(get_dir_type "$seq_type")
    echo "── $seq_id [$dir_type] ──"

    python "${SCRIPT_DIR}/parse_gate_latch_hbond.py" \
        --seq_id "$seq_id" --seq_type "$dir_type" \
        --start-ns "$START_NS" --end-ns "$END_NS"
    py_status=$?

    # parse_gate_latch_hbond.py can exit 0 while still writing an EMPTY
    # crossref CSV (e.g. all 4 gmx hbond .xvg files missing for this
    # sequence -- it warns per-group and just moves on rather than
    # crashing). A zero exit code alone doesn't mean this sequence's
    # result is usable, so also check the output file actually has content.
    crossref_csv="/scratch/alpine/ivta1597/LCA_boltz_models/${dir_type}/${seq_id}/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/${seq_id}_gate_latch_hbond_crossref_$(printf '%d' "$START_NS")_$(printf '%d' "$END_NS")ns.csv"

    if [[ $py_status -eq 0 && -s "$crossref_csv" ]]; then
        ((ok++))
    else
        echo "  FAILED: $seq_id (exit=${py_status}, output=$([[ -s "$crossref_csv" ]] && echo present || echo missing/empty))"
        ((failed++))
    fi

done < "$SEQ_LIST"

echo ""
echo "=== Done ==="
echo "  OK      : $ok"
echo "  Failed  : $failed"
echo "  Skipped : $skipped"
echo ""
echo "Next: python agg_gate_latch_hbond.py --start-ns $START_NS --end-ns $END_NS"
