#!/bin/bash

#SBATCH --job-name=getRscores
#SBATCH --output=output_%j.out
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# Usage:
#   sbatch get_agg_r_scores.sh [start_ns] [end_ns]
#
# Arguments:
#   start_ns - start of analysis window in ns (default: 40)
#   end_ns   - end of analysis window in ns   (default: 500)
#
# Examples:
#   sbatch get_agg_r_scores.sh              # full 40–500 ns
#   sbatch get_agg_r_scores.sh 40 250       # 250 ns window
#   sbatch get_agg_r_scores.sh 40 300       # 300 ns window
#
# Output files:
#   r_scores_all_sequences_{start}_{end}ns.csv
#   dw_scores_all_sequences_{start}_{end}ns.csv
# ============================================================

set -euo pipefail

START_NS=${1:-40}
END_NS=${2:-500}

module purge
module load anaconda
conda activate biosensors

OUT_DIR=/scratch/alpine/ivta1597/LCA_boltz_models/water_analysis_flex

echo "============================================================"
echo "  Aggregating R/D/W scores"
echo "  Window   : ${START_NS}–${END_NS} ns"
echo "  Output   : $OUT_DIR"
echo "============================================================"

python aggregate_r_scores.py \
    --seq_list seq_ids.txt  \
    --out_dir  $OUT_DIR     \
    --start-ns $START_NS    \
    --end-ns   $END_NS

echo ""
echo "Output files written to $OUT_DIR:"
echo "  r_scores_all_sequences_${START_NS}_${END_NS}ns.csv"
echo "  dw_scores_all_sequences_${START_NS}_${END_NS}ns.csv"
