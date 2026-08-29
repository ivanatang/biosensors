#!/bin/bash

#SBATCH --job-name=getRscores
#SBATCH --error=error_%j.err
#SBATCH --output=output_%j.out
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# Usage:
#   sbatch get_agg_r_scores.sh [seq_list] [start_ns] [end_ns]
#
# Arguments:
#   seq_list - path to seq_ids.txt-style sequence list
#              (default: seq_ids_orig.txt in the repo root)
#   start_ns - start of analysis window in ns (default: 40)
#   end_ns   - end of analysis window in ns   (default: 500)
#
# Examples:
#   sbatch get_agg_r_scores.sh                                # seq_ids.txt, full 40–500 ns
#   sbatch get_agg_r_scores.sh seq_ids_missing_water.txt       # custom list, full 40–500 ns
#   sbatch get_agg_r_scores.sh seq_ids.txt 40 250              # 250 ns window
#   sbatch get_agg_r_scores.sh seq_ids.txt 40 300              # 300 ns window
#
# Output files:
#   r_scores_all_sequences_{start}_{end}ns.csv
#   dw_scores_all_sequences_{start}_{end}ns.csv
# ============================================================

set -euo pipefail

SEQ_LIST=${1:-/projects/ivta1597/biosensors/seq_ids_orig.txt}
START_NS=${2:-40}
END_NS=${3:-500}

module purge
module load anaconda
conda activate biosensors

# pandas/scipy (and other compiled-extension) imports need the envs own
# newer libstdc++, not the HPC images old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version GLIBCXX_3.4.29 not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

OUT_DIR=/projects/ivta1597/biosensors/water_analysis

echo "============================================================"
echo "  Aggregating R/D/W scores"
echo "  Seq list : $SEQ_LIST"
echo "  Window   : ${START_NS}–${END_NS} ns"
echo "  Output   : $OUT_DIR"
echo "============================================================"

python aggregate_r_scores.py \
    --seq_list $SEQ_LIST    \
    --out_dir  $OUT_DIR     \
    --start-ns $START_NS    \
    --end-ns   $END_NS

echo ""
echo "Output files written to $OUT_DIR:"
echo "  r_scores_all_sequences_${START_NS}_${END_NS}ns.csv"
echo "  dw_scores_all_sequences_${START_NS}_${END_NS}ns.csv"
