#!/bin/bash

#SBATCH --job-name=em_PYR1_LCA_qfix
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# em_PYR1_LCA_qfix.sh -- EM for the bond-order/charge-fix ("_qfix") pilot
# systems. Identical to em_PYR1_LCA.sh except it reads the _qfix topology
# (deprotonated LCA carboxylate, corrected counterion count) and writes to
# EM_qfix/ so the original production EM output is never touched.
#
# Usage:
#   sbatch em_PYR1_LCA_qfix.sh <ID> <SEQ_TYPE> <PREFIX>

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate biosensors

# pandas/scipy (and other compiled-extension) imports need the envs own
# newer libstdc++, not the HPC images old /lib64/libstdc++.so.6 --
# ImportError: .../libstdc++.so.6: version GLIBCXX_3.4.29 not found
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

# Set some environment variables
DIR=/projects/ivta1597/biosensors
MDP=$DIR/MDP
BASE=/scratch/alpine/ivta1597/LCA_boltz_models

# Get sequence value from command line
ID=$1
SEQ_TYPE=$2 # binders | nonbinders | neg_fail_gate | neg_low_pkt
PREFIX=$3 # pair, bind

if [ "$SEQ_TYPE" == "binders" ]; then
    SUFFIX="binder"
elif [ "$SEQ_TYPE" == "nonbinders" ]; then
    SUFFIX="nb"
elif [ "$SEQ_TYPE" == "neg_fail_gate" ]; then
    SUFFIX="fail_gate"
elif [ "$SEQ_TYPE" == "neg_low_pkt" ]; then
    SUFFIX="low_pkt"
else
    echo "ERROR: Unknown SEQ_TYPE '$SEQ_TYPE'" >&2
    exit 1
fi

# #SBATCH directives are static text parsed before ID/SEQ_TYPE/PREFIX are
# known, so relabel the job now that they are, rather than leaving every
# submission of this script showing the same generic name in squeue.
scontrol update JobId=$SLURM_JOB_ID JobName="${PREFIX}_${ID}_em_qfix" 2>/dev/null

# Navigate to the sequence's scratch directory for trajectory files
SEQ_DIR=$BASE/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}
if [ ! -d "$SEQ_DIR" ]; then
    echo "ERROR: Sequence directory not found: $SEQ_DIR" >&2
    exit 1
fi

TOP=$SEQ_DIR/${PREFIX}_${ID}_${SUFFIX}_dodecahedron_HMR_qfix.top
GRO=$SEQ_DIR/${PREFIX}_${ID}_${SUFFIX}_dodecahedron_HMR_qfix.gro
if [ ! -f "$TOP" ] || [ ! -f "$GRO" ]; then
    echo "ERROR: _qfix topology not found ($TOP / $GRO) -- sync it from OneDrive first" >&2
    exit 1
fi

cd "$SEQ_DIR"
mkdir -p EM_qfix
cd EM_qfix
gmx grompp -f $MDP/em.mdp -c "$GRO" -p "$TOP" -o em.tpr
gmx mdrun -deffnm em
