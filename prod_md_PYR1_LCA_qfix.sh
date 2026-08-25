#!/bin/bash

#SBATCH --job-name=prod_qfix
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# prod_md_PYR1_LCA_qfix.sh -- initial production chunk for the _qfix pilot
# systems. Reads NPT_qfix/npt.gro (from equil_PYR1_LCA_qfix.sh) and the
# _qfix topology. Writes to a _qfix-suffixed run directory so the original
# production trajectory for this sequence is never touched.
#
# For the pilot: submit this alone first and see how far one 24h chunk gets
# (chain xtnd_prod_PYR1_LCA_qfix.sh only if you want to extend further --
# the pilot doesn't need the full 500ns to compare against the original).
#
# Usage:
#   sbatch prod_md_PYR1_LCA_qfix.sh <ID> <SEQ_TYPE> <PREFIX>

export TMPDIR=$SLURM_SCRATCH
export SLURM_EXPORT_ENV=ALL

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate biosensors

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

PME=16
RDD=1.2

D1=6
D2=4
D3=2

# Navigate to the sequence's scratch directory for trajectory files
SEQ_DIR=$BASE/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}
if [ ! -d "$SEQ_DIR" ]; then
    echo "ERROR: Sequence directory not found: $SEQ_DIR" >&2
    exit 1
fi

TOP=$SEQ_DIR/${PREFIX}_${ID}_dodecahedron_HMR_qfix.top
if [ ! -f "$SEQ_DIR/NPT_qfix/npt.gro" ]; then
    echo "ERROR: NPT_qfix/npt.gro not found -- run equil_PYR1_LCA_qfix.sh first" >&2
    exit 1
fi

PROD_DIR=prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd_qfix

cd "$SEQ_DIR"
mkdir -p "$PROD_DIR"
cd "$PROD_DIR"
gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $SEQ_DIR/NPT_qfix/npt.gro -t $SEQ_DIR/NPT_qfix/npt.cpt -p "$TOP" -o prod_md_500ns.tpr
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3
