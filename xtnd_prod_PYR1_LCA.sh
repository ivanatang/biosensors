#!/bin/bash

#SBATCH --job-name=032xt_nb
#SBATCH --output=out_xt_nb_%j.out                  # Output file
#SBATCH --error=err_xt_nb_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Usage:
#   sbatch xtnd_prod_PYR1_LCA.sh <ID> <SEQ_TYPE> <PREFIX> 

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

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

cd "$SEQ_DIR"
cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -s prod_md_500ns.tpr -cpi prod_md_500ns.cpt -append -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3

