#!/bin/bash

#SBATCH --job-name=em_PYR1_LCA
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

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

cd $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}
mkdir EM
cd EM
gmx grompp -f $MDP/em.mdp -c $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/${PREFIX}_${ID}_dodecahedron_HMR.gro -p $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/${PREFIX}_${ID}_dodecahedron_HMR.top -o em.tpr
gmx mdrun -deffnm em

