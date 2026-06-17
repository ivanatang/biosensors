#!/bin/bash

#SBATCH --job-name=em_low_pkt_PYR1_LCA
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
DIR=`pwd`
MDP=$DIR/MDP

# Get sequence value from command line
ID=$1

# Energy minimization - negative low pocket
cd $DIR/neg_low_pkt/pair_${ID}_low_pkt
mkdir EM
cd EM
gmx grompp -f $MDP/em.mdp -c $DIR/neg_low_pkt/pair_${ID}_low_pkt/pair${ID}_dodecahedron_HMR.gro -p $DIR/neg_low_pkt/pair_${ID}_low_pkt/pair${ID}_dodecahedron_HMR.top -o em.tpr
gmx mdrun -deffnm em

