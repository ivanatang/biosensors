#!/bin/bash

#SBATCH --job-name=em_nb_PYR1_LCA
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc3
#SBATCH --partition=amilan
#SBATCH --time=00:20:00
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

conda activate IS_env

# Set some environment variables 
DIR=`pwd`
MDP=$DIR/MDP

# Get sequence value from command line
ID=$1

# Energy minimization - nonbinders
# dodecahedron
#cd $DIR/nonbinders/seq${ID}_nb
#cd $DIR/nonbinders/pair_${ID}_nb
cd $DIR/nonbinders/pair_${ID}_open_nb
mkdir EM
cd EM
#gmx grompp -f $MDP/em.mdp -c $DIR/nonbinders/seq${ID}_nb/HMR/dodecahedron/seq${ID}_nb_dodecahedron_HMR.gro -p $DIR/nonbinders/seq${ID}_nb/HMR/dodecahedron/seq${ID}_nb_dodecahedron_HMR.top -o em.tpr
#gmx grompp -f $MDP/em.mdp -c $DIR/nonbinders/pair_${ID}_nb/pair${ID}_dodecahedron_HMR.gro -p $DIR/nonbinders/pair_${ID}_nb/pair${ID}_dodecahedron_HMR.top -o em.tpr
gmx grompp -f $MDP/em.mdp -c $DIR/nonbinders/pair_${ID}_open_nb/pair${ID}_open_dodecahedron_HMR.gro -p $DIR/nonbinders/pair_${ID}_open_nb/pair${ID}_open_dodecahedron_HMR.top -o em.tpr
gmx mdrun -deffnm em

# cube
#cd $DIR/nonbinders/${SEQ}_nb/HMR
#mkdir EM
#cd EM
#gmx grompp -f $MDP/em.mdp -c $DIR/nonbinders/${SEQ}_nb/HMR/${SEQ}_nb_HMR.gro -p $DIR/nonbinders/${SEQ}_nb/HMR/${SEQ}_nb_HMR.top -o em.tpr
#gmx mdrun -deffnm em

