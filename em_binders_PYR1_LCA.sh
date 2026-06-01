#!/bin/bash

#SBATCH --job-name=em_b_PYR1_LCA
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc3
#SBATCH --partition=amilan
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#export OMPI_MCA_btl="self,openib,vader,tcp"
#export OMPI_MCA_pml="ob1"
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

# Energy minimization - binders
#cd $DIR/binders/seq${ID}_binder
cd $DIR/binders/pair_${ID}_open_binder
mkdir EM
cd EM
#gmx grompp -f $MDP/em.mdp -c $DIR/binders/seq${ID}_binder/HMR/dodecahedron/seq${ID}_b_dodecahedron_HMR.gro -p $DIR/binders/seq${ID}_binder/HMR/dodecahedron/seq${ID}_b_dodecahedron_HMR.top -o em.tpr
#gmx grompp -f $MDP/em.mdp -c $DIR/binders/pair_${ID}_binder/pair${ID}_dodecahedron_HMR.gro -p $DIR/binders/pair_${ID}_binder/pair${ID}_dodecahedron_HMR.top -o em.tpr
gmx grompp -f $MDP/em.mdp -c $DIR/binders/pair_${ID}_open_binder/pair${ID}_open_dodecahedron_HMR.gro -p $DIR/binders/pair_${ID}_open_binder/pair${ID}_open_dodecahedron_HMR.top -o em.tpr
gmx mdrun -deffnm em

