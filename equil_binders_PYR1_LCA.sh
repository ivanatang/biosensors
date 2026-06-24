#!/bin/bash

#SBATCH --job-name=eq_b_PYR1_LCA
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
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

# Set some environment variables
DIR=/projects/ivta1597/biosensors
MDP=$DIR/MDP

# Get sequence value from command line
ID=$1
PREFIX=$2

### Binders - HMR
# dodecahedron unit cell
# NVT
cd $DIR/binders/${PREFIX}_${ID}_binder
mkdir NVT
cd NVT
gmx grompp -f $MDP/nvt.mdp -c $DIR/binders/${PREFIX}_${ID}_binder/EM/em.gro -r $DIR/binders/${PREFIX}_${ID}_binder/EM/em.gro -p $DIR/binders/${PREFIX}_${ID}_binder/${PREFIX}_${ID}_dodecahedron_HMR.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
cd $DIR/binders/${PREFIX}_${ID}_binder
mkdir NPT
cd NPT
gmx grompp -f $MDP/npt.mdp -c $DIR/binders/${PREFIX}_${ID}_binder/NVT/nvt.gro -t $DIR/binders/${PREFIX}_${ID}_binder/NVT/nvt.cpt -p $DIR/binders/${PREFIX}_${ID}_binder/${PREFIX}_${ID}_dodecahedron_HMR.top -r $DIR/binders/${PREFIX}_${ID}_binder/NVT/nvt.gro -o npt.tpr
gmx mdrun -deffnm npt

