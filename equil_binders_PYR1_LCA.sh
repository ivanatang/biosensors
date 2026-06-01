#!/bin/bash

#SBATCH --job-name=eq_b_PYR1_LCA
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc3
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
DIR=`pwd`
MDP=$DIR/MDP

# Get sequence value from command line
ID=$1

### Binders - HMR
# dodecahedron unit cell
# NVT
#cd $DIR/binders/seq${ID}_binder # this is for seqID
#cd $DIR/binders/pair_${ID}_binder # pair_ID
cd $DIR/binders/pair_${ID}_open_binder
mkdir NVT
cd NVT
#gmx grompp -f $MDP/nvt.mdp -c $DIR/binders/seq${ID}_binder/EM/em.gro -r $DIR/binders/seq${ID}_binder/EM/em.gro -p $DIR/binders/seq${ID}_binder/seq${ID}_b_dodecahedron_HMR.top -o nvt.tpr
#gmx grompp -f $MDP/nvt.mdp -c $DIR/binders/pair_${ID}_binder/EM/em.gro -r $DIR/binders/pair_${ID}_binder/EM/em.gro -p $DIR/binders/pair_${ID}_binder/pair${ID}_dodecahedron_HMR.top -o nvt.tpr
gmx grompp -f $MDP/nvt.mdp -c $DIR/binders/pair_${ID}_open_binder/EM/em.gro -r $DIR/binders/pair_${ID}_open_binder/EM/em.gro -p $DIR/binders/pair_${ID}_open_binder/pair${ID}_open_dodecahedron_HMR.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
#cd $DIR/binders/seq${ID}_binder
#cd $DIR/binders/pair_${ID}_binder
cd $DIR/binders/pair_${ID}_open_binder
mkdir NPT
cd NPT
#gmx grompp -f $MDP/npt.mdp -c $DIR/binders/seq${ID}_binder/NVT/nvt.gro -t $DIR/binders/seq${ID}_binder/NVT/nvt.cpt -p $DIR/binders/seq${ID}_binder/seq${ID}_b_dodecahedron_HMR.top -r $DIR/binders/seq${ID}_binder/NVT/nvt.gro -o npt.tpr
#gmx grompp -f $MDP/npt.mdp -c $DIR/binders/pair_${ID}_binder/NVT/nvt.gro -t $DIR/binders/pair_${ID}_binder/NVT/nvt.cpt -p $DIR/binders/pair_${ID}_binder/pair${ID}_dodecahedron_HMR.top -r $DIR/binders/pair_${ID}_binder/NVT/nvt.gro -o npt.tpr
gmx grompp -f $MDP/npt.mdp -c $DIR/binders/pair_${ID}_open_binder/NVT/nvt.gro -t $DIR/binders/pair_${ID}_open_binder/NVT/nvt.cpt -p $DIR/binders/pair_${ID}_open_binder/pair${ID}_open_dodecahedron_HMR.top -r $DIR/binders/pair_${ID}_open_binder/NVT/nvt.gro -o npt.tpr
gmx mdrun -deffnm npt

