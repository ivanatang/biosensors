#!/bin/bash

#SBATCH --job-name=eq_nb_PYR1_LCA
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

### Nonbinders
# dodecahedron
# NVT
#cd $DIR/nonbinders/seq${ID}_nb
#cd $DIR/nonbinders/pair_${ID}_nb
cd $DIR/nonbinders/pair_${ID}_open_nb
mkdir NVT
cd NVT
#gmx grompp -f $MDP/nvt.mdp -c $DIR/nonbinders/seq${ID}_nb/EM/em.gro -r $DIR/nonbinders/seq${ID}_nb/EM/em.gro -p $DIR/nonbinders/seq${ID}_nb/seq${ID}_nb_dodecahedron_HMR.top -o nvt.tpr
#gmx grompp -f $MDP/nvt.mdp -c $DIR/nonbinders/pair_${ID}_nb/EM/em.gro -r $DIR/nonbinders/pair_${ID}_nb/EM/em.gro -p $DIR/nonbinders/pair_${ID}_nb/pair${ID}_dodecahedron_HMR.top -o nvt.tpr
gmx grompp -f $MDP/nvt.mdp -c $DIR/nonbinders/pair_${ID}_open_nb/EM/em.gro -r $DIR/nonbinders/pair_${ID}_open_nb/EM/em.gro -p $DIR/nonbinders/pair_${ID}_open_nb/pair${ID}_open_dodecahedron_HMR.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
#cd $DIR/nonbinders/seq${ID}_nb
#cd $DIR/nonbinders/pair_${ID}_nb
cd $DIR/nonbinders/pair_${ID}_open_nb
mkdir NPT
cd NPT
#gmx grompp -f $MDP/npt.mdp -c $DIR/nonbinders/seq${ID}_nb/NVT/nvt.gro -t $DIR/nonbinders/seq${ID}_nb/NVT/nvt.cpt -p $DIR/nonbinders/seq${ID}_nb/seq${ID}_nb_dodecahedron_HMR.top -r $DIR/nonbinders/seq${ID}_nb/NVT/nvt.gro -o npt.tpr
#gmx grompp -f $MDP/npt.mdp -c $DIR/nonbinders/pair_${ID}_nb/NVT/nvt.gro -t $DIR/nonbinders/pair_${ID}_nb/NVT/nvt.cpt -p $DIR/nonbinders/pair_${ID}_nb/pair${ID}_dodecahedron_HMR.top -r $DIR/nonbinders/pair_${ID}_nb/NVT/nvt.gro -o npt.tpr
gmx grompp -f $MDP/npt.mdp -c $DIR/nonbinders/pair_${ID}_open_nb/NVT/nvt.gro -t $DIR/nonbinders/pair_${ID}_open_nb/NVT/nvt.cpt -p $DIR/nonbinders/pair_${ID}_open_nb/pair${ID}_open_dodecahedron_HMR.top -r $DIR/nonbinders/pair_${ID}_open_nb/NVT/nvt.gro -o npt.tpr
gmx mdrun -deffnm npt

# cube
# NVT
#cd $DIR/nonbinders/${SEQ}_nb/HMR
#mkdir NVT
#cd NVT
#gmx grompp -f $MDP/nvt.mdp -c $DIR/nonbinders/${SEQ}_nb/HMR/EM/em.gro -r $DIR/nonbinders/${SEQ}_nb/HMR/EM/em.gro -p $DIR/nonbinders/${SEQ}_nb/HMR/${SEQ}_nb_HMR.top -o nvt.tpr
#gmx mdrun -deffnm nvt

# NPT
#cd $DIR/nonbinders/${SEQ}_nb/HMR
#mkdir NPT
#cd NPT
#gmx grompp -f $MDP/npt.mdp -c $DIR/nonbinders/${SEQ}_nb/HMR/NVT/nvt.gro -t $DIR/nonbinders/${SEQ}_nb/HMR/NVT/nvt.cpt -p $DIR/nonbinders/${SEQ}_nb/HMR/${SEQ}_nb_HMR.top -r $DIR/nonbinders/${SEQ}_nb/HMR/NVT/nvt.gro -o npt.tpr
#gmx mdrun -deffnm npt
