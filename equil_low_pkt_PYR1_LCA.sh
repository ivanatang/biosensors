#!/bin/bash

#SBATCH --job-name=eq_low_pkt_PYR1_LCA
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

### Binders - HMR
# dodecahedron unit cell
# NVT
#cd $DIR/neg_low_pkt/pair_${ID}_low_pkt
cd $DIR/neg_low_pkt/pair_${ID}_open_low_pkt
mkdir NVT
cd NVT
#gmx grompp -f $MDP/nvt.mdp -c $DIR/neg_low_pkt/pair_${ID}_low_pkt/EM/em.gro -r $DIR/neg_low_pkt/pair_${ID}_low_pkt/EM/em.gro -p $DIR/neg_low_pkt/pair_${ID}_low_pkt/pair${ID}_dodecahedron_HMR.top -o nvt.tpr
gmx grompp -f $MDP/nvt.mdp -c $DIR/neg_low_pkt/pair_${ID}_open_low_pkt/EM/em.gro -r $DIR/neg_low_pkt/pair_${ID}_open_low_pkt/EM/em.gro -p $DIR/neg_low_pkt/pair_${ID}_open_low_pkt/pair${ID}_open_dodecahedron_HMR.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
#cd $DIR/neg_low_pkt/pair_${ID}_low_pkt
cd $DIR/neg_low_pkt/pair_${ID}_open_low_pkt
mkdir NPT
cd NPT
#gmx grompp -f $MDP/npt.mdp -c $DIR/neg_low_pkt/pair_${ID}_low_pkt/NVT/nvt.gro -t $DIR/neg_low_pkt/pair_${ID}_low_pkt/NVT/nvt.cpt -p $DIR/neg_low_pkt/pair_${ID}_low_pkt/pair${ID}_dodecahedron_HMR.top -r $DIR/neg_low_pkt/pair_${ID}_low_pkt/NVT/nvt.gro -o npt.tpr
gmx grompp -f $MDP/npt.mdp -c $DIR/neg_low_pkt/pair_${ID}_open_low_pkt/NVT/nvt.gro -t $DIR/neg_low_pkt/pair_${ID}_open_low_pkt/NVT/nvt.cpt -p $DIR/neg_low_pkt/pair_${ID}_open_low_pkt/pair${ID}_open_dodecahedron_HMR.top -r $DIR/neg_low_pkt/pair_${ID}_open_low_pkt/NVT/nvt.gro -o npt.tpr
gmx mdrun -deffnm npt

