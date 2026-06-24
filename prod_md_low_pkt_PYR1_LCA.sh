#!/bin/bash

#SBATCH --job-name=1413_lp
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc3
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export TMPDIR=$SLURM_SCRATCH
export SLURM_EXPORT_ENV=ALL

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
PME=16
RDD=1.2

D1=6
D2=4
D3=2

# Production simulation - negative low pocket
cd $DIR/neg_low_pkt/pair_${ID}_low_pkt
mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/neg_low_pkt/pair_${ID}_low_pkt/NPT/npt.gro -t $DIR/neg_low_pkt/pair_${ID}_low_pkt/NPT/npt.cpt -p $DIR/neg_low_pkt/pair_${ID}_low_pkt/pair${ID}_dodecahedron_HMR.top -o prod_md_500ns.tpr
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3

