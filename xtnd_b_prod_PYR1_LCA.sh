#!/bin/bash

#SBATCH --job-name=3092_open
#SBATCH --output=output_xtnd_prod_md_b_%j.out                  # Output file
#SBATCH --error=error_xtnd_prod_md_b_%j.err                    # Error file
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
PME=16
RDD=1.2

D1=6
D2=4
D3=2

# Production simulation
# dodecahedron
#cd $DIR/binders/${ID}_binder
#cd $DIR/binders/pair_${ID}_binder
cd $DIR/binders/pair_${ID}_open_binder
cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
#cd prod_md_0p9_cutoff_3dt_64x1
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -s prod_md_500ns.tpr -cpi prod_md_500ns.cpt -append -ntomp $SLURM_CPUS_PER_TASK #-rdd $RDD -npme $PME
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -s prod_md_500ns.tpr -cpi prod_md_500ns.cpt -append -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3

# cube
#cd $DIR/binders/${SEQ}_binder/HMR
#cd prod_md_0p9_cutoff_3dt_64x1
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -s prod_md_500ns.tpr -cpi prod_md_500ns.cpt -append -ntomp $SLURM_CPUS_PER_TASK #-rdd $RDD -npme $PME
