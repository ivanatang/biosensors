#!/bin/bash

#SBATCH --job-name=b_prod
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
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

conda activate biosensors

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
#cd $DIR/binders/${SEQ}_binder
#mkdir prod_md_0p9_cutoff_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
#cd prod_md_0p9_cutoff_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
#gmx_mpi grompp -f $MDP/prod_md.mdp -c $DIR/binders/${SEQ}_binder/NPT/npt.gro -t $DIR/binders/${SEQ}_binder/NPT/npt.cpt -p $DIR/binders/${SEQ}_binder/${SEQ}_b.top -o prod_md_500ns.tpr
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK

# Production simulation - HMR (dodecahedron)
# pair_ file name
cd $DIR/binders/bind_${ID}_binder
mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/binders/bind_${ID}_binder/NPT/npt.gro -t $DIR/binders/bind_${ID}_binder/NPT/npt.cpt -p $DIR/binders/bind_${ID}_binder/bind_${ID}_dodecahedron_HMR.top -o prod_md_500ns.tpr
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3

# seq_ file name
#cd $DIR/binders/${ID}_binder
#mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
#cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
#gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/binders/${SEQ}_binder/NPT/npt.gro -t $DIR/binders/${SEQ}_binder/NPT/npt.cpt -p $DIR/binders/${SEQ}_binder/${SEQ}_b_dodecahedron_HMR.top -o prod_md_500ns.tpr
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3

# Production simulation - HMR (cube)
#cd $DIR/binders/${SEQ}_binder/HMR
#mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
#cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
#gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/binders/${SEQ}_binder/HMR/NPT/npt.gro -t $DIR/binders/${SEQ}_binder/HMR/NPT/npt.cpt -p $DIR/binders/${SEQ}_binder/HMR/${SEQ}_b_HMR.top -o prod_md_500ns.tpr
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK #-rdd $RDD -npme $PME
