#!/bin/bash

#SBATCH --job-name=nb_prod
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

# Production simulation - nonbinder
# pair_ file name
cd $DIR/nonbinders/nonb_${ID}_nb
mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/nonbinders/nonb_${ID}_nb/NPT/npt.gro -t $DIR/nonbinders/nonb_${ID}_nb/NPT/npt.cpt -p $DIR/nonbinders/nonb_${ID}_nb/nonb_${ID}_dodecahedron_HMR.top -o prod_md_500ns.tpr
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3

# seq_ file name
#cd $DIR/nonbinders/${ID}_nb
#mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
#cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
#gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/nonbinders/${ID}_nb/NPT/npt.gro -t $DIR/nonbinders/${ID}_nb/NPT/npt.cpt -p $DIR/nonbinders/${ID}_nb/${ID}_nb_dodecahedron_HMR.top -o prod_md_500ns.tpr
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3

# cube
#cd $DIR/nonbinders/${SEQ}_nb/HMR
#mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
#cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
#gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/nonbinders/${SEQ}_nb/HMR/NPT/npt.gro -t $DIR/nonbinders/${SEQ}_nb/HMR/NPT/npt.cpt -p $DIR/nonbinders/${SEQ}_nb/HMR/${SEQ}_nb_HMR.top -o prod_md_500ns.tpr
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK #-rdd $RDD -npme $PME

#$GMX/gmx grompp -f $MDP/prod_md.mdp -c $DIR/nonbinders/${SEQ}_nb/NPT/npt.gro -t $DIR/nonbinders/${SEQ}_nb/NPT/npt.cpt -p $DIR/nonbinders/${SEQ}_nb/${SEQ}_nb.top -n $DIR/nonbinders/${SEQ}_nb/index.ndx -o prod_md_100ns.tpr
#$GMX/gmx mdrun -deffnm prod_md_100ns -nt 12
