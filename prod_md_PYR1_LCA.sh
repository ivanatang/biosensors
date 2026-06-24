#!/bin/bash

#SBATCH --job-name=3080b_prod
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
DIR=/projects/ivta1597/biosensors
MDP=$DIR/MDP

# Get sequence value from command line
ID=$1
SEQ_TYPE=$2 # binders | nonbinders | neg_fail_gate | neg_low_pkt
PREFIX=$3 # pair, bind

if [ "$SEQ_TYPE" == "binders" ]; then
    SUFFIX="binder"
elif [ "$SEQ_TYPE" == "nonbinders" ]; then
    SUFFIX="nb"
elif [ "$SEQ_TYPE" == "neg_fail_gate" ]; then
    SUFFIX="fail_gate"
elif [ "$SEQ_TYPE" == "neg_low_pkt" ]; then
    SUFFIX="low_pkt"
else
    echo "ERROR: Unknown SEQ_TYPE '$SEQ_TYPE'" >&2
    exit 1
fi

PME=16
RDD=1.2

D1=6
D2=4
D3=2

# Production simulation - HMR (dodecahedron)
cd $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}
mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}_${PME}PME_${D1}${D2}${D3}dd
gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/NPT/npt.gro -t $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/NPT/npt.cpt -p $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/${PREFIX}_${ID}_dodecahedron_HMR.top -o prod_md_500ns.tpr
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK -npme $PME -dd $D1 $D2 $D3

# Production simulation - HMR (cube)
#cd $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/HMR
#mkdir prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
#cd prod_md_0p9_cutoff_3dt_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
#gmx_mpi grompp -f $MDP/prod_md_HMR_3dt.mdp -c $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/HMR/NPT/npt.gro -t $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/HMR/NPT/npt.cpt -p $DIR/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}/HMR/${SEQ}_b_HMR.top -o prod_md_500ns.tpr
#mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK #-rdd $RDD -npme $PME
