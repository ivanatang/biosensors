#!/bin/bash

#SBATCH --job-name=eq_PYR1_LCA_qfix
#SBATCH --output=output_%j.out                  # Output file
#SBATCH --error=error_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# equil_PYR1_LCA_qfix.sh -- NVT/NPT equilibration for the _qfix pilot
# systems. Reads EM_qfix/em.gro (from em_PYR1_LCA_qfix.sh) and the _qfix
# topology, writes to NVT_qfix/ and NPT_qfix/.
#
# Usage:
#   sbatch equil_PYR1_LCA_qfix.sh <ID> <SEQ_TYPE> <PREFIX>

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
BASE=/scratch/alpine/ivta1597/LCA_boltz_models

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

# #SBATCH directives are static text parsed before ID/SEQ_TYPE/PREFIX are
# known, so relabel the job now that they are, rather than leaving every
# submission of this script showing the same generic name in squeue.
scontrol update JobId=$SLURM_JOB_ID JobName="${PREFIX}_${ID}_eq_qfix" 2>/dev/null

# Navigate to the sequence's scratch directory for trajectory files
SEQ_DIR=$BASE/${SEQ_TYPE}/${PREFIX}_${ID}_${SUFFIX}
if [ ! -d "$SEQ_DIR" ]; then
    echo "ERROR: Sequence directory not found: $SEQ_DIR" >&2
    exit 1
fi

TOP=$SEQ_DIR/${PREFIX}_${ID}_${SUFFIX}_dodecahedron_HMR_qfix.top
if [ ! -f "$SEQ_DIR/EM_qfix/em.gro" ]; then
    echo "ERROR: EM_qfix/em.gro not found -- run em_PYR1_LCA_qfix.sh first" >&2
    exit 1
fi

# Position restraints: matches run_genrstr.sh's original protocol (backbone
# restraints, 1000 kJ/mol/nm^2), which patches the standard-pipeline .top
# files but was never run against the _qfix topology. nvt.mdp/npt.mdp both
# reference -DPOSRES, so grompp fails ("macro defined but not used") without
# this. The protein topology is identical between the original and _qfix
# systems (only the ligand changed), but this regenerates from EM_qfix/em.gro
# directly rather than assuming the original posre_protein.itp is reusable.
POSRE=$SEQ_DIR/posre_protein_qfix.itp
if [ ! -f "$POSRE" ]; then
    ( cd "$SEQ_DIR" && echo 4 | gmx genrestr -f EM_qfix/em.gro -o posre_protein_qfix.itp -fc 1000 1000 1000 )
fi
if ! grep -q "POSRES" "$TOP"; then
    LINE=$(grep -in '^#include.*protein' "$TOP" | head -n1 | cut -d: -f1)
    if [ -z "$LINE" ]; then
        echo "ERROR: no '#include ... protein' line found in $TOP" >&2
        exit 1
    fi
    awk -v n="$LINE" '
        { print }
        NR == n {
            print "#ifdef POSRES"
            print "#include \"posre_protein_qfix.itp\""
            print "#endif"
        }
    ' "$TOP" > "${TOP}.tmp" && mv "${TOP}.tmp" "$TOP"
    echo "Inserted POSRES include after line $LINE: $TOP"
fi

# dodecahedron unit cell
# NVT
cd "$SEQ_DIR"
mkdir -p NVT_qfix
cd NVT_qfix
gmx grompp -f $MDP/nvt.mdp -c $SEQ_DIR/EM_qfix/em.gro -r $SEQ_DIR/EM_qfix/em.gro -p "$TOP" -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
cd "$SEQ_DIR"
mkdir -p NPT_qfix
cd NPT_qfix
gmx grompp -f $MDP/npt.mdp -c $SEQ_DIR/NVT_qfix/nvt.gro -t $SEQ_DIR/NVT_qfix/nvt.cpt -p "$TOP" -r $SEQ_DIR/NVT_qfix/nvt.gro -o npt.tpr
gmx mdrun -deffnm npt
