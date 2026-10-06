#!/bin/bash

#SBATCH --job-name=em_PYR1_BA
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# em_PYR1_BA.sh -- two-stage EM for the bile acid (LCA / LCA-3-S) arm set.
# Same inputs/outputs as em_PYR1_LCA_qfix.sh (reads the _qfix topology,
# writes EM_qfix/em.gro for equil_PYR1_LCA_qfix.sh), with a restrained stage
# first:
#   1. EM_qfix/em_restr.*  heavy atoms of protein + ligand restrained
#                          (MDP/em_restrained.mdp, -DEM_RESTRAIN)
#   2. EM_qfix/em.*        unrestrained (MDP/em.mdp), from stage 1's output
#
# Why: 59/228 Boltz-2 arms carry non-polar protein-ligand contacts < 2.8 A
# (needs_restrained_min in md_arm_set.csv). Run on every arm, not just
# those, so the protocol is identical across classes and ligands; the flag is
# more common in the LCA3S arm, so a flagged-only step would correlate
# protocol with ligand.
#
# Usage:
#   sbatch em_PYR1_BA.sh <ID> <SEQ_TYPE> <PREFIX>      e.g. 0004 nonbinders lca3s

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate biosensors

# The /projects GROMACS 2025.3 build (first on PATH) needs the env's newer
# libstdc++, not the HPC image's old /lib64/libstdc++.so.6 -- without this
# gmx fails to start: "/lib64/libstdc++.so.6: version GLIBCXX_3.4.29 not found"
export LD_LIBRARY_PATH="/projects/ivta1597/software/anaconda/envs/biosensors/lib:$LD_LIBRARY_PATH"

DIR=/projects/ivta1597/biosensors
MDP=$DIR/MDP
BASE=${BASE:-/scratch/alpine/ivta1597/BA_boltz_models}
RESTRAINT_FC=1000   # kJ/mol/nm^2

ID=$1
SEQ_TYPE=$2 # binders | nonbinders
PREFIX=$3   # lcam | lca3s

case "$SEQ_TYPE" in
    binders)    SUFFIX="binder" ;;
    nonbinders) SUFFIX="nb" ;;
    *) echo "ERROR: Unknown SEQ_TYPE '$SEQ_TYPE'" >&2; exit 1 ;;
esac

# #SBATCH directives are static text parsed before ID/SEQ_TYPE/PREFIX are
# known, so relabel the job now that they are.
scontrol update JobId=$SLURM_JOB_ID JobName="${PREFIX}_${ID}_em_qfix" 2>/dev/null

NAME=${PREFIX}_${ID}_${SUFFIX}
SEQ_DIR=$BASE/${SEQ_TYPE}/${NAME}
TOP=$SEQ_DIR/${NAME}_dodecahedron_HMR_qfix.top
GRO=$SEQ_DIR/${NAME}_dodecahedron_HMR_qfix.gro
LIG_ITP=$SEQ_DIR/${NAME}_dodecahedron_HMR_qfix_LIG.itp
for f in "$TOP" "$GRO" "$LIG_ITP"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: input not found: $f" >&2
        exit 1
    fi
done

cd "$SEQ_DIR"

# Heavy-atom restraint files. Atom numbers in a position-restraint itp are
# relative to their own moleculetype. Protein is the first molecule in the
# .gro, so genrestr's system numbering on group 2 (Protein-H) equals its
# molecule numbering. The ligand's must come from its own .itp: OpenFF names
# atoms by element, so heavy atoms are those whose name doesn't start with H
# (mass can't be used -- HMR moves mass onto hydrogens).
if [ ! -f restr_em_protein.itp ]; then
    echo 2 | gmx genrestr -f "$GRO" -o restr_em_protein.itp -fc $RESTRAINT_FC $RESTRAINT_FC $RESTRAINT_FC
fi
if [ ! -f restr_em_lig.itp ]; then
    awk -v fc=$RESTRAINT_FC '
        /^\[ *atoms *\]/ { inatoms = 1; print "[ position_restraints ]"; next }
        /^\[/            { inatoms = 0 }
        inatoms && $1 ~ /^[0-9]+$/ && $5 !~ /^H/ { printf "%6d     1  %d  %d  %d\n", $1, fc, fc, fc }
    ' "$LIG_ITP" > restr_em_lig.itp
fi

# Insert each restraint include right after the include of the moleculetype
# it belongs to (the current moleculetype runs until the next [ moleculetype ]).
if ! grep -q "EM_RESTRAIN" "$TOP"; then
    awk '
        { print }
        /^#include.*_protein\.itp"/ { print "#ifdef EM_RESTRAIN"; print "#include \"restr_em_protein.itp\""; print "#endif" }
        /^#include.*_LIG\.itp"/     { print "#ifdef EM_RESTRAIN"; print "#include \"restr_em_lig.itp\"";     print "#endif" }
    ' "$TOP" > "${TOP}.tmp" && mv "${TOP}.tmp" "$TOP"
    [ "$(grep -c EM_RESTRAIN "$TOP")" -eq 2 ] || { echo "ERROR: restraint includes not inserted in $TOP" >&2; exit 1; }
fi

mkdir -p EM_qfix
cd EM_qfix
gmx grompp -f $MDP/em_restrained.mdp -c "$GRO" -r "$GRO" -p "$TOP" -o em_restr.tpr
gmx mdrun -deffnm em_restr
[ -f em_restr.gro ] || { echo "ERROR: restrained EM produced no em_restr.gro" >&2; exit 1; }

gmx grompp -f $MDP/em.mdp -c em_restr.gro -p "$TOP" -o em.tpr
gmx mdrun -deffnm em
[ -f em.gro ] || { echo "ERROR: unrestrained EM produced no em.gro" >&2; exit 1; }
