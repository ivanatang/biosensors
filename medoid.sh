#!/bin/bash      

#SBATCH --job-name=sbmedoid
#SBATCH --output=output_medoid_%j.out                  # Output file
#SBATCH --error=error_medoid_%j.err                    # Error file
#SBATCH --account=ucb351_asc4
#SBATCH --partition=amilan
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=ib
#SBATCH --qos=normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module purge
module load anaconda

conda activate mdanalysis

PYTHON=${PYTHON:-python}

# ----------------------------
# Inputs
# ----------------------------
ID=$1

# binder
#TPR="binders/pair_${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_500ns.tpr"
#XTC="binders/pair_${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc"
#TPR="binders/seq${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_500ns.tpr" # seq_ name
#XTC="binders/seq${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc" # seq_ name
#TPR="binders/seq${ID}_binder/HMR/dodecahedron/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_500ns.tpr" # seq_ name
#XTC="binders/seq${ID}_binder/HMR/dodecahedron/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc" # seq_ name

# nonbinder
#TPR="nonbinders/pair_${ID}_nb/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_500ns.tpr"
#XTC="nonbinders/pair_${ID}_nb/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc"
#TPR="nonbinders/seq${ID}_nb/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_500ns.tpr" # seq_ name
#XTC="nonbinders/seq${ID}_nb/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc" # seq_ name
TPR="nonbinders/seq${ID}_nb/HMR/dodecahedron/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_500ns.tpr" # seq_ name
XTC="nonbinders/seq${ID}_nb/HMR/dodecahedron/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc" # seq_ name

# neg low pkt
#TPR="neg_low_pkt/pair_${ID}_low_pkt/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_500ns.tpr"
#XTC="neg_low_pkt/pair_${ID}_low_pkt/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc"

# neg fail gate
#TPR="neg_fail_gate/pair_${ID}_fail_gate/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_500ns.tpr"
#XTC="neg_fail_gate/pair_${ID}_fail_gate/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/prod_md_40_500ns.xtc"

SCRIPT="find_medoid_250ns.py"

# ----------------------------
# Medoid parameters
# ----------------------------
SELECTION="protein and name CA"
STRIDE=1

# binder
#OUTFILE="binders/pair_${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_CA_250ns.txt" # pair_ name
#OUTFILE="binders/seq${ID}_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_CA_250ns.txt" # seq_ name
#OUTFILE="binders/seq${ID}_binder/HMR/dodecahedron/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_CA_250ns.txt" # seq_ name

# nonbinder
#OUTFILE="nonbinders/pair_${ID}_nb/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_CA_250ns.txt" # pair_ name
#OUTFILE="nonbinders/seq${ID}_nb/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_CA_250ns.txt" # seq_ name
OUTFILE="nonbinders/seq${ID}_nb/HMR/dodecahedron/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_CA_250ns.txt" # seq_ name

# neg low pkt
#OUTFILE="neg_low_pkt/pair_${ID}_low_pkt/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_CA_250ns.txt"

# neg fail gate
#OUTFILE="neg_fail_gate/pair_${ID}_fail_gate/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_CA_250ns.txt"
# ----------------------------
# Run medoid calculation
# ----------------------------
echo "Running medoid calculation"
echo "TPR:      $TPR"
echo "XTC:      $XTC"
echo "Select:   $SELECTION"
echo "Stride:   $STRIDE"

$PYTHON $SCRIPT -s "$TPR" -f "$XTC" --select  "$SELECTION" --stride  "$STRIDE" --start-ns 40 --end-ns 250 --out "$OUTFILE"

echo "Medoid calculation complete"
echo "Output written to $OUTFILE"
