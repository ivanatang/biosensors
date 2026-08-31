#!/bin/bash
# audit_new_ligand_features.sh
# ─────────────────────────────────────────────────────────────────────────────
# Checks every feature-extraction stage's final output for the 4 new-ligand
# test sequences (cdca/glca/lca3s/lca_001_binder) and prints OK/MISSING per
# stage, so you know exactly what still needs to run -- rather than
# discovering gaps one aggregator error at a time.
#
# Usage: bash audit_new_ligand_features.sh
# ─────────────────────────────────────────────────────────────────────────────

SCRATCH="/scratch/alpine/ivta1597/LCA_boltz_models"
ARCHIVE="/pl/active/shirts_archive/IvanaTang/biosensors"
REPO="/projects/ivta1597/biosensors"
SEQ_IDS=(cdca_001_binder glca_001_binder lca3s_001_binder lca_001_binder)
DIR_TYPE="binders"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

check() {
    local label="$1" path="$2"
    if [[ -f "$path" ]]; then
        echo "  OK      $label"
    else
        echo "  MISSING $label  ->  $path"
    fi
}

for seq_id in "${SEQ_IDS[@]}"; do
    echo "============================================================"
    echo "  $seq_id"
    echo "============================================================"

    check "R-scores (water_analysis)" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/water_contacts_40_500ns/${seq_id}_R_scores_40_500ns.csv"

    check "Hbond threshold" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/water_contacts_40_500ns/${seq_id}_thresholds_40_500ns.json"

    check "Water Hbond stability" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/water_contacts_40_500ns/${seq_id}_hbond_summary_40_500ns.csv"

    check "Gate-latch water bridge" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/gate_latch_water_bridge_0_500ns_core/${seq_id}_gate_latch_bridge_0_500ns_core.csv"

    check "Contact-type analysis" \
        "${REPO}/LIG_contacts/contact_type_results_40_500ns/${seq_id}_contact_summary_40_500ns.csv"

    check "Salt bridge" \
        "${ARCHIVE}/${DIR_TYPE}/${seq_id}/${RUNREL}/salt_bridge/saltbridge_occupancy_full.csv"

    check "Hydration (water_spatial)" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/water_density_40_500ns/${seq_id}_hydration_ligand_40_500ns.csv"

    check "Pocket exploration (mdpocket freq_iso)" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/${RUNREL}/mdpocket_${seq_id}_freq_iso_0_5.pdb"

    check "Pocket characterization (mdpocket descriptors)" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/${RUNREL}/mdpocket_${seq_id}_descriptors.txt"

    check "Rg" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/${RUNREL}/Rg_PL.xvg"

    check "SASA" \
        "${SCRATCH}/${DIR_TYPE}/${seq_id}/${RUNREL}/sasa_PL.xvg"

    echo ""
    echo "  In combined RMSD-to-ref summary:"
    if grep -qE "(^|,)${seq_id}(,|$)" "${REPO}/analysis/gate_latch_rmsd_to_ref_summary_500ns.csv" 2>/dev/null; then
        echo "  OK      RMSD-to-ref (row present)"
    else
        echo "  MISSING RMSD-to-ref (no row in gate_latch_rmsd_to_ref_summary_500ns.csv)"
    fi

    if grep -qE "(^|,)${seq_id}(,|$)" "${REPO}/analysis/rmsf_single_residues_per_seq_500ns.csv" 2>/dev/null; then
        echo "  OK      RMSF (row present)"
    else
        echo "  MISSING RMSF (no row in rmsf_single_residues_per_seq_500ns.csv)"
    fi

    echo ""
done

echo "============================================================"
echo "Done. Re-run this after submitting fixes to confirm."
echo "============================================================"
