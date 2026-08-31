#!/bin/bash
# audit_qfix_pilot_features.sh
# ─────────────────────────────────────────────────────────────────────────────
# For each of the 6 qfix pilot sequences, checks every feature family needed
# by build_new_ligand_feat_table.py's schema (dw_pocket whole+core+tail,
# contact-type, salt bridge, hydration ligand+pocket, gate-latch water
# bridge, plus the Phase-3 post-processing prerequisite that RMSD-to-ref/RMSF
# depend on) for BOTH the standard run and the _qfix run, side by side --
# so it's clear at a glance how much of the qfix post-processing pipeline
# still needs to run before a feature comparison or model-swap eval is
# possible.
#
# Usage: bash audit_qfix_pilot_features.sh
# ─────────────────────────────────────────────────────────────────────────────

SCRATCH="/scratch/alpine/ivta1597/LCA_boltz_models"
ARCHIVE="/pl/active/shirts_archive/IvanaTang/biosensors"
REPO="/projects/ivta1597/biosensors"
RUNREL="prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"

declare -A DIR_TYPE=(
    [bind_022_binder]=binders  [bind_019_binder]=binders  [bind_020_binder]=binders
    [nonb_006_nb]=nonbinders   [nonb_008_nb]=nonbinders   [nonb_009_nb]=nonbinders
)
SEQ_IDS=(bind_022_binder bind_019_binder bind_020_binder nonb_006_nb nonb_008_nb nonb_009_nb)

check() {
    local label="$1" std_path="$2" qfix_path="$3"
    local std="MISS" qfix="MISS"
    [[ -f "$std_path" ]] && std="OK"
    [[ -f "$qfix_path" ]] && qfix="OK"
    printf "  %-32s std=%-4s qfix=%-4s\n" "$label" "$std" "$qfix"
}

for seq_id in "${SEQ_IDS[@]}"; do
    dir_type="${DIR_TYPE[$seq_id]}"
    echo "============================================================"
    echo "  $seq_id  ($dir_type)"
    echo "============================================================"

    # Phase-3 post-processing prerequisite (feeds RMSD-to-ref/RMSF/many others)
    check "Phase 3 (medoid_PL.pdb)" \
        "${SCRATCH}/${dir_type}/${seq_id}/${RUNREL}/500ns/medoid_PL.pdb" \
        "${SCRATCH}/${dir_type}/${seq_id}/${RUNREL}_qfix/500ns/medoid_PL.pdb"

    check "RMSF (rmsf_PL.xvg)" \
        "${SCRATCH}/${dir_type}/${seq_id}/${RUNREL}/500ns/rmsf_PL.xvg" \
        "${SCRATCH}/${dir_type}/${seq_id}/${RUNREL}_qfix/500ns/rmsf_PL.xvg"

    # extract_gate_latch_rmsd_feats.py reads/writes directly against the
    # PetaLibrary archive (BASE=archive in that script), not scratch --
    # unlike every other per-sequence output checked here.
    check "RMSD-to-ref (gate_rmsd_to_ref.xvg)" \
        "${ARCHIVE}/${dir_type}/${seq_id}/${RUNREL}/500ns/gate_rmsd_to_ref.xvg" \
        "${ARCHIVE}/${dir_type}/${seq_id}/${RUNREL}_qfix/500ns/gate_rmsd_to_ref.xvg"

    check "R-scores (whole)" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_contacts_40_500ns/${seq_id}_R_scores_40_500ns.csv" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_contacts_40_500ns_qfix/${seq_id}_R_scores_40_500ns.csv"

    check "R-scores (core)" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_contacts_40_500ns_core/${seq_id}_R_scores_40_500ns_core.csv" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_contacts_40_500ns_core_qfix/${seq_id}_R_scores_40_500ns_core.csv"

    check "R-scores (tail)" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_contacts_40_500ns_tail/${seq_id}_R_scores_40_500ns_tail.csv" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_contacts_40_500ns_tail_qfix/${seq_id}_R_scores_40_500ns_tail.csv"

    check "Contact-type" \
        "${REPO}/LIG_contacts/contact_type_results_40_500ns/${seq_id}_contact_summary_40_500ns.csv" \
        "${REPO}/LIG_contacts/contact_type_results_40_500ns_qfix/${seq_id}_contact_summary_40_500ns.csv"

    check "Salt bridge" \
        "${ARCHIVE}/${dir_type}/${seq_id}/${RUNREL}/salt_bridge/saltbridge_occupancy_full.csv" \
        "${ARCHIVE}/${dir_type}/${seq_id}/${RUNREL}_qfix/salt_bridge/saltbridge_occupancy_full.csv"

    check "Hydration (ligand)" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_density_40_500ns/${seq_id}_hydration_ligand_40_500ns.csv" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_density_40_500ns_qfix/${seq_id}_hydration_ligand_40_500ns.csv"

    check "Hydration (pocket)" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_density_40_500ns/${seq_id}_hydration_pocket_40_500ns.csv" \
        "${SCRATCH}/${dir_type}/${seq_id}/water_density_40_500ns_qfix/${seq_id}_hydration_pocket_40_500ns.csv"

    check "Gate-latch water bridge" \
        "${SCRATCH}/${dir_type}/${seq_id}/gate_latch_water_bridge_0_500ns_core/${seq_id}_gate_latch_bridge_0_500ns_core.csv" \
        "${SCRATCH}/${dir_type}/${seq_id}/gate_latch_water_bridge_0_500ns_core_qfix/${seq_id}_gate_latch_bridge_0_500ns_core.csv"

    echo ""
done

echo "============================================================"
echo "std = standard (original ligand parameterization, already in feat_table_500ns.xlsx)"
echo "qfix = bond-order-fixed reparameterization, being evaluated as a replacement"
echo "============================================================"
