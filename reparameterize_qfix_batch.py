#!/usr/bin/env python3
"""Backfills the _qfix topology on sequences missed by earlier qfix batches.

Standalone (non-notebook) reparameterization script. Mirrors
protein_ligand_prep_pipeline.ipynb's Step 2 (bond-order fix) and Step 3
(parameterize/solvate/export) exactly, driven by a plain ID list instead of
manual notebook cell edits.

Reuses each sequence's existing protein_<prefix><id>_fixed_H.pdb (protein
topology is unaffected by the ligand bond-order fix, so Step 1/PDBFixer is
skipped), reads the existing ligand_<prefix><id>.pdb, applies the same
AssignBondOrdersFromTemplate fix, and writes GROMACS output with the
project's established _qfix naming convention (verified against the
already-completed bind_022_binder pilot sequence):
    ligand_qfix.sdf
    packmol_solv_<box_shape>_qfix/
    packed_<box_shape>_<prefix>_<id>_<suffix>_qfix.pdb
    <prefix>_<id>_<suffix>_<box_shape>_HMR_qfix.{top,gro,*.itp,*.mdp}

Writes to the same OneDrive-mounted path the notebook does, not the git
repo, since that's the shared location production Alpine jobs sync from.

Usage:
    python reparameterize_qfix_batch.py --ids 3069 3070 --seq-type binder
"""
import argparse
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit
from openff.interchange import Interchange
from openff.interchange.components._packmol import RHOMBIC_DODECAHEDRON, pack_box

ONEDRIVE_BASE = "/Users/ivanatang/Library/CloudStorage/OneDrive-UCB-O365/Shirts Lab/LCA_boltz_models"
PREFIX = "pair"
BOX_SHAPE = "dodecahedron"

SEQ_TYPE_INFO = {
    "binder":       ("binders",      "binder"),
    "nonbinder":    ("nonbinders",   "nb"),
    "neg_low_pkt":  ("neg_low_pkt",  "low_pkt"),
    "neg_fail_gate":("neg_fail_gate","fail_gate"),
}

# Reference chemistry for the bond-order fix, identical to
# protein_ligand_prep_pipeline.ipynb's LCA_TEMPLATE. A PDB file stores only
# atom positions and raw CONECT connectivity, never bond order or formal
# charge, so RDKit's PDB parser defaults every bond to single and fills
# remaining valence with implicit H. For LCA's C24 carboxylic acid this
# silently produces a neutral geminal diol (-CH(OH)(OH)) instead of the
# deprotonated carboxylate (LCA-, pKa ~5, so deprotonated at the modeled pH
# 7.4 -- the same ionization state fed into the Boltz-2 structure
# prediction as SMILES). AssignBondOrdersFromTemplate maps the correct bond
# orders/charge from this reference SMILES onto the PDB-derived 3D
# coordinates, leaving the coordinates themselves untouched.
# Source: PubChem CID 9903 (lithocholic acid), deprotonated at the acid.
LCA_TEMPLATE = Chem.MolFromSmiles(
    "C[C@H](CCC(=O)[O-])[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)C"
)


def parse_args():
    """Parses CLI arguments for the qfix backfill batch."""
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ids", nargs="+", required=True, help="Sequence ID numbers, e.g. 3069 3070")
    p.add_argument("--seq-type", required=True, choices=list(SEQ_TYPE_INFO))
    return p.parse_args()


def reparameterize_one(seq_id):
    """Applies the qfix bond-order fix and re-runs parameterize/solvate/export.

    Uses the module-level `args.seq_type` to resolve the sequence's
    directory (see parse_args).

    Args:
        seq_id (str): Sequence ID number (e.g. "3069").

    Raises:
        FileNotFoundError: The sequence directory or a required input file
            (fixed protein PDB, ligand PDB) doesn't exist.
    """
    folder_name, suffix = SEQ_TYPE_INFO[args.seq_type]
    seq_dir = os.path.join(ONEDRIVE_BASE, folder_name, f"{PREFIX}_{seq_id}_{suffix}")
    if not os.path.isdir(seq_dir):
        raise FileNotFoundError(f"Sequence directory not found: {seq_dir}")

    pdb_file = os.path.join(seq_dir, f"protein_{PREFIX}{seq_id}_fixed_H.pdb")
    ligand_pdb = os.path.join(seq_dir, f"ligand_{PREFIX}{seq_id}.pdb")
    sdf_file = os.path.join(seq_dir, "ligand_qfix.sdf")
    for f in (pdb_file, ligand_pdb):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Required input not found: {f}")

    print(f"\n=== Processing {PREFIX}_{seq_id}_{suffix} ===")

    # --- Step 2: bond-order fix ---
    mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
    mol = AllChem.AssignBondOrdersFromTemplate(LCA_TEMPLATE, mol)
    w = Chem.SDWriter(sdf_file)
    w.write(mol)
    w.close()
    print(f"  Bond-order fix applied (charge={Chem.GetFormalCharge(mol):+d}) -> {sdf_file}")

    # --- Step 3: parameterize, solvate, export ---
    ligand = Molecule.from_file(sdf_file)
    ligand.name = "LIG"
    ligand_intrcg = Interchange.from_smirnoff(
        force_field=ForceField("openff_unconstrained-2.0.0.offxml"),
        topology=[ligand],
    )

    protein_full = Topology.from_pdb(pdb_file)
    protein = protein_full.molecule(0)
    protein.name = "protein"
    ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")
    protein_intrcg = Interchange.from_smirnoff(
        force_field=ff14sb,
        topology=protein.to_topology(),
    )

    docked_intrcg = protein_intrcg.combine(ligand_intrcg)

    total_charge = round(sum(docked_intrcg["Electrostatics"].charges.values()), 3)
    assert total_charge == protein.total_charge + ligand.total_charge, (
        f"Total charge mismatch: {total_charge} vs {protein.total_charge + ligand.total_charge}"
    )
    print(f"  Total charge: {total_charge}")

    total_charge_e = float(total_charge.m)
    if total_charge_e < 0:
        counterion = Molecule.from_smiles("[Na+]")
        counterion.name = "NA"
        n_counterions = int(round(abs(total_charge_e)))
    elif total_charge_e > 0:
        counterion = Molecule.from_smiles("[Cl-]")
        counterion.name = "CL"
        n_counterions = int(round(abs(total_charge_e)))
    else:
        counterion = None
        n_counterions = 0

    water = Molecule.from_smiles("O")
    water.name = "SOL"
    water.generate_conformers(n_conformers=1)

    xyz = protein.conformers[0].to(unit.nanometer).m
    centroid = xyz.mean(axis=0)
    protein_radius_nm = np.sqrt(((xyz - centroid) ** 2).sum(axis=1).max())
    buffer_nm = 2.0
    scale_nm = 2.0 * protein_radius_nm + buffer_nm

    box_vectors = (scale_nm * RHOMBIC_DODECAHEDRON) * unit.nanometer
    box_nm = box_vectors.to(unit.nanometer).m
    V_box_nm3 = abs(np.linalg.det(box_nm))
    V_solute_nm3 = (4.0 / 3.0) * np.pi * (protein_radius_nm ** 3)
    n_water = int(33.4 * max(V_box_nm3 - V_solute_nm3, 0.0))
    packmol_dir = os.path.join(seq_dir, f"packmol_solv_{BOX_SHAPE}_qfix")
    print(f"  Protein radius: {protein_radius_nm:.3f} nm, scale: {scale_nm:.3f} nm, box vol: {V_box_nm3:.2f} nm^3")
    print(f"  Waters: {n_water}, counterions: {n_counterions}")

    molecules = [water]
    number_of_copies = [n_water]
    if counterion is not None and n_counterions > 0:
        molecules.append(counterion)
        number_of_copies.append(n_counterions)

    # docked_intrcg.topology is an InterchangeTopology, which this
    # openff-interchange version deliberately stores no positions on (use
    # Interchange.positions instead, per its own NoPositionsError message).
    # pack_box's solute needs an actual positioned Topology, so build one
    # directly from the two Molecules, each already carrying the correct
    # docked conformer from its own file.
    docked_topology = Topology.from_molecules([protein, ligand])
    packed_topology = pack_box(
        solute=docked_topology,
        molecules=molecules,
        number_of_copies=number_of_copies,
        box_vectors=box_vectors,
        center_solute=True,
        tolerance=2.0 * unit.angstrom,
        working_directory=packmol_dir,
        retain_working_files=True,
    )
    print(f"  Total molecules packed: {packed_topology.n_molecules}")

    packed_pdb = os.path.join(seq_dir, f"packed_{BOX_SHAPE}_{PREFIX}_{seq_id}_{suffix}_qfix.pdb")
    packed_topology.to_file(packed_pdb)

    topology_molecules = [water] * n_water
    if counterion is not None:
        topology_molecules += [counterion] * n_counterions

    water_intrcg = Interchange.from_smirnoff(
        force_field=ForceField("openff_unconstrained-2.0.0.offxml"),
        topology=topology_molecules,
    )
    system_intrcg = docked_intrcg.combine(water_intrcg)
    system_intrcg.positions = packed_topology.get_positions()
    system_intrcg.box = packed_topology.box_vectors

    cwd = os.getcwd()
    os.chdir(seq_dir)
    try:
        out_prefix = f"{PREFIX}_{seq_id}_{suffix}_{BOX_SHAPE}_HMR_qfix"
        system_intrcg.to_gromacs(prefix=out_prefix, decimal=3, hydrogen_mass=3, monolithic=False)
        print(f"  GROMACS files written: {seq_dir}/{out_prefix}.{{top,gro}}")
    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    args = parse_args()
    failed = []
    for seq_id in args.ids:
        try:
            reparameterize_one(seq_id)
        except Exception as e:
            print(f"  FAILED {seq_id}: {e}")
            failed.append(seq_id)

    print(f"\n=== Done: {len(args.ids) - len(failed)}/{len(args.ids)} succeeded ===")
    if failed:
        print(f"Failed: {failed}")
