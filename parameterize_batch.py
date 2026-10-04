#!/usr/bin/env python3
"""Parameterizes bile acid (BA) arm complexes into HMR-ready GROMACS inputs.

Generalizes reparameterize_qfix_batch.py (kept unchanged as the record of the
LCA qfix backfill) to multiple ligands and to starting from a Boltz-2 complex
PDB, so the full prep runs as one script:

    1. Split the complex (chain A protein, resname LIG ligand + CONECT) into
       protein_<prefix><id>.pdb and ligand_<prefix><id>.pdb.
    2. PDBFixer: add missing atoms and hydrogens at pH 7.0 (the
       protein_ligand_prep_pipeline.ipynb Step 1 settings).
    3. Ligand bond-order/charge fix from a per-ligand template (see
       LIGAND_TEMPLATES), with formal-charge and stereochemistry checks.
    4. Parameterize/solvate/export exactly as the qfix pipeline does (Sage
       2.0.0 ligand with AM1-BCC, ff14SB protein, TIP3P, counterions only,
       rhombic dodecahedron with 2 nm buffer, packmol, HMR export).

Outputs keep the qfix naming, so em/equil/prod_*_qfix.sh run unchanged on
them (they look for <seq_id>_dodecahedron_HMR_qfix.{top,gro}):
    <out_base>/<subdir>/<seq_id>/
        protein_<prefix><id>.pdb, protein_<prefix><id>_fixed_H.pdb
        ligand_<prefix><id>.pdb, ligand_qfix.sdf
        packmol_solv_dodecahedron_qfix/
        packed_dodecahedron_<seq_id>_qfix.pdb
        <seq_id>_dodecahedron_HMR_qfix.{top,gro,*.itp,*.mdp}

Needs the openff-qfix conda env (AmberTools sqm for AM1-BCC; not
installable in the Python 3.13 biosensors env). Activate it rather than
calling its python by path: OpenFF finds sqm/packmol on PATH, and without
the env's bin/ there AM1-BCC fails with "No registered toolkits can provide
the capability assign_partial_charges".

Usage:
    conda activate openff-qfix
    python parameterize_batch.py                       # all arms in seq_ids_ba.txt
    python parameterize_batch.py --only lcam_0004_binder lca3s_0004_nb
"""
import argparse
import csv
import os

import numpy as np
from openff.interchange import Interchange
from openff.interchange.components._packmol import RHOMBIC_DODECAHEDRON, pack_box
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = os.path.dirname(os.path.abspath(__file__))
ONEDRIVE_BA_BASE = "/Users/ivanatang/Library/CloudStorage/OneDrive-UCB-O365/Shirts Lab/BA_boltz_models"
BOX_SHAPE = "dodecahedron"
SUFFIX_SUBDIR = {"binder": "binders", "nb": "nonbinders"}

# Per-ligand reference chemistry for the bond-order fix. A PDB stores only
# coordinates and raw CONECT connectivity, so RDKit's PDB parser makes every
# bond single and fills valence with implicit H (for LCA's C24 acid that is a
# neutral geminal diol; see reparameterize_qfix_batch.py).
# AssignBondOrdersFromTemplate copies bond orders and formal charges from
# these SMILES onto the PDB coordinates without moving them.
#   lcam:  PubChem CID 9903 (lithocholic acid, 3alpha-OH), carboxylate
#          deprotonated (pKa ~5), net -1. "LCAM" is the collaborator's name
#          for this monoanion; same 27 heavy atoms as LCA.
#   lca3s: the same SMILES with the C3 O replaced by a sulfate monoester
#          (pKa < 1, fully ionized), net -2. Only the substituent on O
#          changes, so C3's chiral tag, written relative to unchanged
#          neighbor order, still encodes 3alpha.
LIGAND_TEMPLATES = {
    "lcam": ("C[C@H](CCC(=O)[O-])[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)C", -1),
    "lca3s": ("C[C@H](CCC(=O)[O-])[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)OS(=O)(=O)[O-])C)C", -2),
}


def parse_args() -> argparse.Namespace:
    """Parses CLI arguments."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seq-list", default=os.path.join(REPO, "seq_ids_ba.txt"))
    p.add_argument("--meta", default=os.path.join(REPO, "ba_arm_metadata.csv"),
                   help="Arm table with seq_id -> source complex PDB (from make_ba_seq_ids.py).")
    p.add_argument("--out-base", default=ONEDRIVE_BA_BASE)
    p.add_argument("--only", nargs="+", help="Restrict to these seq_ids.")
    p.add_argument("--overwrite", action="store_true",
                   help="Re-run arms whose _HMR_qfix.top already exists.")
    return p.parse_args()


def split_seq_id(seq_id: str) -> tuple[str, str, str]:
    """Splits a BA seq_id into its parts.

    Args:
        seq_id: Run name of the form <prefix>_<id>_<suffix>, e.g.
            "lca3s_0004_nb".

    Returns:
        (prefix, id, suffix), e.g. ("lca3s", "0004", "nb").

    Raises:
        ValueError: The prefix has no ligand template or the suffix is unknown.
    """
    prefix, short_id, suffix = seq_id.split("_")
    if prefix not in LIGAND_TEMPLATES or suffix not in SUFFIX_SUBDIR:
        raise ValueError(f"Unrecognized BA seq_id: {seq_id}")
    return prefix, short_id, suffix


def split_complex(complex_pdb: str, protein_pdb: str, ligand_pdb: str) -> None:
    """Writes the protein (chain A) and ligand (resname LIG) of a complex PDB.

    CONECT records are kept for the ligand only; RDKit needs them to build
    the ligand graph, while PDBFixer/OpenFF rebuild protein bonds from
    residue templates.

    Args:
        complex_pdb: Boltz-2 complex PDB.
        protein_pdb: Output path for the protein.
        ligand_pdb: Output path for the ligand.

    Raises:
        ValueError: The complex has no LIG atoms or no ligand CONECT records.
    """
    with open(complex_pdb) as f:
        lines = f.readlines()
    prot = [l for l in lines if l.startswith("ATOM")]
    lig = [l for l in lines if l.startswith("HETATM") and l[17:20] == "LIG"]
    lig_serials = {l[6:11].strip() for l in lig}
    conect = [l for l in lines if l.startswith("CONECT") and l[6:11].strip() in lig_serials]
    if not lig or not conect:
        raise ValueError(f"No LIG atoms/CONECT records in {complex_pdb}")
    with open(protein_pdb, "w") as f:
        f.writelines(prot + ["TER\n", "END\n"])
    with open(ligand_pdb, "w") as f:
        f.writelines(lig + conect + ["END\n"])


def fix_protein(protein_pdb: str, fixed_pdb: str) -> None:
    """Adds missing atoms and pH 7.0 hydrogens with PDBFixer.

    Args:
        protein_pdb: Heavy-atom protein PDB.
        fixed_pdb: Output path for the repaired, protonated protein.
    """
    fixer = PDBFixer(filename=protein_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    with open(fixed_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


def ligand_to_sdf(ligand_pdb: str, prefix: str, sdf_file: str) -> Chem.Mol:
    """Applies the ligand template, checks charge and stereo, and writes SDF.

    Stereo is perceived from the Boltz 3D coordinates (the template's chiral
    tags are not transferred) and every center must match the template. This
    guards against shipping a 3beta-epimer: the handoff README reports C3
    inversions in 103 of 274 structures of the previous set. md_arm_set.csv
    structures already passed that gate; this re-checks it on the exact
    atoms being parameterized.

    Args:
        ligand_pdb: Heavy-atom ligand PDB with CONECT records.
        prefix: Ligand key in LIGAND_TEMPLATES ("lcam" or "lca3s").
        sdf_file: Output SDF path.

    Returns:
        The bond-order-corrected RDKit molecule (heavy atoms, PDB atom order).

    Raises:
        ValueError: Formal charge or any stereocenter disagrees with the
            template.
    """
    smiles, charge = LIGAND_TEMPLATES[prefix]
    template = Chem.MolFromSmiles(smiles)
    mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    if Chem.GetFormalCharge(mol) != charge:
        raise ValueError(f"Formal charge {Chem.GetFormalCharge(mol)} != template {charge}")

    Chem.AssignStereochemistryFrom3D(mol)
    match = mol.GetSubstructMatch(template)  # template atom i -> mol atom match[i]
    want = dict(Chem.FindMolChiralCenters(template, useLegacyImplementation=False))
    got = dict(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))
    bad = [(mol.GetAtomWithIdx(match[t]).GetPDBResidueInfo().GetName().strip(), cip, got.get(match[t]))
           for t, cip in want.items() if got.get(match[t]) != cip]
    if bad:
        raise ValueError("Stereocenter mismatch vs template (atom, expected, found): " + str(bad))

    w = Chem.SDWriter(sdf_file)
    w.write(mol)
    w.close()
    return mol


def parameterize(seq_dir: str, seq_id: str, fixed_pdb: str, sdf_file: str) -> dict:
    """Builds, solvates, neutralizes, and exports the system to GROMACS.

    Identical physics to reparameterize_qfix_batch.reparameterize_one(); see
    that function for the packmol/Interchange workarounds.

    Args:
        seq_dir: Run directory (outputs are written here).
        seq_id: Run name, used for output file names.
        fixed_pdb: PDBFixer-protonated protein PDB.
        sdf_file: Bond-order-corrected ligand SDF.

    Returns:
        Summary with ligand/protein/total charge, water and counterion counts.

    Raises:
        AssertionError: Summed partial charges differ from the formal charges.
    """
    ligand = Molecule.from_file(sdf_file)
    ligand.name = "LIG"
    sage = ForceField("openff_unconstrained-2.0.0.offxml")
    ligand_intrcg = Interchange.from_smirnoff(force_field=sage, topology=[ligand])

    protein = Topology.from_pdb(fixed_pdb).molecule(0)
    protein.name = "protein"
    protein_intrcg = Interchange.from_smirnoff(
        force_field=ForceField("ff14sb_off_impropers_0.0.3.offxml"),
        topology=protein.to_topology(),
    )
    docked_intrcg = protein_intrcg.combine(ligand_intrcg)

    total_charge = round(sum(docked_intrcg["Electrostatics"].charges.values()), 3)
    assert total_charge == protein.total_charge + ligand.total_charge, (
        f"Total charge mismatch: {total_charge} vs {protein.total_charge + ligand.total_charge}"
    )
    total_charge_e = float(total_charge.m)
    n_counterions = int(round(abs(total_charge_e)))
    counterion = None
    if total_charge_e < 0:
        counterion = Molecule.from_smiles("[Na+]")
        counterion.name = "NA"
    elif total_charge_e > 0:
        counterion = Molecule.from_smiles("[Cl-]")
        counterion.name = "CL"

    water = Molecule.from_smiles("O")
    water.name = "SOL"
    water.generate_conformers(n_conformers=1)

    xyz = protein.conformers[0].to(unit.nanometer).m
    protein_radius_nm = np.sqrt(((xyz - xyz.mean(axis=0)) ** 2).sum(axis=1).max())
    scale_nm = 2.0 * protein_radius_nm + 2.0
    box_vectors = (scale_nm * RHOMBIC_DODECAHEDRON) * unit.nanometer
    V_box_nm3 = abs(np.linalg.det(box_vectors.to(unit.nanometer).m))
    V_solute_nm3 = (4.0 / 3.0) * np.pi * protein_radius_nm ** 3
    n_water = int(33.4 * max(V_box_nm3 - V_solute_nm3, 0.0))

    molecules, copies = [water], [n_water]
    if counterion is not None:
        molecules.append(counterion)
        copies.append(n_counterions)
    packed_topology = pack_box(
        solute=Topology.from_molecules([protein, ligand]),
        molecules=molecules,
        number_of_copies=copies,
        box_vectors=box_vectors,
        center_solute=True,
        tolerance=2.0 * unit.angstrom,
        working_directory=os.path.join(seq_dir, f"packmol_solv_{BOX_SHAPE}_qfix"),
        retain_working_files=True,
    )
    packed_topology.to_file(os.path.join(seq_dir, f"packed_{BOX_SHAPE}_{seq_id}_qfix.pdb"))

    solvent = [water] * n_water + ([counterion] * n_counterions if counterion else [])
    system_intrcg = docked_intrcg.combine(Interchange.from_smirnoff(force_field=sage, topology=solvent))
    system_intrcg.positions = packed_topology.get_positions()
    system_intrcg.box = packed_topology.box_vectors

    cwd = os.getcwd()
    os.chdir(seq_dir)
    try:
        system_intrcg.to_gromacs(prefix=f"{seq_id}_{BOX_SHAPE}_HMR_qfix", decimal=3,
                                 hydrogen_mass=3, monolithic=False)
    finally:
        os.chdir(cwd)

    return {
        "ligand_charge": ligand.total_charge.m,
        "protein_charge": protein.total_charge.m,
        "system_charge": total_charge_e,
        "counterion": counterion.name if counterion else "",
        "n_counterions": n_counterions,
        "n_water": n_water,
    }


def process_arm(seq_id: str, complex_pdb: str, out_base: str, overwrite: bool) -> dict | None:
    """Runs the full prep for one arm.

    Args:
        seq_id: BA run name, e.g. "lca3s_0004_nb".
        complex_pdb: Source Boltz-2 complex PDB.
        out_base: Root containing binders/ and nonbinders/.
        overwrite: Re-run even if the HMR topology already exists.

    Returns:
        The parameterize() summary, or None if skipped as already done.
    """
    prefix, short_id, suffix = split_seq_id(seq_id)
    seq_dir = os.path.join(out_base, SUFFIX_SUBDIR[suffix], seq_id)
    if not overwrite and os.path.isfile(os.path.join(seq_dir, f"{seq_id}_{BOX_SHAPE}_HMR_qfix.top")):
        print(f"SKIP {seq_id}: already parameterized")
        return None
    os.makedirs(seq_dir, exist_ok=True)
    print(f"\n=== {seq_id} ===")

    protein_pdb = os.path.join(seq_dir, f"protein_{prefix}{short_id}.pdb")
    fixed_pdb = os.path.join(seq_dir, f"protein_{prefix}{short_id}_fixed_H.pdb")
    ligand_pdb = os.path.join(seq_dir, f"ligand_{prefix}{short_id}.pdb")
    sdf_file = os.path.join(seq_dir, "ligand_qfix.sdf")

    split_complex(complex_pdb, protein_pdb, ligand_pdb)
    fix_protein(protein_pdb, fixed_pdb)
    mol = ligand_to_sdf(ligand_pdb, prefix, sdf_file)
    print(f"  ligand: {mol.GetNumAtoms()} heavy atoms, charge {Chem.GetFormalCharge(mol):+d}, stereo OK")
    summary = parameterize(seq_dir, seq_id, fixed_pdb, sdf_file)
    print(f"  protein {summary['protein_charge']:+.0f}, system {summary['system_charge']:+.0f} -> "
          f"{summary['n_counterions']} {summary['counterion']}, {summary['n_water']} waters")
    return summary


def main() -> None:
    """Parameterizes every requested arm and writes a per-run summary CSV."""
    args = parse_args()
    with open(args.meta, newline="") as f:
        pdb_of = {r["seq_id"]: os.path.join(REPO, r["pdb"]) for r in csv.DictReader(f)}
    with open(args.seq_list) as f:
        seq_ids = [l.split("\t")[0].strip() for l in f if l.strip() and not l.startswith("#")]
    if args.only:
        unknown = set(args.only) - set(seq_ids)
        if unknown:
            raise SystemExit(f"Not in {args.seq_list}: {sorted(unknown)}")
        seq_ids = [s for s in seq_ids if s in args.only]

    results, failed = [], []
    for seq_id in seq_ids:
        try:
            summary = process_arm(seq_id, pdb_of[seq_id], args.out_base, args.overwrite)
            if summary:
                results.append({"seq_id": seq_id, **summary})
        except Exception as e:
            print(f"  FAILED {seq_id}: {type(e).__name__}: {e}")
            failed.append(seq_id)

    if results:
        log = os.path.join(args.out_base, "parameterize_summary.csv")
        new = not os.path.isfile(log)
        with open(log, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0]))
            if new:
                w.writeheader()
            w.writerows(results)
        print(f"\nSummary appended to {log}")
    print(f"=== Done: {len(results)} parameterized, {len(failed)} failed ===")
    if failed:
        print(f"Failed: {failed}")


if __name__ == "__main__":
    main()
