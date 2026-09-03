#!/usr/bin/env python
"""Renders seq10_binder with 5 pocket residues each in a distinct color.

Single publication-quality figure of the seq10_binder medoid structure
(PYR1+LCA biosensor), for a research talk. Single structure, single camera
(not a multi-panel comparison like medoid_comparison.py in this directory).
The 5 residues below are individual positions of interest, not loop ranges
like the gate (resi 84-90) / latch (resi 114-118) in other pymol_renders/
scripts -- don't confuse these with that unrelated figure.

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/five_residue_highlight.py

Produces one ray-traced PNG in pymol_renders/output/:

    seq10_binder_5residue_highlight.png

Each residue is shown as a colored cartoon segment (reads even when the
sidechain points away from camera) plus sticks for the full residue
(backbone + sidechain, not just sidechain) in the same color. The rest of
the protein is neutral gray cartoon; the ligand (chain B, resn LIG) is
always shown in full as yellow-carbon sticks. Background is plain white
(cmd.bg_color("white")) -- this project moved away from tinted backgrounds
in an earlier iteration; keep it plain white.

Camera: one view, oriented and zoomed on the union of all 5 residues' atoms
plus the full ligand (cmd.orient then cmd.zoom with generous buffer), so
nothing is cropped. Deliberately wider than the tight gate/latch close-ups
elsewhere in pymol_renders/, since these 5 residues are spread across
different parts of the pocket (some near the ligand's steroid core, some
near its tail/carboxylate end).

No in-scene labels: this script only renders the colored structure.
Residue labels (e.g. "Phe61") are added afterward by a separate
pixel-space post-processing script (PIL-based), consistent with this
project's decision to abandon in-scene PyMOL labels (see
medoid_comparison.py docstring for why) in favor of post-processing that
can actually be previewed/verified.

To re-run on a different structure/residue set, edit PDB_PATH and
RESIDUES below.
"""

import os

from pymol import cmd, util

# --------------------------------------------------------------------------
# Config - edit these to re-point the script at a different structure/set
# --------------------------------------------------------------------------
REPO_ROOT = "/Users/ivanatang/Developer/biosensors"
OUT_DIR = os.path.join(REPO_ROOT, "pymol_renders", "output")

PDB_PATH = os.path.join(
    REPO_ROOT,
    "binders/seq10_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_PL.pdb",
)
OUT_PNG = os.path.join(OUT_DIR, "seq10_binder_5residue_highlight.png")

OBJ_NAME = "seq10_binder"

PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
LIGAND_RESN = "LIG"
PROTEIN_GRAY = "gray80"

# The 5 highlighted residues, each a distinct colorblind-safe color.
# These are single residue positions (not loop ranges) -- see module
# docstring.
RESIDUES = [
    {"resi": 61, "hex": "#E69F00"},   # orange
    {"resi": 62, "hex": "#56B4E9"},   # sky blue
    {"resi": 81, "hex": "#009E73"},   # bluish green
    {"resi": 116, "hex": "#B2182B"},  # strong crimson/red (was vermillion #D55E00,
                                       # too close in hue to resi 61's orange at
                                       # ray-traced stick scale)
    {"resi": 164, "hex": "#CC79A7"},  # reddish purple
]

IMG_WIDTH = 1800
IMG_HEIGHT = 1400
PNG_DPI = 300

# Camera framing knobs
ZOOM_BUFFER = 6.0  # padding (Angstrom) around the framed selection


# --------------------------------------------------------------------------
# Small color helper (pure python, no numpy dependency)
# --------------------------------------------------------------------------
def hex_to_rgb01(hex_code):
    """Converts a "#RRGGBB" hex color to a 0-1 RGB tuple for PyMOL.

    Args:
        hex_code (str): Hex color, with or without leading "#".

    Returns:
        tuple[float, float, float]: (r, g, b), each in [0, 1].
    """
    hex_code = hex_code.lstrip("#")
    return tuple(int(hex_code[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


# --------------------------------------------------------------------------
# Scene construction
# --------------------------------------------------------------------------
def setup_global_render_settings():
    """Resets the PyMOL session and applies the global render/style settings."""
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 0)
    cmd.set("ray_shadows", 1)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("specular", 0.25)
    cmd.set("ambient", 0.35)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("ray_trace_fog", 0)
    cmd.set("depth_cue", 0)


def load_and_style():
    """Loads the medoid structure and builds the base cartoon/sticks styling.

    Strips water/ions/hydrogens, then styles the protein cartoon, ligand
    sticks, and the 5 highlighted residues.

    Returns:
        list[str]: PyMOL selection names for the 5 highlighted residues,
        used later to build the camera framing selection.
    """
    cmd.load(PDB_PATH, OBJ_NAME)
    cmd.remove(f"{OBJ_NAME} and solvent")
    cmd.remove(f"{OBJ_NAME} and hydro")

    protein_sel = f"{OBJ_NAME} and chain {PROTEIN_CHAIN} and polymer"
    ligand_sel = f"{OBJ_NAME} and chain {LIGAND_CHAIN} and resn {LIGAND_RESN}"
    cmd.select(f"{OBJ_NAME}_ligand", ligand_sel)

    # base protein cartoon, neutral gray
    cmd.hide("everything", OBJ_NAME)
    cmd.show("cartoon", protein_sel)
    cmd.color(PROTEIN_GRAY, protein_sel)

    # Global cartoon transparency (applied to the whole protein cartoon,
    # not per-residue): with 5 highlighted residues spread around a
    # compact globular fold, at least one (e.g. resi 164, buried under a
    # beta-sheet) is not reachable by any single non-occluding camera
    # angle. Making the cartoon semi-transparent lets buried highlighted
    # sticks show through while the fold itself is still readable.
    cmd.set("cartoon_transparency", 0.5, protein_sel)

    # ligand: sticks, yellow-carbon scheme, always fully shown
    cmd.show("sticks", ligand_sel)
    cmd.set("stick_radius", 0.22, ligand_sel)
    util.cbay(ligand_sel)

    # 5 highlighted residues: colored cartoon segment + full-residue sticks
    # (backbone + sidechain, not just sidechain) in the same distinct color
    residue_sel_names = []
    for res in RESIDUES:
        resi = res["resi"]
        color_name = f"resi{resi}_color"
        cmd.set_color(color_name, list(hex_to_rgb01(res["hex"])))

        sel_name = f"{OBJ_NAME}_resi{resi}"
        cmd.select(
            sel_name,
            f"{OBJ_NAME} and chain {PROTEIN_CHAIN} and resi {resi}",
        )
        residue_sel_names.append(sel_name)

        cmd.color(color_name, sel_name)
        cmd.show("sticks", f"{sel_name} and not name H*")
        cmd.set("stick_radius", 0.2, sel_name)

    cmd.deselect()
    return residue_sel_names


def frame_camera(residue_sel_names):
    """Orients and zooms the single camera on the 5 residues plus ligand.

    Uses generous padding so nothing is cropped, since the 5 residues are
    spread across different parts of the pocket.

    Args:
        residue_sel_names (list[str]): Selection names from load_and_style.
    """
    zoom_sel = " or ".join(residue_sel_names) + f" or {OBJ_NAME}_ligand"
    cmd.orient(zoom_sel)
    cmd.zoom(zoom_sel, buffer=ZOOM_BUFFER)


def render():
    """Ray-traces the scene and writes OUT_PNG."""
    cmd.bg_color("white")
    cmd.ray(IMG_WIDTH, IMG_HEIGHT)
    cmd.png(OUT_PNG, dpi=PNG_DPI, ray=0)
    print(f"[five_residue_highlight] wrote {OUT_PNG}")


def main():
    """Builds the scene, frames the camera, and renders the figure."""
    os.makedirs(OUT_DIR, exist_ok=True)
    setup_global_render_settings()

    residue_sel_names = load_and_style()
    frame_camera(residue_sel_names)
    render()


main()
