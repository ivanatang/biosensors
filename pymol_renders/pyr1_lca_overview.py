#!/usr/bin/env python
"""Renders a "hero shot" of PYR1 (lca_001_binder) bound to LCA.

Publication-quality title/transition-slide image for a research talk, not
an analysis figure: unlike gate_latch_water_network.py /
medoid_comparison.py / five_residue_highlight.py, this script deliberately
adds no in-scene text labels, no distance lines, and no multi-panel
comparison -- just one clean, well-lit view of the bound complex. The gate
loop (resi 84-90) and latch loop (resi 114-118) still get distinct cartoon
colors (this project's central mechanistic story, and an all-gray cartoon
would undersell it), but purely as a visual accent -- no legend/labels are
drawn for them here.

Input is two PDB files sharing one coordinate frame (same MD complex,
already split into a protein file and a ligand file, both protonated):

    binders/lca_001_binder/protein_lca001_fixed_H.pdb  (chain A, resi 1-181)
    binders/lca_001_binder/ligand_lca001.pdb           (resn LIG, chain B)

Both load as separate PyMOL objects into the same session/scene; no
alignment step is needed since they already share coordinates from the
same source frame.

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/pyr1_lca_overview.py

Produces one ray-traced PNG in pymol_renders/output/:

    pyr1_lca_overview_hero.png

Background is plain white (cmd.bg_color("white")), per this project's
convention (see five_residue_highlight.py / medoid_comparison.py /
gate_latch_water_network.py).

Protein cartoon is pale gray-blue -- neutral, not one of CLAUDE.md's
categorical GROUP_COLOR entries, since that palette distinguishes
binder/nonbinder/etc. groups across figures, which doesn't apply to a
single-structure title slide. Gate (orange, #FE6100) and latch (purple,
#785EF0) use this repo's standard landmark colors (see CLAUDE.md "Key
domain conventions" and medoid_comparison.py) for consistency with the
rest of the deck, even though they aren't labeled here. The ligand is
sticks colored by element with green carbons -- reads clearly against
both the gray protein and the orange/purple landmarks, and is
deliberately distinct from the yellow-carbon convention used elsewhere in
pymol_renders/, so this title-slide ligand pops as the star of the image
rather than blending into the established analysis-figure visual
language.

Camera -- v2, what changed and why:

  v1 used a plain cmd.orient() on the gate+latch+ligand union, then
  cmd.zoom() on the whole protein+ligand with a wide buffer. The render
  looked face-on into the central beta-sheet, with its strands occluding
  most of the ligand, and the wide zoom made the ligand a small feature
  lost in the overall fold -- not usable as a hero shot.

  v2 fixes this two ways:

  1. Deterministic "look into the pocket opening" camera instead of a
     plain cmd.orient(). PYR1's pocket opens on the gate/latch end of the
     fold, so the camera should sit on the gate/latch side of the protein,
     looking back through the gate/latch into the ligand, not face-on to
     the beta-sheet core on the opposite side. Rather than guess-and-check
     with cmd.turn() (the FLIP_VERTICAL/TURN_Y_DEG pattern in
     gate_latch_water_network.py, which needs renders to tune against),
     the view is built directly from selection geometry:
       - protein_core_centroid = centroid of all protein CA atoms
       - gate_latch_centroid   = centroid of gate+latch CA atoms
       - outward_dir = normalize(gate_latch_centroid - protein_core_centroid),
         the direction from the protein's center of mass toward the
         gate/latch "bulge" -- the pocket-opening axis.
     A camera basis is built with the standard gluLookAt derivation
     (forward = -outward_dir, i.e. camera on the outward/gate-latch side
     looking back toward the ligand; right = normalize(cross(forward,
     up_ref)); up = cross(right, forward)) and written directly into
     PyMOL's view rotation (cmd.get_view()[0:9]) via cmd.set_view(). Only
     the rotation is set this way; cmd.zoom() (called right after, in
     frame_camera()) computes the translation/scale/clipping needed to fit
     the framed selection, since zoom preserves whatever rotation is
     already set and only dollies/recenters -- the same "set_view now,
     zoom later" pattern medoid_comparison.py uses for its shared camera.
     If gate/latch/ligand atoms are missing (shouldn't happen given
     load_structures()'s atom-count checks), this falls back to a plain
     cmd.orient() rather than crashing.
  2. Re-derived zoom target + occlusion mitigations:
       - zoom target is now "lig or gate or latch" (ligand + the two
         accent loops), not the whole protein, so the ligand fills a
         meaningful fraction of the frame instead of being a small feature
         in a wide whole-fold shot. Some surrounding protein cartoon is
         still visible from buffer padding and perspective falloff, so the
         shot still reads as "ligand nestled in a folded protein," just
         ligand-centric rather than fold-centric.
       - protein cartoon transparency raised to CARTOON_TRANSPARENCY so
         any sheet/strand still in front of part of the ligand lets the
         sticks show through rather than fully occluding them (same fix
         five_residue_highlight.py and gate_latch_water_network.py use for
         buried-residue/network occlusion; gate/latch are reset to fully
         opaque afterward so the accent colors stay vivid).
       - ligand stick_radius increased (LIGAND_STICK_RADIUS) so it reads
         clearly even where a strand edge crosses in front of it.

  v3 widened the zoom target back to the whole protein+ligand (see
  frame_camera()) so the shot reads as an establishing/title-slide view of
  the whole fold with the ligand visible in the pocket, not a tight
  mechanistic close-up -- v2's lig+gate+latch-only zoom cropped out most
  of the protein. The v2 pocket-opening rotation still points the
  ligand-facing side toward the viewer, so widening the zoom target keeps
  the ligand visible while restoring whole-fold context.

  v4 removed gate/latch accent highlighting from the cartoon coloring (see
  style_scene()); the "gate"/"latch" selections are kept only so
  frame_camera()'s pocket-opening axis calculation still has them.

Re-running after seeing an actual render: if the deterministic axis above
isn't quite right (e.g. up/down or left/right flipped), the cheapest fix
is a fixed cmd.turn("y", <deg>) / cmd.turn("z", <deg>) tweak in
frame_camera() after set_pocket_opening_view() -- the same
one-constant-toggle pattern as FLIP_VERTICAL/TURN_Y_DEG in
gate_latch_water_network.py, left out here since it needs tuning against
a real render.
"""

import os

from pymol import cmd, util

# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
REPO_ROOT = "/Users/ivanatang/Developer/biosensors"
OUT_DIR = os.path.join(REPO_ROOT, "pymol_renders", "output")

PROTEIN_PDB = os.path.join(
    REPO_ROOT, "binders/lca_001_binder/protein_lca001_fixed_H.pdb"
)
LIGAND_PDB = os.path.join(
    REPO_ROOT, "binders/lca_001_binder/ligand_lca001.pdb"
)
OUT_PNG = os.path.join(OUT_DIR, "pyr1_lca_overview_hero.png")

PROTEIN_OBJ = "pyr1"
LIGAND_OBJ = "lca_ligand"

PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
LIGAND_RESN = "LIG"

# Structural landmarks (see CLAUDE.md "Key domain conventions")
GATE_RESI = "84-90"
LATCH_RESI = "114-118"

# Colors
PROTEIN_COLOR_HEX = "#648FFF"     # CLAUDE.md GROUP_COLOR "Binder" blue --
                                   # single flat color for the whole cartoon
                                   # (v4: gate/latch accent highlighting
                                   # removed; gate/latch selections are
                                   # still used below only to compute the
                                   # pocket-opening camera axis, not color)
LIGAND_CARBON_COLOR = "green"     # deliberately distinct from the
                                   # yellow-carbon convention used elsewhere
                                   # in pymol_renders/, see module docstring

# Occlusion mitigation (v2, see module docstring "CAMERA" section)
CARTOON_TRANSPARENCY = 0.28       # applied to the whole protein cartoon
LIGAND_STICK_RADIUS = 0.32        # bumped up from an earlier 0.22 so the
                                   # ligand pops even where a strand crosses
                                   # in front of it

IMG_WIDTH = 2400
IMG_HEIGHT = 1800
PNG_DPI = 300

# Camera framing knob -- tighter than v1's whole-protein zoom, now anchored
# on ligand+gate+latch (see module docstring point 2)
ZOOM_BUFFER = 5.0


# --------------------------------------------------------------------------
# Small color / vector helpers (pure python, no numpy dependency)
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


def _sub(a, b):
    """Returns the 3-vector difference a - b."""
    return tuple(a[i] - b[i] for i in range(3))


def _cross(a, b):
    """Returns the 3D cross product a x b."""
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a, b):
    """Returns the dot product of two 3-vectors."""
    return sum(a[i] * b[i] for i in range(3))


def _norm(v):
    """Normalizes a 3-vector; returns (0, 0, 1) if v is near-zero length."""
    length = sum(c * c for c in v) ** 0.5
    if length < 1e-6:
        return (0.0, 0.0, 1.0)
    return tuple(c / length for c in v)


def _centroid(model):
    """Returns the mean atomic coordinate of a PyMOL chempy model.

    Args:
        model: PyMOL chempy Indexed/Storable model (has a .atom list).

    Returns:
        tuple[float, float, float]: (x, y, z) centroid.
    """
    n = len(model.atom)
    sx = sum(a.coord[0] for a in model.atom) / n
    sy = sum(a.coord[1] for a in model.atom) / n
    sz = sum(a.coord[2] for a in model.atom) / n
    return (sx, sy, sz)


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


def load_structures():
    """Loads the protein and ligand PDBs as two separate objects, one session.

    They already share one coordinate frame (same MD complex, pre-split
    into two files), so no alignment step is needed for them to render
    together as a single bound complex.

    Raises:
        SystemExit: No protein polymer atoms, or no ligand atoms, loaded.
    """
    cmd.load(PROTEIN_PDB, PROTEIN_OBJ)
    cmd.load(LIGAND_PDB, LIGAND_OBJ)

    cmd.remove(f"{PROTEIN_OBJ} and hydro")
    cmd.remove(f"{PROTEIN_OBJ} and solvent")
    # Ligand hydrogens are kept (already protonated) -- a title slide of
    # "the ligand" benefits from a complete, chemically faithful stick
    # model, unlike a tight mechanistic close-up where H's just add clutter.

    n_protein = cmd.count_atoms(f"{PROTEIN_OBJ} and polymer")
    n_ligand = cmd.count_atoms(f"{LIGAND_OBJ} and resn {LIGAND_RESN}")
    print(
        f"[pyr1_lca_overview] loaded protein: {n_protein} atoms "
        f"(chain {PROTEIN_CHAIN} expected); ligand: {n_ligand} atoms "
        f"(resn {LIGAND_RESN})"
    )
    if n_protein == 0:
        raise SystemExit(
            f"[pyr1_lca_overview] ABORT: no protein polymer atoms loaded "
            f"from {PROTEIN_PDB}"
        )
    if n_ligand == 0:
        raise SystemExit(
            f"[pyr1_lca_overview] ABORT: no resn {LIGAND_RESN} atoms "
            f"loaded from {LIGAND_PDB}"
        )


def style_scene():
    """Styles the protein cartoon (flat color) and ligand (element-colored sticks)."""
    protein_sel = f"{PROTEIN_OBJ} and chain {PROTEIN_CHAIN} and polymer"
    gate_sel = f"{protein_sel} and resi {GATE_RESI}"
    latch_sel = f"{protein_sel} and resi {LATCH_RESI}"
    ligand_sel = f"{LIGAND_OBJ} and resn {LIGAND_RESN}"

    cmd.select("gate", gate_sel)
    cmd.select("latch", latch_sel)
    cmd.select("lig", ligand_sel)

    # Base protein cartoon, single flat color across the whole chain (v4:
    # no gate/latch accent highlighting -- "gate"/"latch" selections above
    # are kept only so frame_camera()'s pocket-opening axis calculation
    # still has them to compute against).
    cmd.hide("everything", PROTEIN_OBJ)
    cmd.show("cartoon", protein_sel)
    cmd.set_color("protein_color", list(hex_to_rgb01(PROTEIN_COLOR_HEX)))
    cmd.color("protein_color", protein_sel)

    # Occlusion mitigation (v2): make the whole protein cartoon translucent
    # so ligand sticks that end up behind a beta strand still read through.
    cmd.set("cartoon_transparency", CARTOON_TRANSPARENCY, protein_sel)

    # ligand: full sticks, colored by element, green carbons so it pops
    # as the visual focal point of the shot. Stick radius bumped up (v2)
    # so it stays legible even where a strand edge crosses in front.
    cmd.hide("everything", LIGAND_OBJ)
    cmd.show("sticks", ligand_sel)
    cmd.set("stick_radius", LIGAND_STICK_RADIUS, ligand_sel)
    util.cbag(ligand_sel)  # color-by-atom with green carbons

    cmd.deselect()


def set_pocket_opening_view():
    """Builds a camera rotation looking into the pocket from the gate/latch end.

    Rather than face-on into the opposite beta-sheet core (see module
    docstring "Camera" section for why plain cmd.orient() picked the wrong
    side in v1). outward_dir = normalize(gate_latch_centroid -
    protein_core_centroid) is the direction from the protein's overall
    center of mass toward the gate/latch bulge, i.e. the pocket-opening
    axis. A camera is placed conceptually on the outward_dir side, looking
    back in (standard gluLookAt basis construction), and only the 3x3
    rotation block of PyMOL's view matrix (cmd.get_view()[0:9]) is
    overwritten via cmd.set_view() -- translation/scale/clipping are left
    for cmd.zoom() to (re)compute in frame_camera(), since zoom preserves
    whatever rotation is already set and only dollies/recenters (same
    "set_view now, zoom later" pattern medoid_comparison.py uses).

    Returns:
        bool: True if the deterministic view was set; False if it fell
        back (missing atoms in one of the three selections it needs --
        shouldn't happen given load_structures()'s atom-count checks, but
        guarded rather than assumed).
    """
    protein_ca_sel = (
        f"{PROTEIN_OBJ} and chain {PROTEIN_CHAIN} and polymer and name CA"
    )
    gate_latch_ca_sel = f"({protein_ca_sel}) and (resi {GATE_RESI} or resi {LATCH_RESI})"

    protein_model = cmd.get_model(protein_ca_sel)
    gate_latch_model = cmd.get_model(gate_latch_ca_sel)

    if not protein_model.atom or not gate_latch_model.atom:
        print(
            "[pyr1_lca_overview] WARNING: could not compute deterministic "
            "pocket-opening camera (missing CA atoms in protein or "
            "gate+latch selection); falling back to cmd.orient()"
        )
        return False

    protein_core_centroid = _centroid(protein_model)
    gate_latch_centroid = _centroid(gate_latch_model)

    outward_dir = _norm(_sub(gate_latch_centroid, protein_core_centroid))
    print(
        f"[pyr1_lca_overview] pocket-opening axis "
        f"(protein core -> gate/latch): {outward_dir}"
    )

    # Reference "up" vector for the look-at basis -- arbitrary world axis,
    # just needs to not be near-parallel to outward_dir (checked below).
    up_ref = (0.0, 0.0, 1.0)
    if abs(_dot(up_ref, outward_dir)) > 0.9:
        up_ref = (0.0, 1.0, 0.0)

    # Standard gluLookAt-style basis: forward = direction the camera looks
    # (from the outward/gate-latch side back toward the pocket core), right
    # = cross(forward, up_ref), up = cross(right, forward). PyMOL's view
    # rotation rows are (right, up, out-of-screen-toward-camera); the third
    # row (out-of-screen-toward-camera) equals -forward = outward_dir by
    # construction.
    forward = tuple(-c for c in outward_dir)
    right = _norm(_cross(forward, up_ref))
    up = _cross(right, forward)

    view = list(cmd.get_view())
    view[0:3] = right
    view[3:6] = up
    view[6:9] = outward_dir
    cmd.set_view(view)
    return True


def frame_camera():
    """Sets the pocket-opening camera rotation, then zooms on the whole complex.

    Uses set_pocket_opening_view(), falling back to cmd.orient(zoom_sel) if
    it can't be computed. Zoom target is the whole complex (v3, see module
    docstring), not just lig+gate+latch (v2 cropped out most of the
    protein) -- the pocket-opening rotation still points the ligand-facing
    side toward the viewer, so widening the zoom target keeps the ligand
    clearly visible while restoring whole-fold context.
    """
    zoom_sel = "lig or gate or latch"
    got_deterministic_view = set_pocket_opening_view()
    if not got_deterministic_view:
        cmd.orient(zoom_sel)

    cmd.zoom(f"{PROTEIN_OBJ} or {LIGAND_OBJ}", buffer=ZOOM_BUFFER)


def render():
    """Ray-traces the scene and writes OUT_PNG."""
    cmd.bg_color("white")
    cmd.ray(IMG_WIDTH, IMG_HEIGHT)
    cmd.png(OUT_PNG, dpi=PNG_DPI, ray=0)
    print(f"[pyr1_lca_overview] wrote {OUT_PNG}")


def main():
    """Runs the full pipeline: load, style, frame camera, render."""
    os.makedirs(OUT_DIR, exist_ok=True)
    setup_global_render_settings()

    load_structures()
    style_scene()
    frame_camera()
    render()


main()
