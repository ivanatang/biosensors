#!/usr/bin/env python
"""Renders a rich, poster-quality "hero shot" of PYR1 (lca_001_binder) bound to LCA.

This is a showcase/gift image for a coworker, not an analysis figure and not
a title slide -- unlike pyr1_lca_overview.py (flat pale-blue cartoon + green
sticks on plain white, deliberately restrained for a talk transition slide),
this script leans into visual complexity: a cartoon + translucent molecular
surface hybrid, a dark gradient backdrop, richer lighting, and a bolder
ligand treatment, so it reads as a polished "hero render" rather than a lab
figure. It is a new, independent file; pyr1_lca_overview.py and its output
PNG are left untouched. The small vector/color helpers and the
pocket-opening camera derivation are intentionally duplicated here (not
imported) so this script has no load-bearing dependency on the other file's
internals.

Input is the same two coordinate-matched PDBs used by pyr1_lca_overview.py
(same MD complex frame, already split into protein and ligand files, both
protonated):

    binders/lca_001_binder/protein_lca001_fixed_H.pdb  (chain A, resi 1-181)
    binders/lca_001_binder/ligand_lca001.pdb           (resn LIG, chain B)

They load as two separate PyMOL objects sharing one coordinate frame; no
alignment step is needed.

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/pyr1_lca_hero_complex.py

Produces one ray-traced PNG in pymol_renders/output/:

    pyr1_lca_biosensor_hero.png

Styling choices:

  Background -- a dark navy-to-near-black vertical gradient
  (bg_gradient/bg_rgb_top/bg_rgb_bottom), not this repo's usual plain white.
  CLAUDE.md's plain-white convention is for analysis figures meant to sit in
  a paper/notebook grid; this is explicitly a one-off showcase image, and a
  dark backdrop is standard practice for molecular "hero shots" because it
  lets ray-traced specular highlights and the translucent surface shell read
  clearly (both wash out against white).

  Protein -- pale cool grey-blue cartoon (#D7DEE8, close to white so it
  reads clearly against the dark background without competing with the
  ligand or the gate/latch accents) wrapped in a translucent molecular
  surface (surface_transparency 0.7) built over the *same* per-atom colors
  as the cartoon. Because surface and cartoon share atom colors, the gate
  loop (resi 84-90, #FE6100 orange) and latch loop (resi 114-118, #785EF0
  purple) -- this project's central mechanistic story, see CLAUDE.md "Key
  domain conventions" -- show through as colored patches on the glassy
  outer shell as well as on the ribbon beneath it, so the accent reads at
  both structural levels instead of only on the cartoon.

  Ligand -- sticks plus small spheres (a classic "ball-and-stick" look,
  richer than sticks alone) with a custom bright cyan carbon color
  (#00E5FF, set via cmd.color + util.cnc rather than a stock util.cb*
  preset so the hex is exact) chosen to sit far from both the warm
  gate-orange and latch-purple accents on the color wheel, so the ligand
  unambiguously reads as the focal point rather than blending with the
  mechanistic accents.

  Lighting -- ray_shadows on, three-point-style light_count, moderate
  specular/reflect/shininess for a glossy-but-not-plastic look, perspective
  camera (orthoscopic off, unlike pyr1_lca_overview.py's orthoscopic=1) for
  a more dynamic hero-shot feel, and a subtle depth_cue/fog matched to the
  background's dark end so depth falls off gently toward the back of the
  frame.

Camera -- reuses pyr1_lca_overview.py's deterministic "look into the pocket
opening" derivation (protein-core-centroid -> gate/latch-centroid outward
axis, gluLookAt-style basis written into cmd.get_view()[0:9] via
cmd.set_view(), then cmd.zoom() computes translation/scale/clipping), since
that technique is already validated against a real render and PYR1's pocket
reliably opens on the gate/latch side of the fold. See
set_pocket_opening_view() below.

Iteration log (recorded here per CLAUDE.md's "why" documentation convention,
same reasoning as pyr1_lca_overview.py's v1-v4 history -- first-pass camera
and occlusion choices are rarely right until checked against a real render):

  v1: initial version as described above -- dark gradient background,
  cartoon+surface hybrid with shared gate/latch accent coloring, cyan
  ball-and-stick ligand, perspective camera reusing the pocket-opening axis,
  zoom target the whole protein+ligand union with a generous buffer.

  v2: two fixes after inspecting the v1 render. (1) The surface shell
  shared per-atom coloring with the cartoon, so the gate/latch accent was
  rendered twice (opaque cartoon color seen through a translucent surface
  layer of the same hue), reading as a muddy brown/liver tone instead of a
  clean orange/purple accent -- fixed by decoupling surface color from
  cartoon color via PyMOL's surface_color setting (SURFACE_HEX, a uniform
  pale silvery-blue "glass" tone), so gate/latch stay crisp on the ribbon
  while the shell reads as neutral. (2) ZOOM_BUFFER of 4.0 left the
  molecule small with a lot of dead gradient background around it --
  tightened to 1.5. Ligand stick/sphere size also bumped slightly so it
  holds its own as the focal point.

  v3: the v2 fix made the accent colors clean, but the surface shell still
  rendered as dark charcoal rather than glass -- at SURFACE_TRANSPARENCY
  0.7 over a near-black background, only ~30% of the surface's own color
  reaches the eye, so v2's merely-pale SURFACE_HEX wasn't light enough.
  Pushed SURFACE_HEX much lighter and brought SURFACE_TRANSPARENCY down to
  0.55, and nudged ambient up slightly, so the shell reads as pale glass
  rather than rock while staying translucent enough to see the ribbon and
  ligand through it.

  v4: switched the background from v1-v3's dark indigo-to-black gradient to
  plain white, per user request, matching this repo's usual analysis-figure
  convention after all. This meant re-tuning BASE_CARTOON_HEX and
  SURFACE_HEX, both of which had been pushed pale specifically to read
  against a dark backdrop (see v1, v3) -- unchanged, they would have nearly
  vanished against white. Both were darkened to medium-saturation blues
  that keep enough contrast on white while staying visually distinct from
  the gate/latch accent colors. BG_TOP_HEX/BG_BOTTOM_HEX (the gradient
  endpoints) were replaced by a single flat BG_HEX.

  v5: latch accent changed from CLAUDE.md's usual landmark purple
  (#785EF0) to pink (#DC267F), by request, for this render only -- see
  LATCH_HEX.
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
OUT_PNG = os.path.join(OUT_DIR, "pyr1_lca_biosensor_hero.png")

PROTEIN_OBJ = "pyr1"
LIGAND_OBJ = "lca_ligand"

PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
LIGAND_RESN = "LIG"

# Structural landmarks (see CLAUDE.md "Key domain conventions")
GATE_RESI = "84-90"
LATCH_RESI = "114-118"

# Colors
BASE_CARTOON_HEX = "#5B7FBE"   # v4: darkened from #D7DEE8 (v1-v3, a pale
                                # grey-blue chosen to read against a dark
                                # background) now that the background is
                                # white -- that pale a color would nearly
                                # vanish against white, so this uses a
                                # medium-saturation blue close to CLAUDE.md's
                                # "Binder" categorical color (#648FFF) for
                                # enough contrast while staying visually
                                # distinct from the warm gate/latch accents
GATE_HEX = "#FE6100"           # CLAUDE.md gate loop landmark color
LATCH_HEX = "#DC267F"          # v5: pink, by request, for this render only
                                # -- CLAUDE.md's landmark color for latch is
                                # purple (#785EF0); this hex is instead
                                # CLAUDE.md's GROUP_COLOR "False Positive"
                                # pink, reused here for its IBM
                                # colorblind-safe palette membership even
                                # though that entry's usual meaning
                                # (classifier false positives) doesn't apply
                                # to this structural-landmark context
LIGAND_CARBON_HEX = "#00E5FF"  # bright cyan, chosen to sit far from the
                                # warm gate/latch hues so the ligand reads
                                # as the unambiguous focal point
SURFACE_HEX = "#B9CBE8"        # v4: darkened from #EAF1FA (v3, tuned to be
                                # visible against near-black) -- that near-
                                # white shell would be nearly invisible
                                # against a white background, so this uses a
                                # pale but clearly-toned blue that still
                                # reads as a translucent "glass" halo around
                                # the cartoon rather than disappearing into
                                # the page

BG_HEX = "#FFFFFF"             # v4: plain white, matching this repo's
                                # analysis-figure convention (see CLAUDE.md
                                # and pyr1_lca_overview.py) -- v1-v3 used a
                                # dark gradient instead; no longer gradient,
                                # so a single flat color replaces
                                # BG_TOP_HEX/BG_BOTTOM_HEX

# Surface / cartoon hybrid
SURFACE_TRANSPARENCY = 0.55    # v3: brought down from 0.7 (v1) -- combined
                                # with the brighter SURFACE_HEX, lets enough
                                # of the shell's own pale color reach the
                                # eye to read as glass rather than dark
                                # charcoal, while still staying translucent
                                # enough for the cartoon ribbon and ligand
                                # to show through clearly
CARTOON_TRANSPARENCY = 0.0     # cartoon itself stays fully opaque; only the
                                # surface layer is translucent, so ribbon
                                # detail stays crisp under the shell

# Ligand geometry
LIGAND_STICK_RADIUS = 0.3      # v2: bumped from 0.26 so the ligand holds
                                # its own as the focal point against the
                                # now-busier surface+cartoon shell
LIGAND_SPHERE_SCALE = 0.26     # small spheres at atom centers on top of
                                # sticks -- classic "ball-and-stick" look

IMG_WIDTH = 2400
IMG_HEIGHT = 1800
PNG_DPI = 300

# Camera framing
ZOOM_BUFFER = 1.5              # v2: tightened from 4.0 -- v1's wide buffer
                                # left the molecule small in frame with a lot
                                # of dead gradient background around it

# Camera fine-tune knobs (see module docstring "Iteration log"; left at 0
# for v1, adjusted after inspecting a real render if the deterministic
# pocket-opening axis needs a small correction).
TURN_Y_DEG = 0
TURN_X_DEG = 0


# --------------------------------------------------------------------------
# Small color / vector helpers (pure python, no numpy dependency)
#
# Duplicated from pyr1_lca_overview.py rather than imported -- this script
# is meant to stand alone, see module docstring.
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
    """Resets the PyMOL session and applies global render, lighting, and background settings.

    Perspective camera (orthoscopic off) is a deliberate departure from this
    repo's orthoscopic analysis-figure convention (see pyr1_lca_overview.py)
    -- this is a one-off showcase render, not a figure meant to sit in a
    uniform grid, so a more dynamic/photographic camera is used instead. The
    background itself is plain white (v4; v1-v3 used a dark indigo-to-black
    gradient), matching this repo's usual convention after all -- see
    module docstring "v4" for why the base cartoon/surface colors also
    needed darkening to stay legible against it.
    """
    cmd.reinitialize()

    cmd.set_color("bg_flat", list(hex_to_rgb01(BG_HEX)))
    cmd.bg_color("bg_flat")
    cmd.set("bg_gradient", 0)
    cmd.set("ray_opaque_background", 1)

    cmd.set("ray_trace_mode", 0)
    cmd.set("ray_shadows", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 0)
    cmd.set("field_of_view", 28)

    # Lighting balance -- moderate specular/reflect for a glossy but not
    # plastic-looking surface shell; light_count > 1 softens shadow edges.
    cmd.set("light_count", 3)
    cmd.set("ambient", 0.28)
    cmd.set("direct", 0.65)
    cmd.set("reflect", 0.45)
    cmd.set("specular", 0.4)
    cmd.set("shininess", 60)
    cmd.set("spec_power", 150)
    cmd.set("spec_direct", 0.25)
    cmd.set("spec_direct_power", 55)

    # Subtle depth fog so distance falls off gently rather than looking flat
    # under a perspective camera. PyMOL has no separate fog_color setting;
    # fog always blends toward the current bg_color (white, set above).
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.5)

    # Transparent-surface correctness: real depth-sorted transparency
    # (mode 2) is needed for the translucent surface shell to composite
    # correctly against the opaque cartoon and ligand beneath it, and
    # two-sided lighting keeps the surface's inner faces lit properly where
    # it curves back on itself around the pocket.
    cmd.set("transparency_mode", 2)
    cmd.set("two_sided_lighting", 1)
    cmd.set("ray_interior_color", "grey20")

    cmd.set("surface_quality", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)


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
    # Ligand hydrogens are kept (already protonated) -- a hero shot of the
    # bound ligand benefits from a complete, chemically faithful model.

    n_protein = cmd.count_atoms(f"{PROTEIN_OBJ} and polymer")
    n_ligand = cmd.count_atoms(f"{LIGAND_OBJ} and resn {LIGAND_RESN}")
    print(
        f"[pyr1_lca_hero_complex] loaded protein: {n_protein} atoms "
        f"(chain {PROTEIN_CHAIN} expected); ligand: {n_ligand} atoms "
        f"(resn {LIGAND_RESN})"
    )
    if n_protein == 0:
        raise SystemExit(
            f"[pyr1_lca_hero_complex] ABORT: no protein polymer atoms "
            f"loaded from {PROTEIN_PDB}"
        )
    if n_ligand == 0:
        raise SystemExit(
            f"[pyr1_lca_hero_complex] ABORT: no resn {LIGAND_RESN} atoms "
            f"loaded from {LIGAND_PDB}"
        )


def style_scene():
    """Styles the protein as a cartoon+surface hybrid and the ligand as ball-and-stick.

    Cartoon uses per-atom coloring (base grey-blue, gate orange, latch
    purple) for the mechanistic accent. Surface uses a single uniform color
    (surface_color setting, decoupled from atom color) so the translucent
    shell reads as neutral glass rather than double-applying the gate/latch
    hues -- see module docstring "v2" iteration note.
    """
    protein_sel = f"{PROTEIN_OBJ} and chain {PROTEIN_CHAIN} and polymer"
    gate_sel = f"{protein_sel} and resi {GATE_RESI}"
    latch_sel = f"{protein_sel} and resi {LATCH_RESI}"
    ligand_sel = f"{LIGAND_OBJ} and resn {LIGAND_RESN}"

    cmd.select("gate", gate_sel)
    cmd.select("latch", latch_sel)
    cmd.select("lig", ligand_sel)

    # Shared atom coloring for cartoon + surface (both representations
    # inherit whatever color is set on the atom unless overridden per-rep,
    # so one set of color commands styles both layers at once).
    cmd.set_color("base_color", list(hex_to_rgb01(BASE_CARTOON_HEX)))
    cmd.set_color("gate_color", list(hex_to_rgb01(GATE_HEX)))
    cmd.set_color("latch_color", list(hex_to_rgb01(LATCH_HEX)))
    cmd.color("base_color", protein_sel)
    cmd.color("gate_color", gate_sel)
    cmd.color("latch_color", latch_sel)

    cmd.hide("everything", PROTEIN_OBJ)
    cmd.show("cartoon", protein_sel)
    cmd.set("cartoon_transparency", CARTOON_TRANSPARENCY, protein_sel)

    cmd.show("surface", protein_sel)
    cmd.set("transparency", SURFACE_TRANSPARENCY, protein_sel)
    cmd.set_color("surface_neutral", list(hex_to_rgb01(SURFACE_HEX)))
    cmd.set("surface_color", "surface_neutral", protein_sel)

    # Ligand: ball-and-stick with a custom cyan carbon color (cmd.color +
    # util.cnc rather than a stock util.cb* preset, since those only offer
    # a fixed palette of carbon colors and this hex is a deliberate choice,
    # see module docstring "Ligand").
    cmd.hide("everything", LIGAND_OBJ)
    cmd.show("sticks", ligand_sel)
    cmd.show("spheres", ligand_sel)
    cmd.set("stick_radius", LIGAND_STICK_RADIUS, ligand_sel)
    cmd.set("sphere_scale", LIGAND_SPHERE_SCALE, ligand_sel)
    cmd.set_color("lig_carbon_color", list(hex_to_rgb01(LIGAND_CARBON_HEX)))
    cmd.color("lig_carbon_color", ligand_sel)
    util.cnc(ligand_sel)  # color non-carbon elements by element, keep carbons cyan

    cmd.deselect()


def set_pocket_opening_view():
    """Builds a camera rotation looking into the pocket from the gate/latch end.

    Reused from pyr1_lca_overview.py (see that script's module docstring
    "Camera" section for the full derivation rationale and v1 failure mode
    it fixes -- a plain cmd.orient() looks face-on into the central
    beta-sheet and occludes the ligand). Summary: outward_dir =
    normalize(gate_latch_centroid - protein_core_centroid) is the
    pocket-opening axis; a gluLookAt-style basis is built from it and
    written into cmd.get_view()[0:9] via cmd.set_view(). Only rotation is
    set here; cmd.zoom() in frame_camera() computes translation/scale.

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
            "[pyr1_lca_hero_complex] WARNING: could not compute "
            "deterministic pocket-opening camera (missing CA atoms in "
            "protein or gate+latch selection); falling back to cmd.orient()"
        )
        return False

    protein_core_centroid = _centroid(protein_model)
    gate_latch_centroid = _centroid(gate_latch_model)

    outward_dir = _norm(_sub(gate_latch_centroid, protein_core_centroid))
    print(
        f"[pyr1_lca_hero_complex] pocket-opening axis "
        f"(protein core -> gate/latch): {outward_dir}"
    )

    up_ref = (0.0, 0.0, 1.0)
    if abs(_dot(up_ref, outward_dir)) > 0.9:
        up_ref = (0.0, 1.0, 0.0)

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

    Falls back to cmd.orient() on "lig or gate or latch" if the
    deterministic view can't be computed. TURN_Y_DEG/TURN_X_DEG are
    post-hoc fine-tune knobs (see module docstring "Iteration log") applied
    after the deterministic rotation but before zoom, in case the axis
    needs a small correction once checked against a real render.
    """
    zoom_sel = "lig or gate or latch"
    got_deterministic_view = set_pocket_opening_view()
    if not got_deterministic_view:
        cmd.orient(zoom_sel)

    if TURN_Y_DEG:
        cmd.turn("y", TURN_Y_DEG)
    if TURN_X_DEG:
        cmd.turn("x", TURN_X_DEG)

    cmd.zoom(f"{PROTEIN_OBJ} or {LIGAND_OBJ}", buffer=ZOOM_BUFFER)


def render():
    """Ray-traces the scene and writes OUT_PNG."""
    cmd.ray(IMG_WIDTH, IMG_HEIGHT)
    cmd.png(OUT_PNG, dpi=PNG_DPI, ray=0)
    print(f"[pyr1_lca_hero_complex] wrote {OUT_PNG}")


def main():
    """Runs the full pipeline: load, style, frame camera, render."""
    os.makedirs(OUT_DIR, exist_ok=True)
    setup_global_render_settings()

    load_structures()
    style_scene()
    frame_camera()
    render()


main()
