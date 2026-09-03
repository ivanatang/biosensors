#!/usr/bin/env python
"""Renders one real MD frame showing a gate-latch-ligand bridging water.

Publication-quality figure of one real MD frame (pair_3085_binder) showing
a single bridging water molecule simultaneously in contact with the gate
loop (resi 84-90), the latch loop (resi 114-118), and the ligand's
carboxylate/core oxygen, for a research talk figure.

Mechanistic story being illustrated: this repo's water-mediated contact
analysis (water_analysis/gate_latch_water_bridge.py) found that a single
water bridging gate + latch + ligand simultaneously occurs significantly
more often in true binders than false positives. This script renders one
concrete, real example frame of that geometry (not synthetic/idealized),
extracted from an actual trajectory, with the network made explicit via
dashed distance lines from the bridging water to (a) the nearest gate
heavy atom, (b) the nearest latch heavy atom, and (c) the specific
bridging ligand oxygen -- all three determined by selection query at
render time, not assumed in advance (see find_nearest_landmark_atom; the
actual closest atoms are printed to stdout when this script runs).

Input is a standalone single-frame PDB with protein + ligand + exactly one
water, all in chain A (components are distinguished by resn/resi, not
chain -- don't add chain-based selections for ligand/water). The current
input (v4, see PDB_PATH below) re-extracts a fix for a periodic-boundary
"molecule not whole" issue in the protein chain that an earlier v3 version
had (v3 only had the water's position PBC-fixed, not the protein, which
produced spurious long "bonds" between sequence-adjacent atoms wrapped to
different periodic images -- fixed via mdtraj's
image_molecules(make_whole=True), see the PDB_PATH comment below).
Topology/atom order is unchanged v3->v4; only coordinates were re-imaged.

The ligand (resn LIG, resi 182) has generic atom names (all "C"/"O", not
unique per atom), so the specific bridging ligand oxygen is selected by
PDB atom serial number (`id 2873`), not by name. This is verified at
runtime (see verify_bridge_oxygen) before any styling happens; the script
aborts with a clear error rather than guessing if that atom isn't a single
LIG oxygen.

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/gate_latch_water_network.py

Produces one ray-traced PNG in pymol_renders/output/:

    pair_3085_binder_water_network.png

Background is plain white (cmd.bg_color("white")), per this project's
convention (see five_residue_highlight.py / medoid_comparison.py).

Gate (orange, #FE6100) and latch (purple, #785EF0) use this repo's
standard landmark colors (see CLAUDE.md "Key domain conventions" and
medoid_comparison.py). The ligand is full yellow-carbon sticks (util.cbay)
including the bridging oxygen (id 2873), which keeps its default canonical
oxygen color (red) rather than a custom highlight (an earlier pink
override read as confusing, reverted); it's still identifiable as "the
interacting atom" by being rendered as an enlarged sphere, distinguished
by size/position rather than color. The bridging water (resn HOH, resi
2269) is an enlarged red/white ball-and-stick (radii bumped up after an
earlier render showed it partly lost against nearby side-chain sticks) so
it can't be mistaken for background noise (there's no other water in this
file, since only this one bridging water was extracted into the PDB).

Occlusion: the base protein cartoon (protein_sel, which includes gate/
latch since they're part of the backbone) renders at cartoon_transparency
0.5 -- the same fix as five_residue_highlight.py's buried-residue problem
-- so cartoon/sticks in front of the water and its dashed network lines
don't block them; gate/latch cartoon is then reset to fully opaque (0.0)
so the translucency fix doesn't wash out those landmark colors.

Network lines: dashed cmd.distance lines run from the bridging water's
oxygen to the nearest gate heavy atom, nearest latch heavy atom, and the
bridging ligand oxygen -- this is what makes the figure read as "a
network," not three separately-colored blobs. The nearest gate/latch atom
is found by a within-3.5-A polar-atom (O/N) query, widened stepwise (6.0 A
polar, then any heavy atom) only if that query comes up empty; whichever
tier matched is printed. Per this project's convention (see
medoid_comparison.py docstring for why in-scene text labels were
abandoned in favor of pixel-space post-processing), the numeric distance
value labels PyMOL would normally draw on each dashed line are hidden --
they clutter this close a zoom without adding info beyond what's already
printed to stdout -- but the dashed lines stay, since those convey the
network geometry in-scene.

Camera: single close-up view, oriented (cmd.orient) and zoomed (cmd.zoom,
buffer=ZOOM_BUFFER) tightly on the union of gate, latch, the bridging
ligand oxygen, and the bridging water -- same "orient on the landmark
union, no explicit up/down logic" pattern as gate_latch_movie.py's
set_fixed_camera(). Since plain cmd.orient() isn't guaranteed to put
gate/latch at the top of frame, FLIP_VERTICAL (near the top of this file)
is a one-constant toggle that calls cmd.turn("x", 180) right after
orient/zoom if the render comes out upside-down/backwards -- deliberately
left as a manual toggle to check against the actual render, not solved
analytically.

No in-scene text labels ("Gate"/"Latch"/etc): this script only renders the
colored/connected structure. Region text labels are added afterward by a
separate pixel-space PIL post-processing script, per this project's
pattern (see five_residue_highlight.py / medoid_comparison.py docstrings).

To re-run on a different frame/water, edit PDB_PATH, BRIDGE_LIG_ATOM_ID,
and WATER_RESI below.
"""

import os

from pymol import cmd, util

# --------------------------------------------------------------------------
# Config - edit these to re-point the script at a different frame/water
# --------------------------------------------------------------------------
REPO_ROOT = "/Users/ivanatang/Developer/biosensors"
OUT_DIR = os.path.join(REPO_ROOT, "pymol_renders", "output")

PDB_PATH = os.path.join(
    REPO_ROOT,
    "pymol_renders/scratch_traj/pair_3085_binder_bridge_frame_v4.pdb",
)
# v4 fixes a periodic-boundary "molecule not whole" issue in the PROTEIN
# chain (the earlier v3 file only had the water's position fixed; PyMOL
# drew spurious long "bonds" between sequence-adjacent atoms wrapped to
# different periodic images). v4 was re-extracted with mdtraj's
# image_molecules(make_whole=True), verified 0 CA-CA gaps > 10 A, and the
# water-ligand bridge distance re-confirmed at 3.45 A. Topology/atom order
# (and therefore PDB serial numbers / BRIDGE_LIG_ATOM_ID / WATER_RESI below)
# is unchanged from v3 -- only coordinates were re-imaged.
OUT_PNG = os.path.join(OUT_DIR, "pair_3085_binder_water_network.png")

OBJ_NAME = "pair_3085_bridge"

# Structural landmarks (see CLAUDE.md "Key domain conventions")
GATE_RESI = "84-90"
LATCH_RESI = "114-118"

LIGAND_RESN = "LIG"
LIGAND_RESI = 182
BRIDGE_LIG_ATOM_ID = 2873  # PDB atom serial number of the bridging ligand O
# (Was mistakenly 2872 -- an off-by-one between the 0-based mdtraj array
# index used during frame extraction and the 1-based PDB serial number in
# the saved file. Serial 2872 is a DIFFERENT, unrelated ligand oxygen
# ~15.6 A from the water; 2873 is the true bridging atom at 3.45 A,
# confirmed against the periodic-aware mdtraj distance calculation.)

WATER_RESN = "HOH"
WATER_RESI = 2269

# Colors
PROTEIN_GRAY = "gray80"
GATE_COLOR_HEX = "#FE6100"     # orange, this repo's standard gate color
LATCH_COLOR_HEX = "#785EF0"    # purple, this repo's standard latch color

# Nearest-landmark-atom search cutoffs (Angstrom), see
# find_nearest_landmark_atom() docstring for the widening logic
POLAR_CUTOFF_A = 3.5
POLAR_CUTOFF_WIDENED_A = 6.0

IMG_WIDTH = 1800
IMG_HEIGHT = 1400
PNG_DPI = 300

# Camera framing knobs
ZOOM_BUFFER = 4.0  # padding (Angstrom) around the framed network selection

# cmd.orient() on the gate+latch+ligand+water union is not guaranteed to
# put gate/latch at the top of frame for every structure (it happened to
# for gate_latch_movie.py's pair_3101_binder selection, with no explicit
# "flip to top" logic). If this frame's orientation renders
# upside-down/backwards, set this to True to rotate the camera 180 degrees
# about x right after orient/zoom -- a one-constant toggle, not solved
# analytically.
FLIP_VERTICAL = True

# Extra rotation around the vertical (screen Y) axis, applied after the
# orient/zoom/flip above, to reveal which latch residue the water actually
# contacts (His115's ring was partly edge-on/overlapping the water in the
# straight-on view). Degrees, positive turns the view left. One-constant
# toggle, tuned empirically against the actual render, same pattern as
# FLIP_VERTICAL.
TURN_Y_DEG = -25


# --------------------------------------------------------------------------
# Small color / geometry helpers (pure python, no numpy dependency)
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


def _dist(a, b):
    """Returns the Euclidean distance between two 3-coordinates."""
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


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


def load_structure():
    """Loads the single-frame PDB and defines the protein/gate/latch/etc selections.

    Strips hydrogens from the protein and ligand (they clutter a close-up
    stick render), but keeps the water's hydrogens -- the bridging water's
    O and both H's need to be visible.
    """
    cmd.load(PDB_PATH, OBJ_NAME)
    cmd.remove(f"{OBJ_NAME} and hydro and not resn {WATER_RESN}")

    cmd.select(f"{OBJ_NAME}_protein", f"{OBJ_NAME} and polymer.protein")
    cmd.select(
        f"{OBJ_NAME}_gate", f"{OBJ_NAME}_protein and resi {GATE_RESI}"
    )
    cmd.select(
        f"{OBJ_NAME}_latch", f"{OBJ_NAME}_protein and resi {LATCH_RESI}"
    )
    cmd.select(
        f"{OBJ_NAME}_ligand",
        f"{OBJ_NAME} and resn {LIGAND_RESN} and resi {LIGAND_RESI}",
    )
    cmd.select(
        f"{OBJ_NAME}_bridge_o", f"{OBJ_NAME} and id {BRIDGE_LIG_ATOM_ID}"
    )
    cmd.select(
        f"{OBJ_NAME}_water",
        f"{OBJ_NAME} and resn {WATER_RESN} and resi {WATER_RESI}",
    )
    cmd.select(
        f"{OBJ_NAME}_water_o", f"{OBJ_NAME}_water and elem O"
    )


def verify_bridge_oxygen():
    """Verifies `id {BRIDGE_LIG_ATOM_ID}` is a single LIG oxygen, before styling.

    Raises:
        SystemExit: The id selects other than exactly one atom, or that
            atom isn't an oxygen belonging to resn LIG.
    """
    sel = f"{OBJ_NAME}_bridge_o"
    n = cmd.count_atoms(sel)
    if n != 1:
        raise SystemExit(
            f"[gate_latch_water_network] ABORT: `id {BRIDGE_LIG_ATOM_ID}` "
            f"selected {n} atoms, expected exactly 1. Check "
            f"BRIDGE_LIG_ATOM_ID / the input PDB."
        )
    model = cmd.get_model(sel)
    atom = model.atom[0]
    if atom.symbol.strip().upper() != "O" or atom.resn != LIGAND_RESN:
        raise SystemExit(
            f"[gate_latch_water_network] ABORT: id {BRIDGE_LIG_ATOM_ID} "
            f"is {atom.resn}/{atom.name} (element {atom.symbol!r}), "
            f"expected an oxygen belonging to resn {LIGAND_RESN}."
        )
    print(
        f"[gate_latch_water_network] verified bridging ligand atom: "
        f"id {BRIDGE_LIG_ATOM_ID} = {atom.resn}{atom.resi}/{atom.name} "
        f"(element O) -- OK"
    )


def verify_water():
    """Checks the bridging water resolves to exactly one O + two H, printing a warning if not."""
    n_total = cmd.count_atoms(f"{OBJ_NAME}_water")
    n_o = cmd.count_atoms(f"{OBJ_NAME}_water and elem O")
    n_h = cmd.count_atoms(f"{OBJ_NAME}_water and elem H")
    print(
        f"[gate_latch_water_network] bridging water resn {WATER_RESN} "
        f"resi {WATER_RESI}: {n_total} atoms total ({n_o} O, {n_h} H)"
    )
    if n_total != 3 or n_o != 1 or n_h != 2:
        print(
            "[gate_latch_water_network] WARNING: expected exactly 1 O + "
            "2 H for the bridging water; got something else. Proceeding, "
            "but double-check the render."
        )


def style_scene():
    """Styles the protein cartoon, gate/latch, ligand, and bridging water.

    Base protein cartoon (neutral gray) + gate/latch highlight (cartoon +
    full-residue sticks) + ligand sticks with the bridging oxygen
    distinguished + the bridging water as an enlarged ball-and-stick.
    """
    protein_sel = f"{OBJ_NAME}_protein"
    gate_sel = f"{OBJ_NAME}_gate"
    latch_sel = f"{OBJ_NAME}_latch"
    ligand_sel = f"{OBJ_NAME}_ligand"
    bridge_o_sel = f"{OBJ_NAME}_bridge_o"
    water_sel = f"{OBJ_NAME}_water"

    cmd.hide("everything", OBJ_NAME)

    # Base protein cartoon, neutral gray, nothing else highlighted.
    # Semi-transparent (same fix as five_residue_highlight.py's buried-
    # residue problem): at this tight a zoom, cartoon/sticks in front of the
    # water and its dashed distance lines occluded the network in an
    # earlier render.
    cmd.show("cartoon", protein_sel)
    cmd.color(PROTEIN_GRAY, protein_sel)
    cmd.set("cartoon_transparency", 0.5, protein_sel)

    # gate / latch: distinct colors, thicker loop tube + full-residue
    # sticks (backbone + sidechain) so whichever atom ends up nearest the
    # water (which may be a backbone carbonyl O, not just a sidechain) is
    # visible
    cmd.set_color("gate_color", list(hex_to_rgb01(GATE_COLOR_HEX)))
    cmd.set_color("latch_color", list(hex_to_rgb01(LATCH_COLOR_HEX)))

    cmd.color("gate_color", gate_sel)
    cmd.color("latch_color", latch_sel)
    # Keep gate/latch cartoon fully opaque even though the rest of the
    # protein is now translucent (protein_sel above includes gate/latch
    # since they're part of the protein backbone) -- they're landmark
    # colors, not part of the occlusion problem being fixed.
    cmd.set("cartoon_transparency", 0.0, f"{gate_sel} or {latch_sel}")

    cmd.set("cartoon_loop_radius", 0.4, gate_sel)
    cmd.set("cartoon_loop_radius", 0.4, latch_sel)
    cmd.set("cartoon_tube_radius", 0.4, gate_sel)
    cmd.set("cartoon_tube_radius", 0.4, latch_sel)

    cmd.show("sticks", f"({gate_sel} or {latch_sel}) and not hydro")
    cmd.set("stick_radius", 0.2, f"{gate_sel} or {latch_sel}")

    # ligand: yellow-carbon sticks throughout, EXCEPT the bridging oxygen
    cmd.show("sticks", ligand_sel)
    cmd.set("stick_radius", 0.22, ligand_sel)
    util.cbay(ligand_sel)

    # Left at its default canonical oxygen color (red, from util.cbay above)
    # rather than a custom highlight color -- an earlier pink override here
    # read as confusing rather than clarifying. Still shown as an enlarged
    # sphere below so it's identifiable as "the specific interacting atom"
    # by size/position, not by an unusual color.
    cmd.show("spheres", bridge_o_sel)
    cmd.set("sphere_scale", 0.35, bridge_o_sel)

    # bridging water: enlarged red O / white H ball-and-stick, deliberately
    # more prominent than default atom scale so it can't be overlooked.
    # Bumped up further (was 0.25/0.3) after an earlier render showed the
    # water partly lost against nearby side-chain sticks even with the
    # cartoon now translucent.
    cmd.show("sticks", water_sel)
    cmd.show("spheres", water_sel)
    cmd.set("stick_radius", 0.3, water_sel)
    cmd.set("sphere_scale", 0.4, water_sel)
    cmd.color("red", f"{water_sel} and elem O")
    cmd.color("white", f"{water_sel} and elem H")

    cmd.deselect()


# --------------------------------------------------------------------------
# Network geometry: find nearest gate/latch atom, draw dashed distances
# --------------------------------------------------------------------------
def find_nearest_landmark_atom(region_sel, region_label, water_o_coord):
    """Finds the landmark (gate or latch) atom nearest the bridging water's oxygen.

    Doesn't assume which atom it is in advance. Tries, in order: (1) polar
    heavy atoms (O/N) within POLAR_CUTOFF_A of the water O, (2) polar heavy
    atoms (O/N) within the widened POLAR_CUTOFF_WIDENED_A, (3) any heavy
    (non-H) atom in the region regardless of distance -- uses the first
    non-empty tier.

    Args:
        region_sel (str): PyMOL selection for the gate or latch region.
        region_label (str): "gate" or "latch", used in print statements.
        water_o_coord (list[float]): [x, y, z] of the bridging water's oxygen.

    Returns:
        tuple: (atom, distance) for the nearest matching atom.

    Raises:
        SystemExit: `region_sel` contains no atoms at all.
    """
    tiers = [
        (
            f"({region_sel}) and (elem O+N) within {POLAR_CUTOFF_A} of "
            f"({OBJ_NAME}_water_o)",
            f"polar (O/N) atom within {POLAR_CUTOFF_A} A",
        ),
        (
            f"({region_sel}) and (elem O+N) within "
            f"{POLAR_CUTOFF_WIDENED_A} of ({OBJ_NAME}_water_o)",
            f"polar (O/N) atom within widened {POLAR_CUTOFF_WIDENED_A} A "
            f"(none found within {POLAR_CUTOFF_A} A)",
        ),
        (
            f"({region_sel}) and not hydro",
            f"nearest heavy atom, unrestricted distance (no polar atom "
            f"found within {POLAR_CUTOFF_WIDENED_A} A)",
        ),
    ]

    for candidate_sel, note in tiers:
        if cmd.count_atoms(candidate_sel) == 0:
            continue
        model = cmd.get_model(candidate_sel)
        best_atom, best_d = None, None
        for a in model.atom:
            d = _dist(a.coord, water_o_coord)
            if best_d is None or d < best_d:
                best_atom, best_d = a, d
        print(
            f"[gate_latch_water_network] nearest {region_label} atom to "
            f"bridging water O: {best_atom.resn}{best_atom.resi}/"
            f"{best_atom.name} at {best_d:.2f} A ({note})"
        )
        return best_atom, best_d

    raise SystemExit(
        f"[gate_latch_water_network] ABORT: found no atoms at all in "
        f"region {region_label!r} ({region_sel}); check the resi range."
    )


def draw_network_distances():
    """Draws dashed distance lines from the bridging water to gate/latch/ligand.

    Returns:
        tuple: (gate_atom, latch_atom) chempy Atom objects, for use by
        place_position_markers.
    """
    water_o_sel = f"{OBJ_NAME}_water_o"
    water_o_coord = cmd.get_model(water_o_sel).atom[0].coord

    gate_atom, _ = find_nearest_landmark_atom(
        f"{OBJ_NAME}_gate", "gate", water_o_coord
    )
    latch_atom, _ = find_nearest_landmark_atom(
        f"{OBJ_NAME}_latch", "latch", water_o_coord
    )

    gate_atom_sel = f"{OBJ_NAME} and id {gate_atom.id}"
    latch_atom_sel = f"{OBJ_NAME} and id {latch_atom.id}"
    bridge_o_sel = f"{OBJ_NAME}_bridge_o"

    dist_names = [
        cmd.distance(f"{OBJ_NAME}_dist_gate", water_o_sel, gate_atom_sel),
        cmd.distance(f"{OBJ_NAME}_dist_latch", water_o_sel, latch_atom_sel),
        cmd.distance(f"{OBJ_NAME}_dist_ligand", water_o_sel, bridge_o_sel),
    ]
    print(
        f"[gate_latch_water_network] water O -> bridging ligand O "
        f"(id {BRIDGE_LIG_ATOM_ID}) distance = {dist_names[2]:.2f} A"
    )

    for dist_obj in (
        f"{OBJ_NAME}_dist_gate",
        f"{OBJ_NAME}_dist_latch",
        f"{OBJ_NAME}_dist_ligand",
    ):
        # Keep the dashed lines (they convey the network geometry) but hide
        # the numeric value labels -- clutter at this close a zoom without
        # adding info beyond what's printed to stdout above (see module
        # docstring for the project's in-scene-label convention).
        # dash_radius/dash_gap tightened (was width-only at 4/0.3) after an
        # earlier render showed the lines blending into nearby side-chain
        # sticks -- thicker, more continuous dashes read more clearly as
        # "connections" at this zoom.
        cmd.hide("labels", dist_obj)
        cmd.color("gray20", dist_obj)
        cmd.set("dash_width", 5, dist_obj)
        cmd.set("dash_radius", 0.08, dist_obj)
        cmd.set("dash_gap", 0.15, dist_obj)
        cmd.set("dash_length", 0.15, dist_obj)

    return gate_atom, latch_atom


def place_position_markers(gate_atom, latch_atom):
    """Places cyan marker pseudoatoms at the found contact-residue positions.

    Used only to locate where to put text labels in a separate pixel-space
    post-processing pass (same HSV-color-detection technique as
    add_water_network_labels.py): rendered once for that purpose, then
    hidden before the real output PNG is rendered, so no dot/marker symbol
    appears in the final image, per this project's "text-only, no pointer
    symbols" labeling convention.

    Args:
        gate_atom: chempy Atom nearest the water on the gate, from
            draw_network_distances.
        latch_atom: chempy Atom nearest the water on the latch, from
            draw_network_distances.
    """
    cmd.pseudoatom("gate_contact_marker", pos=list(gate_atom.coord))
    cmd.pseudoatom("latch_contact_marker", pos=list(latch_atom.coord))
    cmd.set_color("marker_color", list(hex_to_rgb01("#00FFFF")))  # unique cyan
    cmd.color("marker_color", "gate_contact_marker or latch_contact_marker")
    cmd.show("spheres", "gate_contact_marker or latch_contact_marker")
    cmd.set("sphere_scale", 0.15, "gate_contact_marker or latch_contact_marker")


# --------------------------------------------------------------------------
# Camera + render
# --------------------------------------------------------------------------
def frame_camera():
    """Orients and zooms tightly on gate, latch, ligand oxygen, and water.

    Applies the FLIP_VERTICAL / TURN_Y_DEG manual toggles afterward if set.
    """
    network_sel = (
        f"{OBJ_NAME}_gate or {OBJ_NAME}_latch or {OBJ_NAME}_bridge_o "
        f"or {OBJ_NAME}_water"
    )
    cmd.orient(network_sel)
    cmd.zoom(network_sel, buffer=ZOOM_BUFFER)

    if FLIP_VERTICAL:
        cmd.turn("x", 180)

    if TURN_Y_DEG:
        cmd.turn("y", TURN_Y_DEG)


MARKER_PNG = os.path.join(OUT_DIR, "pair_3085_binder_water_network_MARKERS_tmp.png")


def render():
    """Renders two passes: a marker pass (for the labeling script) and the real output.

    Both passes are explicitly re-ray-traced at the same size. ray=0 in
    cmd.png() just dumps whatever is in the current offscreen buffer, and
    that buffer is invalidated by any scene change (e.g. the hide() below)
    -- skipping the second cmd.ray() silently falls back to the small
    default OpenGL viewport instead of IMG_WIDTH x IMG_HEIGHT. This bit us
    once (marker pass came out 1800x1400 but the real pass came out
    640x480, breaking the marker-to-real pixel coordinate mapping the
    labeling script depends on).
    """
    cmd.bg_color("white")

    # Pass 1: markers visible, for the labeling script to locate by color.
    cmd.ray(IMG_WIDTH, IMG_HEIGHT)
    cmd.png(MARKER_PNG, dpi=PNG_DPI, ray=0)
    print(f"[gate_latch_water_network] wrote marker-position pass {MARKER_PNG}")

    # Pass 2: hide markers, re-ray, render the real (clean) output.
    cmd.hide("everything", "gate_contact_marker or latch_contact_marker")
    cmd.ray(IMG_WIDTH, IMG_HEIGHT)
    cmd.png(OUT_PNG, dpi=PNG_DPI, ray=0)
    print(f"[gate_latch_water_network] wrote {OUT_PNG}")


def main():
    """Runs the full pipeline: load, verify, style, network lines, camera, render."""
    os.makedirs(OUT_DIR, exist_ok=True)
    setup_global_render_settings()

    load_structure()
    verify_bridge_oxygen()
    verify_water()

    style_scene()
    gate_atom, latch_atom = draw_network_distances()
    frame_camera()
    place_position_markers(gate_atom, latch_atom)
    render()

    print(
        f"[gate_latch_water_network] contact residues: "
        f"gate={gate_atom.resn}{gate_atom.resi}  latch={latch_atom.resn}{latch_atom.resi}"
    )


main()
