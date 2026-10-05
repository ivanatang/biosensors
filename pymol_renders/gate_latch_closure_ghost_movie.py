#!/usr/bin/env python
"""Renders the 3D half of a gate-closure movie with a fixed "ghost" of the
open-state gate and latch loops and a live gate-to-partner distance dash.

Which trajectory window and which distance are set in
gate_closure_movie_config.py (GATE_MOVIE=<key>; see that file for each
movie's data-selection rationale and how its window xtc was cut):
  - pair3087_lig3oh (default): Ca88 -> ligand 3-OH O, pair_3087_binder
    108-128 ns, the gate closing onto the ligand (~1.2 -> ~0.6 nm).
  - pair3088_ca116: Ca88 -> Ca116, pair_3088_binder 440-460 ns
    (~1.45 -> ~0.80 nm gate-latch closure).
pair_3101_binder (gate_latch_pocket_oriented_movie.py) was dropped as a
source: it only closes ~1.5 A on Ca88-Ca116 and reopens by ~35 ns.

The 3D readout is recomputed by PyMOL from the trajectory each frame. For
pair3088_ca116 it was checked against gate_latch116_timeseries.xvg (states
1/2/3/266 = 1.551/1.547/1.506/0.778 nm, matching xvg rows 2/4/6/532), which
also confirmed the load_traj frame offset described in the config.

Pipeline (two scripts; PyMOL's bundled python has no matplotlib):

    # 1. 3D frames, chunked to stay under this 8 GB machine's memory ceiling
    #    (see gate_latch_pocket_oriented_movie.py "Chunked rendering").
    for start in 1 91 181; do
      stop=$((start + 89)); [ $stop -gt 266 ] && stop=266
      GATE_MOVIE=pair3087_lig3oh GATE_LATCH_FRAME_START=$start GATE_LATCH_FRAME_STOP=$stop \
        /opt/homebrew/bin/pymol -cq pymol_renders/gate_latch_closure_ghost_movie.py
    done
    # 2. Trace strip + vertical stack + ffmpeg stitch.
    GATE_MOVIE=pair3087_lig3oh ~/miniforge3/bin/python pymol_renders/gate_latch_trace_panel.py

Styling, drift removal (intra_fit on core CA, excluding gate/latch/Lb7a5/
C-term), and the rotate-the-model-not-the-camera orientation (gate upper-
left, latch upper-right) are carried over from
gate_latch_pocket_oriented_movie.py; see its docstring for rationale. The
orientation basis is computed from state 1, which here is the open state.

Time and the ghost legend are drawn in the trace strip, not in-scene.

Ghost: gate (84-90) and latch (114-118) from the window's most-open state
(see pick_open_state()),
copied after drift-fit + reorientation into a single-state object, drawn
as an opaque light-gray tube (see GHOST_HEX for why not translucent). With
PyMOL's default static_singletons=1 a single-state object is drawn in every
frame, so it stays fixed while the live loops pull away from it.
"""

import os
import sys

from pymol import cmd, util

# Hard-coded rather than derived from __file__: under `pymol -cq script.py`
# __file__ doesn't point at this script, so the sibling import fails.
sys.path.insert(0, "/Users/ivanatang/Developer/biosensors/pymol_renders")
from gate_closure_movie_config import REPO_ROOT, TIME_SPACING_NS, get_movie  # noqa: E402

MOVIE = get_movie()
OUT_DIR = os.path.join(REPO_ROOT, "pymol_renders", "output")
PDB = MOVIE["pdb"]
XTC = MOVIE["xtc"]
RAW_N_FRAMES_TOTAL = MOVIE["n_raw_frames"]
T0_NS = MOVIE["t0_ns"]

PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
LIGAND_RESN = "LIG"
GATE_RESI = "84-90"
LATCH_RESI = "114-118"
FLEXIBLE_EXCLUDE_RESI = "84-90+114-118+148-166"
DIST_GATE_RESI = "88"

GATE_COLOR_HEX = "#FE6100"
LATCH_COLOR_HEX = "#DC267F"
BASE_PROTEIN_HEX = "#5B7FBE"
LIGAND_CARBON_HEX = "#00E5FF"
SURFACE_HEX = "#B9CBE8"
BG_HEX = "#FFFFFF"
# Opaque light gray, not cartoon_transparency: a translucent ghost behind the
# translucent surface gets dropped by transparency_mode 2's sorting (tested:
# visible with the protein disabled, invisible with it shown). Thinner than
# the live loops (0.4) so the colored loops win wherever the two overlap.
GHOST_HEX = "#C4C4C4"
GHOST_TUBE_RADIUS = 0.28

DIST_DASH_COLOR = "gray20"
DIST_OBJ_NAME = "gate_latch_dist"
DIST_LABEL_OBJ_NAME = "gate_latch_dist_label"
DIST_LABEL_Y_OFFSET = 1.2
DIST_LABEL_SIZE = 24

SURFACE_TRANSPARENCY = 0.6
# Tighter than gate_latch_pocket_oriented_movie.py's 6.0: that buffer was
# headroom for an in-scene time label, which here lives in the trace strip
# (gate_latch_trace_panel.py). Zoom target is the loops alone (see main()),
# so this buffer is what keeps the ligand's ring end and pocket in frame.
ZOOM_BUFFER = 5.0
CAMERA_Y_SHIFT = 7.0  # Angstroms; see main()

# 1600x1000 so the 1600x400 trace strip stacks underneath into 1600x1400.
VIEWPORT_W = 1600
VIEWPORT_H = 1000

STRIDE = int(os.environ.get("GATE_LATCH_STRIDE", "2"))
RAW_STOP = RAW_N_FRAMES_TOTAL
# load_traj(start=1, interval=N) skips raw frame 1 and keeps frames N, 2N, ...
# (its own log prints "skipping set 1 ... read set 2 into state 1"), so state
# s is raw frame s*STRIDE (see gate_closure_movie_config.py).
N_OUTPUT_FRAMES = RAW_N_FRAMES_TOTAL // STRIDE
FRAME_START = int(os.environ.get("GATE_LATCH_FRAME_START", "1"))
FRAME_STOP = int(os.environ.get("GATE_LATCH_FRAME_STOP", str(N_OUTPUT_FRAMES)))

OUT_NAME = MOVIE["out_name"]
FRAME_DIR = os.path.join(OUT_DIR, f".{OUT_NAME}_frames3d")


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


def _add(a, b, scale=1.0):
    """Returns a + b * scale for 3-vectors a, b."""
    return tuple(a[i] + b[i] * scale for i in range(3))


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
    length = _dot(v, v) ** 0.5
    if length < 1e-6:
        return (0.0, 0.0, 1.0)
    return tuple(c / length for c in v)


def _centroid(model):
    """Returns the mean atomic coordinate of a PyMOL chempy model.

    Args:
        model: PyMOL chempy model (has a .atom list).

    Returns:
        tuple[float, float, float]: (x, y, z) centroid.
    """
    n = len(model.atom)
    return tuple(sum(a.coord[i] for a in model.atom) / n for i in range(3))


def setup_scene():
    """Resets PyMOL, loads topology + trajectory, and defines selections.

    Returns:
        int: Number of trajectory states loaded.

    Raises:
        RuntimeError: A required selection is empty or the Ca88/partner
            selections don't resolve to exactly one atom each.
    """
    cmd.reinitialize()
    cmd.set_color("bg_flat", list(hex_to_rgb01(BG_HEX)))
    cmd.bg_color("bg_flat")
    cmd.set("ray_opaque_background", 1)
    cmd.set("orthoscopic", 0)
    cmd.set("field_of_view", 28)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("ambient", 0.32)
    cmd.set("direct", 0.65)
    cmd.set("specular", 0.3)
    cmd.set("light_count", 3)
    cmd.set("two_sided_lighting", 1)
    cmd.set("transparency_mode", 2)
    cmd.set("surface_quality", 1)
    cmd.set("defer_builds_mode", 3)

    cmd.load(PDB, "mol")
    cmd.load_traj(XTC, "mol", state=1, start=1, stop=RAW_STOP, interval=STRIDE, format="xtc")
    cmd.remove("mol and hydro")
    n_states = cmd.count_states("mol")
    print(f"[gate_latch_closure_ghost_movie] {MOVIE['key']}: loaded {n_states} states (stride={STRIDE})")

    protein = f"mol and chain {PROTEIN_CHAIN} and polymer"
    cmd.select("protein_sel", protein)
    cmd.select("gate_sel", f"{protein} and resi {GATE_RESI}")
    cmd.select("latch_sel", f"{protein} and resi {LATCH_RESI}")
    cmd.select("ligand_sel", f"mol and chain {LIGAND_CHAIN} and resn {LIGAND_RESN}")
    cmd.select("landmark_sel", "gate_sel or latch_sel or ligand_sel")
    cmd.select("core_fit_sel", f"{protein} and name CA and not resi {FLEXIBLE_EXCLUDE_RESI}")
    cmd.select("gate_ca88_sel", f"{protein} and resi {DIST_GATE_RESI} and name CA")
    cmd.select("partner_sel", MOVIE["partner_sel"])

    for sel in ("core_fit_sel", "gate_sel", "latch_sel", "ligand_sel"):
        if cmd.count_atoms(sel) == 0:
            raise RuntimeError(f"[gate_latch_closure_ghost_movie] {sel} is empty")
    for sel in ("gate_ca88_sel", "partner_sel"):
        if cmd.count_atoms(sel) != 1:
            raise RuntimeError(f"[gate_latch_closure_ghost_movie] {sel} is not exactly one atom")
    return n_states


def style_scene():
    """Applies the cartoon+surface styling of gate_latch_pocket_oriented_movie.py."""
    cmd.hide("everything", "mol")
    cmd.set_color("base_color", list(hex_to_rgb01(BASE_PROTEIN_HEX)))
    cmd.set_color("gate_color", list(hex_to_rgb01(GATE_COLOR_HEX)))
    cmd.set_color("latch_color", list(hex_to_rgb01(LATCH_COLOR_HEX)))
    cmd.color("base_color", "protein_sel")
    cmd.color("gate_color", "gate_sel")
    cmd.color("latch_color", "latch_sel")

    cmd.show("cartoon", "protein_sel")
    cmd.set("cartoon_loop_radius", 0.4, "gate_sel or latch_sel")
    cmd.set("cartoon_tube_radius", 0.4, "gate_sel or latch_sel")

    cmd.show("surface", "protein_sel")
    cmd.set("transparency", SURFACE_TRANSPARENCY, "protein_sel")
    cmd.set_color("surface_neutral", list(hex_to_rgb01(SURFACE_HEX)))
    cmd.set("surface_color", "surface_neutral", "protein_sel")

    cmd.show("sticks", "(gate_sel or latch_sel) and not name C+N+O")
    cmd.set("stick_radius", 0.18, "gate_sel or latch_sel")
    cmd.color("gate_color", "gate_sel and elem C")
    cmd.color("latch_color", "latch_sel and elem C")

    cmd.set_color("lig_carbon_color", list(hex_to_rgb01(LIGAND_CARBON_HEX)))
    cmd.show("sticks", "ligand_sel")
    cmd.show("spheres", "ligand_sel")
    cmd.set("stick_radius", 0.3, "ligand_sel")
    cmd.set("sphere_scale", 0.26, "ligand_sel")
    cmd.color("lig_carbon_color", "ligand_sel")
    util.cnc("ligand_sel")
    cmd.deselect()


def orient_gate_left_latch_right():
    """Rotates every state so gate reads upper-left and latch upper-right.

    Same basis and model-rotation workaround as
    gate_latch_pocket_oriented_movie.py's function of the same name.
    """
    gate_c = _centroid(cmd.get_model("gate_sel and name CA", state=1))
    latch_c = _centroid(cmd.get_model("latch_sel and name CA", state=1))
    core_c = _centroid(cmd.get_model("core_fit_sel", state=1))
    pocket_c = _centroid(cmd.get_model("landmark_sel", state=1))

    gate_dir = _norm(_sub(gate_c, pocket_c))
    latch_dir = _norm(_sub(latch_c, pocket_c))
    up_axis = _norm(_add(gate_dir, latch_dir))
    if _dot(up_axis, _sub(pocket_c, core_c)) < 0:
        up_axis = tuple(-c for c in up_axis)
    right_raw = _sub(latch_dir, gate_dir)
    right_axis = _norm(_add(right_raw, up_axis, -_dot(right_raw, up_axis)))
    forward_axis = _cross(right_axis, up_axis)

    rot = list(right_axis) + [0] + list(up_axis) + [0] + list(forward_axis) + [0, 0, 0, 0, 1]
    cmd.transform_selection("mol", rot, homogenous=1, state=0)


def pick_open_state(n_states: int) -> int:
    """Returns the state with the largest ~1 ns rolling-mean gate-to-partner distance.

    The window's first frame isn't necessarily its most open one (pair_3087
    starts at ~0.95 nm but peaks at ~1.3 nm just before closing), and the
    ghost should show the widest opening. A rolling mean, not the raw max,
    so a one-frame spike can't pick the ghost.

    Args:
        n_states: Number of loaded trajectory states.

    Returns:
        1-indexed state number.
    """
    d = [cmd.get_distance("gate_ca88_sel", "partner_sel", state=s) for s in range(1, n_states + 1)]
    half = max(1, int(round(0.5 / (STRIDE * TIME_SPACING_NS))))
    best_s, best_mean = 1, -1.0
    for i in range(half, len(d) - half):
        mean = sum(d[i - half:i + half + 1]) / (2 * half + 1)
        if mean > best_mean:
            best_s, best_mean = i + 1, mean
    return best_s


def add_open_state_ghost(open_state: int):
    """Creates a fixed, opaque light-gray copy of the gate and latch at open_state.

    Must run after drift-fit and orientation so the ghost sits in the same
    frame as the live loops. Backbone-only cartoon (no side chains) so it
    reads as an outline of where the loops were, not a second structure.
    Writes the ghost's time to FRAME_DIR/ghost_time_ns.txt for the trace
    strip's legend (gate_latch_trace_panel.py).

    Args:
        open_state: 1-indexed state to copy (see pick_open_state()).
    """
    cmd.create("ghost_open", "gate_sel or latch_sel", source_state=open_state, target_state=1)
    cmd.hide("everything", "ghost_open")
    cmd.set_color("ghost_color", list(hex_to_rgb01(GHOST_HEX)))
    cmd.color("ghost_color", "ghost_open")
    cmd.show("cartoon", "ghost_open")
    cmd.set("cartoon_loop_radius", GHOST_TUBE_RADIUS, "ghost_open")
    cmd.set("cartoon_tube_radius", GHOST_TUBE_RADIUS, "ghost_open")

    ghost_ns = T0_NS + (open_state * STRIDE - 1) * TIME_SPACING_NS
    os.makedirs(FRAME_DIR, exist_ok=True)
    with open(os.path.join(FRAME_DIR, "ghost_time_ns.txt"), "w") as f:
        f.write(f"{ghost_ns:.4f}\n")
    print(f"[gate_latch_closure_ghost_movie] ghost = state {open_state} (t={ghost_ns:.2f} ns)")


def render_frames(frame_start, frame_stop):
    """Renders OpenGL PNGs for states [frame_start, frame_stop], 1-indexed inclusive.

    Redraws the Ca88-to-partner dash and its midpoint label from scratch each
    frame (see gate_latch_pocket_oriented_movie.py render_frames()). The
running time is drawn by gate_latch_trace_panel.py, not here.

    Args:
        frame_start: First state to render.
        frame_stop: Last state to render.
    """
    cmd.viewport(VIEWPORT_W, VIEWPORT_H)
    for s in range(frame_start, frame_stop + 1):
        cmd.frame(s)
        sim_ns = T0_NS + (s * STRIDE - 1) * TIME_SPACING_NS
        g = cmd.get_atom_coords("gate_ca88_sel", state=s)
        l = cmd.get_atom_coords("partner_sel", state=s)
        dist_nm = _dot(_sub(g, l), _sub(g, l)) ** 0.5 / 10.0

        cmd.delete(DIST_OBJ_NAME)
        cmd.distance(DIST_OBJ_NAME, "gate_ca88_sel", "partner_sel", state=s)
        cmd.hide("labels", DIST_OBJ_NAME)
        cmd.color(DIST_DASH_COLOR, DIST_OBJ_NAME)
        cmd.set("dash_width", 5, DIST_OBJ_NAME)
        cmd.set("dash_radius", 0.08, DIST_OBJ_NAME)
        cmd.set("dash_gap", 0.15, DIST_OBJ_NAME)
        cmd.set("dash_length", 0.15, DIST_OBJ_NAME)

        mid = tuple((g[i] + l[i]) / 2.0 for i in range(3))
        cmd.delete(DIST_LABEL_OBJ_NAME)
        cmd.pseudoatom(DIST_LABEL_OBJ_NAME, pos=[mid[0], mid[1] + DIST_LABEL_Y_OFFSET, mid[2]])
        cmd.hide("everything", DIST_LABEL_OBJ_NAME)
        cmd.set("label_size", DIST_LABEL_SIZE, DIST_LABEL_OBJ_NAME)
        cmd.set("label_color", "black", DIST_LABEL_OBJ_NAME)
        cmd.set("label_font_id", 7, DIST_LABEL_OBJ_NAME)
        cmd.set("label_outline_color", "white", DIST_LABEL_OBJ_NAME)
        cmd.label(DIST_LABEL_OBJ_NAME, '"%.2f nm"' % dist_nm)

        cmd.png(os.path.join(FRAME_DIR, f"frame_{s - 1:05d}.png"), ray=0, quiet=1)
        if s == frame_start or s % 25 == 0 or s == frame_stop:
            print(f"[gate_latch_closure_ghost_movie] frame {s} t={sim_ns:.2f} ns Ca88-{MOVIE['partner_label']}={dist_nm:.3f} nm")


def main():
    """Loads, fits, styles, orients, adds the ghost, and renders this chunk's frames."""
    n_states = setup_scene()
    rms = cmd.intra_fit("core_fit_sel", state=1)
    print(f"[gate_latch_closure_ghost_movie] intra_fit max RMSD to state 1 = {max(rms):.3f} A")
    style_scene()
    orient_gate_left_latch_right()
    add_open_state_ghost(pick_open_state(n_states))
    cmd.frame(1)
    # Center on the two loops, not landmark_sel: the ligand's tail reaches
    # well below the pocket, so centering on gate+latch+ligand left the top
    # third of the frame empty. The ligand's ring end stays in view below.
    cmd.zoom("gate_sel or latch_sel", buffer=ZOOM_BUFFER, state=1)
    # The loops are the top edge of the protein, so a loop-centered view is
    # half empty above them; slide the view so they sit in the upper part of
    # the frame and the pocket/ligand fill the rest. Model is pre-rotated
    # (world +Y = up), so a plain y move is enough.
    cmd.move("y", CAMERA_Y_SHIFT)

    os.makedirs(FRAME_DIR, exist_ok=True)
    render_frames(FRAME_START, min(FRAME_STOP, n_states))


main()
