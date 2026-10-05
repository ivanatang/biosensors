#!/usr/bin/env python
"""Renders a zoomed-into-the-pocket movie of gate/latch closure for pair_3101_binder,
oriented so the gate loop sits upper-left and the latch loop sits upper-right.

Same source data, drift-removal, and pyr1_lca_biosensor.png-matching palette
as pymol_renders/gate_latch_movie_hero.py (pair_3101_binder's t=0-30 ns
closure window; blue protein, orange gate, pink latch, cyan ball-and-stick
ligand, white background -- see that script's docstring for the full
mechanistic/data-selection rationale, not repeated here). This is a new,
independent file because the camera composition is a deliberate departure
from gate_latch_movie_hero.py's pocket-opening-axis view: by request, this
script derives its own camera so the gate loop (resi 84-90) lands in the
screen's upper-left and the latch loop (resi 114-118) lands in the upper-
right, with the ligand roughly centered between them, rather than reusing
the whole-protein pocket-opening derivation.

Also by request, this style keeps a translucent molecular surface (dropped
in gate_latch_movie_hero.py's v3 for its tighter landmark-only zoom -- see
that script's docstring -- because at THAT camera's close range the surface
filled the frame like fog). This script's camera looks straight down the
pocket opening (see "Camera" below), so the surface reads as a clean
"porthole" around the ligand instead of a fog wall -- confirmed against a
real preview render before committing to a full movie run.

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/gate_latch_pocket_oriented_movie.py

Smoke test (renders only a handful of frames, fast, writes to a
"_smoketest" suffixed mp4 so it never clobbers the real deliverable):

    GATE_LATCH_MAX_FRAMES=25 /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/gate_latch_pocket_oriented_movie.py

Other environment variable overrides (all optional, same meaning as
gate_latch_movie_hero.py):
    GATE_LATCH_STRIDE (default 2), GATE_LATCH_MAX_FRAMES,
    GATE_LATCH_FPS (default 10), GATE_LATCH_KEEP_FRAMES ("1" to keep temp
    frames), FFMPEG_BIN.

Camera -- gate upper-left, latch upper-right. Built from state-1 (post
drift-fit) CA centroids of the gate and latch selections, relative to a
pocket center (centroid of gate+latch+ligand):

    gate_dir  = normalize(gate_centroid  - pocket_center)
    latch_dir = normalize(latch_centroid - pocket_center)
    up_axis    = normalize(gate_dir + latch_dir)   # both loops read "upper"
    right_axis = normalize(orthogonalize(latch_dir - gate_dir, up_axis))
                                                     # latch reads "right",
                                                     # gate reads "left"
    forward_axis = cross(right_axis, up_axis)       # completes a
                                                     # right-handed frame

up_axis is flipped if it points toward the stable core rather than away
from it (pocket_center - core_centroid), so gate/latch land "upper" rather
than "lower" regardless of which way that sum vector happens to point.

Implementation note -- rotate the model, not the camera. The natural way
to apply this basis is cmd.set_view() with (right_axis, up_axis,
forward_axis) written into the view matrix (the approach
pyr1_lca_hero_complex.py and gate_latch_movie_hero.py both use for their
own, differently-derived, camera bases). That approach was tried here first
and produced visibly wrong, seemingly rotated/mirrored compositions --
isolated tests (placing simple marker points at known offsets along a
candidate basis with no other geometry) showed cmd.get_view() faithfully
echoing back the exact matrix passed to cmd.set_view(), yet the actual
OpenGL render not matching that matrix, for generic (non axis-aligned)
bases specifically -- simple axis-aligned bases rendered correctly. This
points to a matrix-to-internal-rotation conversion issue in this PyMOL
build for generic view matrices, not a bug in the basis derivation itself
(confirmed correct via dot-product sign checks against gate/latch
positions). The workaround, verified against the same isolated marker
test: physically rotate the model's coordinates (cmd.transform_selection())
so that (right_axis, up_axis, forward_axis) become the new (X, Y, Z) axes,
then render with PyMOL's default/identity camera, which was confirmed
reliable. The rotation is applied once, across every loaded trajectory
state (state=0), so it composes cleanly with the existing rigid-body-drift
fit: drift-fit first (removes whole-protein tumbling relative to state 1,
see design decision 1 below), then this one fixed reorientation (state=0,
applied identically to all states), so the fixed default camera only ever
shows internal conformational change, exactly like the other movies in
this repo -- just via a rotated *model* instead of a rotated *camera*.

Design decisions carried over unchanged from gate_latch_movie_hero.py (see
that script's own docstring for the full rationale):

  1. Rigid-body-drift removal via cmd.intra_fit() on core_fit_sel (protein
     CA atoms excluding gate, latch, and the Lb7a5/C-terminal recoil loop),
     fit to state 1, before any camera work.
  2. Fixed camera / no per-frame orient-zoom-turn calls -- only cmd.frame()
     changes between rendered frames.
  3. On-screen running time label via a coordinate-anchored pseudoatom
     (re-derived below for this script's own default-camera setup, since
     there is no longer a camera "up vector" to place it along -- the
     rotated coordinate frame's up is simply world +Y).
  4. OpenGL (non-ray-traced) rendering for speed.

Distance measurement (added by request, v2/v3): a dashed line plus a live
"Ca88-Ca116: X.XX nm" readout in the HUD text.

v2 first tried Ca88 (gate) to LCA's ring hydroxyl oxygen (the "OH group at
the top of LCA" from the original request, identified by coordinate/bonding
analysis as distinct from the ligand's C24 carboxylic-acid oxygens, and
confirmed to sit near the gate rather than the solvent-facing tail).
Checked against the real trajectory before committing to a full render,
that pair does NOT show a clean closing trend -- it fluctuates ~3.3-8.0 A
throughout, with windowed means of 5.19 A (t=0-2ns) -> 4.54 A (t=2-4ns) ->
4.88 A (t=4-30ns), i.e. a dip during the transition window that doesn't
stay closed afterward, unlike the metric this project already uses
elsewhere for "the gate closes" claims.

v3 (current) switches to that already-established metric instead, per user
request after being shown the v2 finding: Ca of gate residue 88 ("P88") to
Ca of latch residue 116, as built by run_gate_latch.sh (`GATE_RES=88` /
`LATCH_RES=116`, index gate_latch116.ndx) and archived per sequence as
gate_latch116_timeseries.xvg. Don't confuse it with gate_latch_timeseries.xvg
(index gate_latch.ndx), an older Ca88-Ca117 measurement: for pair_3101 its
atom pair is 1383/1835, and 1835 is Ca of LEU117 in prod_md_500ns.gro (the
116 index is 1383/1811). An earlier version of this docstring mixed the two
up and quoted the Ca117 file's ~1.02 nm at t=0.

pair_3101_binder's archived gate_latch116_timeseries.xvg: mean 0.71 nm over
t=0-0.5 ns, 0.59 nm over 2-4 ns, 0.55 nm (range 0.47-0.72) over 4-30 ns.
The movie computes Ca88-Ca116 fresh every frame from its own topology/
trajectory (DIST_GATE_RESI="88", DIST_LATCH_RESI="116", plain resi+name-CA
selections; no atom-`id` hack needed since Ca atoms have real per-residue
names, unlike this trajectory's anonymized LIG atoms) and gets the same
values (0.70 nm at t=0 -> ~0.53-0.62 nm from ~4 ns on), so what's drawn and
what's printed agree with the archived file. Note the change is small
(~1.5 A), and the gate reopens to ~0.8-1.2 nm after ~30 ns (mean 0.81 nm
over 30-60 ns), outside this movie's window; see
gate_latch_closure_ghost_movie.py for windows with larger closures.

The dashed line itself is recreated every frame (delete + cmd.distance(...,
state=s) inside render_frames(), not created once with state=0) for the
same "don't rely on undocumented per-version PyMOL auto-update behavior"
caution as this script's other camera/geometry choices (see "Implementation
note" above) -- explicit per-frame recomputation is simple and guaranteed
correct regardless. Dash styling (color gray20, hidden per-object numeric
label, dash_width/radius/gap/length) matches
pymol_renders/gate_latch_water_network.py's established convention for
in-scene distance dashes in this repo -- the dash object's own numeric
label is always hidden.

v4 (current): two fixes after watching the v3 render. (1) Both the time
text and the appended "Ca88-Ca116: X.XX nm" text were clipped at the top
frame edge -- ZOOM_BUFFER (3.5) left too little headroom above
landmark_sel for the label's actual pixel height once rendered for real,
not just for its anchor point. Bumped to 6.0, paired with a larger anchor
margin (add_time_label_anchor()). (2) By request, the distance number moved
out of the fixed top HUD into its own label (DIST_LABEL_OBJ_NAME),
recreated every frame at the dashed line's own current midpoint (offset by
DIST_LABEL_Y_OFFSET along the rotated model's +Y so it clears the line
itself) -- so the number visibly sits next to, and shortens along with, the
line it's measuring, rather than being read off a static corner counter.
time_label_anchor reverted to "N ns" only.
"""

import glob
import os
import shutil
import subprocess

from pymol import cmd, util

# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
REPO_ROOT = "/Users/ivanatang/Developer/biosensors"
OUT_DIR = os.path.join(REPO_ROOT, "pymol_renders", "output")

DATA_DIR = (
    "/Users/ivanatang/Library/CloudStorage/OneDrive-UCB-O365/Shirts Lab/"
    "LCA_boltz_models/binders/pair_3101_binder/"
    "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
)
PDB = os.path.join(DATA_DIR, "medoid_PL.pdb")
XTC = os.path.join(
    REPO_ROOT,
    "pymol_renders",
    "scratch_traj",
    "gate_closure_pair3101_0_30ns_PL.xtc",
)

# 801 raw frames at 37.5 ps spacing, spanning simulation time 0-30 ns.
RAW_N_FRAMES_TOTAL = 801
TIME_SPACING_NS = 0.0375
T0_NS = 0.0

PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
LIGAND_RESN = "LIG"

GATE_RESI = "84-90"
LATCH_RESI = "114-118"
FLEXIBLE_EXCLUDE_RESI = "84-90+114-118+148-166"

# Gate-closure distance measurement (by request): live Ca88 (gate) to Ca116
# (latch) distance, so the closing motion has a number attached to it, not
# just a visual impression. This is this project's already-established
# gate-latch distance metric (see module docstring "Distance measurement"
# section and run_gate_latch.sh's GATE_RES=88/LATCH_RES=116), not a
# from-scratch choice.
DIST_GATE_RESI = "88"           # single residue within the gate loop
                                 # (84-90)
DIST_LATCH_RESI = "116"         # single residue within the latch loop
                                 # (114-118)

# Colors -- match pyr1_lca_biosensor.png / gate_latch_movie_hero.py
GATE_COLOR_HEX = "#FE6100"     # orange
LATCH_COLOR_HEX = "#DC267F"    # pink
BASE_PROTEIN_HEX = "#5B7FBE"   # medium blue
LIGAND_CARBON_HEX = "#00E5FF"  # bright cyan -- distinct from protein/gate/latch
SURFACE_HEX = "#B9CBE8"        # pale glassy blue, matches pyr1_lca_biosensor.png
BG_HEX = "#FFFFFF"             # white, per request

# Gate-closure distance dash styling -- matches
# pymol_renders/gate_latch_water_network.py's established convention for
# in-scene distance dashes (gray20, numeric label hidden in favor of an
# externally-placed value -- see module docstring "Distance measurement").
DIST_DASH_COLOR = "gray20"
DIST_DASH_WIDTH = 5
DIST_DASH_RADIUS = 0.08
DIST_DASH_GAP = 0.15
DIST_DASH_LENGTH = 0.15
DIST_OBJ_NAME = "gate_latch_dist"

# Distance number placement -- by request, parked next to the dashed line's
# own (per-frame) midpoint rather than folded into the fixed top HUD text,
# so it visibly moves and updates alongside the line (see render_frames()).
# Smaller than the time label (DIST_LABEL_SIZE < the time anchor's 30) since
# it's a close-up, in-scene annotation rather than a page-level HUD.
DIST_LABEL_OBJ_NAME = "gate_latch_dist_label"
DIST_LABEL_Y_OFFSET = 1.2  # Angstroms, lifts the label clear of the dash
DIST_LABEL_SIZE = 24

CARTOON_TRANSPARENCY = 0.0
SURFACE_TRANSPARENCY = 0.6

ZOOM_BUFFER = 6.0  # tight, pocket-focused zoom (see request: "zoomed into
                 # the binding pocket"); originally tuned to 3.5 (buffer=7
                 # and buffer=4 left too much surrounding fold in frame,
                 # buffer=2.5 clipped the time label at the top edge), but
                 # 3.5 turned out to still clip the label's top edge once
                 # rendered for real (the zoom target is landmark_sel alone,
                 # not the label anchor -- see set_fixed_camera() -- so the
                 # anchor's own margin has to stay safely inside this
                 # buffer, which 3.5 didn't leave enough room for). Bumped
                 # to 6.0, paired with a larger anchor margin (see
                 # add_time_label_anchor()), to leave real headroom for the
                 # label's pixel height rather than just its anchor point.

LIGAND_STICK_RADIUS = 0.3
LIGAND_SPHERE_SCALE = 0.26

VIEWPORT_W = 1600
VIEWPORT_H = 1200

FFMPEG_BIN = os.environ.get("FFMPEG_BIN", "/opt/homebrew/bin/ffmpeg")

# --------------------------------------------------------------------------
# Env-var-driven run parameters
# --------------------------------------------------------------------------
STRIDE = int(os.environ.get("GATE_LATCH_STRIDE", "2"))
FPS = int(os.environ.get("GATE_LATCH_FPS", "10"))
KEEP_FRAMES = os.environ.get("GATE_LATCH_KEEP_FRAMES", "0") == "1"
_max_frames_env = os.environ.get("GATE_LATCH_MAX_FRAMES", "").strip()

FULL_N_OUTPUT_FRAMES = (RAW_N_FRAMES_TOTAL - 1) // STRIDE + 1

if _max_frames_env:
    IS_SMOKE_TEST = True
    N_OUTPUT_FRAMES = min(int(_max_frames_env), FULL_N_OUTPUT_FRAMES)
else:
    IS_SMOKE_TEST = False
    N_OUTPUT_FRAMES = FULL_N_OUTPUT_FRAMES

RAW_STOP = min(RAW_N_FRAMES_TOTAL, (N_OUTPUT_FRAMES - 1) * STRIDE + 1)

OUT_NAME = "gate_latch_pocket_oriented_pair3101_binder_0-30ns"
if IS_SMOKE_TEST:
    OUT_NAME += "_smoketest"
OUT_MP4 = os.path.join(OUT_DIR, OUT_NAME + ".mp4")

# Chunked rendering (added after a full single-process run got OOM-killed on
# this 8 GB machine under heavy load from other running apps -- with
# defer_builds_mode=3 the cartoon/surface geometry for 400 states isn't all
# resident at once, but the process's cumulative footprint across ~400
# frames of surface+cartoon+distance-object churn was still enough to get
# killed partway through (138/400 frames rendered) when free memory was
# down to ~170 MB. Splitting the render into several independent PyMOL
# process invocations -- each exits when its chunk is done, so the OS fully
# reclaims its memory before the next chunk starts -- sidesteps this
# without touching render quality. Frames from every chunk land in one
# fixed, non-temp directory (FRAME_DIR, not tempfile.mkdtemp() -- a fresh
# temp dir per invocation would scatter each chunk's frames into a
# different, unstitchable directory) keyed by GLOBAL output-frame index, so
# a final GATE_LATCH_STITCH_ONLY=1 invocation (skips loading/rendering
# entirely) can ffmpeg-stitch the complete set in one pass. Default
# behavior (no chunk env vars set) is unchanged: one process renders every
# frame and stitches, exactly as before.
#
# Usage for a chunked run (bash), 3 chunks of ~134 frames each:
#   for start in 1 135 269; do
#     stop=$((start + 133)); [ $stop -gt 400 ] && stop=400
#     GATE_LATCH_FRAME_START=$start GATE_LATCH_FRAME_STOP=$stop GATE_LATCH_SKIP_STITCH=1 \
#       /opt/homebrew/bin/pymol -cq gate_latch_pocket_oriented_movie.py
#   done
#   GATE_LATCH_STITCH_ONLY=1 /opt/homebrew/bin/pymol -cq gate_latch_pocket_oriented_movie.py
FRAME_START = int(os.environ.get("GATE_LATCH_FRAME_START", "1"))
FRAME_STOP = int(os.environ.get("GATE_LATCH_FRAME_STOP", str(N_OUTPUT_FRAMES)))
SKIP_STITCH = os.environ.get("GATE_LATCH_SKIP_STITCH", "0") == "1"
STITCH_ONLY = os.environ.get("GATE_LATCH_STITCH_ONLY", "0") == "1"
FRAME_DIR = os.path.join(OUT_DIR, f".{OUT_NAME}_frames")


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
    length = (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) ** 0.5
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
    # Many-state memory/perf: build cartoon/surface geometry per-frame
    # instead of caching every state's geometry simultaneously.
    cmd.set("defer_builds_mode", 3)


def load_topology_and_trajectory():
    """Loads the topology PDB and trajectory, strips hydrogens.

    Returns:
        int: Number of trajectory states loaded.
    """
    cmd.load(PDB, "mol")
    n_atoms_pdb = cmd.count_atoms("mol")
    print(f"[gate_latch_pocket_oriented_movie] loaded topology: {n_atoms_pdb} atoms from {PDB}")

    cmd.load_traj(
        XTC,
        "mol",
        state=1,
        start=1,
        stop=RAW_STOP,
        interval=STRIDE,
        format="xtc",
    )
    n_states = cmd.count_states("mol")
    print(
        f"[gate_latch_pocket_oriented_movie] loaded {n_states} states from trajectory "
        f"(stride={STRIDE}, raw stop frame={RAW_STOP} of {RAW_N_FRAMES_TOTAL})"
    )

    cmd.remove("mol and hydro")
    print(f"[gate_latch_pocket_oriented_movie] after removing hydrogens: {cmd.count_atoms('mol')} atoms")

    return n_states


def define_selections():
    """Creates and validates the protein/gate/latch/ligand/core-fit selections.

    Raises:
        RuntimeError: Any of the required selections is empty.
    """
    protein_sel = f"mol and chain {PROTEIN_CHAIN} and polymer"
    gate_sel = f"{protein_sel} and resi {GATE_RESI}"
    latch_sel = f"{protein_sel} and resi {LATCH_RESI}"
    ligand_sel = f"mol and chain {LIGAND_CHAIN} and resn {LIGAND_RESN}"
    core_fit_sel = f"{protein_sel} and name CA and not resi {FLEXIBLE_EXCLUDE_RESI}"

    gate_ca88_sel = f"{protein_sel} and resi {DIST_GATE_RESI} and name CA"
    latch_ca116_sel = f"{protein_sel} and resi {DIST_LATCH_RESI} and name CA"

    cmd.select("protein_sel", protein_sel)
    cmd.select("gate_sel", gate_sel)
    cmd.select("latch_sel", latch_sel)
    cmd.select("ligand_sel", ligand_sel)
    cmd.select("landmark_sel", "gate_sel or latch_sel or ligand_sel")
    cmd.select("core_fit_sel", core_fit_sel)
    cmd.select("gate_ca88_sel", gate_ca88_sel)
    cmd.select("latch_ca116_sel", latch_ca116_sel)

    n_core = cmd.count_atoms("core_fit_sel")
    n_gate = cmd.count_atoms("gate_sel")
    n_latch = cmd.count_atoms("latch_sel")
    n_lig = cmd.count_atoms("ligand_sel")
    n_gate_ca88 = cmd.count_atoms("gate_ca88_sel")
    n_latch_ca116 = cmd.count_atoms("latch_ca116_sel")
    print(
        f"[gate_latch_pocket_oriented_movie] selections: core_fit={n_core} atoms, "
        f"gate={n_gate} atoms, latch={n_latch} atoms, ligand={n_lig} atoms, "
        f"gate_ca88={n_gate_ca88} atoms, latch_ca116={n_latch_ca116} atoms"
    )
    if n_core == 0 or n_gate == 0 or n_latch == 0 or n_lig == 0:
        raise RuntimeError(
            "[gate_latch_pocket_oriented_movie] one or more required selections is "
            "empty; check chain/resi/resn conventions before rendering."
        )
    if n_gate_ca88 != 1 or n_latch_ca116 != 1:
        raise RuntimeError(
            "[gate_latch_pocket_oriented_movie] gate_ca88_sel or latch_ca116_sel did "
            f"not resolve to exactly one atom (got {n_gate_ca88}, {n_latch_ca116}); "
            "check DIST_GATE_RESI/DIST_LATCH_RESI against this trajectory's topology."
        )


def remove_rigid_body_drift():
    """Fits every loaded state onto state 1 using stable-core CA atoms only.

    See module docstring design decision 1.
    """
    rms_list = cmd.intra_fit("core_fit_sel", state=1)
    print(
        f"[gate_latch_pocket_oriented_movie] intra_fit on core_fit_sel done; "
        f"max per-state RMSD to state 1 = {max(rms_list):.3f} A"
    )


def style_scene():
    """Applies cartoon+surface hybrid styling matching pyr1_lca_biosensor.png."""
    cmd.hide("everything", "mol")

    cmd.set_color("base_color", list(hex_to_rgb01(BASE_PROTEIN_HEX)))
    cmd.set_color("gate_color", list(hex_to_rgb01(GATE_COLOR_HEX)))
    cmd.set_color("latch_color", list(hex_to_rgb01(LATCH_COLOR_HEX)))
    cmd.color("base_color", "protein_sel")
    cmd.color("gate_color", "gate_sel")
    cmd.color("latch_color", "latch_sel")

    cmd.show("cartoon", "protein_sel")
    cmd.set("cartoon_transparency", CARTOON_TRANSPARENCY, "protein_sel")
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
    cmd.set("stick_radius", LIGAND_STICK_RADIUS, "ligand_sel")
    cmd.set("sphere_scale", LIGAND_SPHERE_SCALE, "ligand_sel")
    cmd.color("lig_carbon_color", "ligand_sel")
    util.cnc("ligand_sel")

    cmd.deselect()


def orient_gate_left_latch_right():
    """Rotates the model so gate reads upper-left and latch reads upper-right.

    See module docstring "Camera" and "Implementation note" for the basis
    derivation and why this rotates the *model* (cmd.transform_selection,
    state=0 -- every loaded trajectory state) rather than the camera.
    Computed once, from state-1 (already drift-fit) CA/ligand positions;
    the same fixed rotation is then applied identically to every state, so
    it composes with the drift-fit exactly like a fixed camera would.

    Returns:
        str: Name of a pseudoatom placed above the (now-rotated) landmark
        region, for the time-label HUD (see add_time_label_anchor()).
    """
    gate_model = cmd.get_model("gate_sel and name CA", state=1)
    latch_model = cmd.get_model("latch_sel and name CA", state=1)
    core_model = cmd.get_model("core_fit_sel", state=1)
    landmark_model = cmd.get_model("landmark_sel", state=1)

    gate_c = _centroid(gate_model)
    latch_c = _centroid(latch_model)
    core_c = _centroid(core_model)
    pocket_c = _centroid(landmark_model)

    gate_dir = _norm(_sub(gate_c, pocket_c))
    latch_dir = _norm(_sub(latch_c, pocket_c))

    up_axis = _norm(_add(gate_dir, latch_dir))
    if _dot(up_axis, _sub(pocket_c, core_c)) < 0:
        up_axis = tuple(-c for c in up_axis)

    right_raw = _sub(latch_dir, gate_dir)
    up_component = _dot(right_raw, up_axis)
    right_axis = _norm(_sub(right_raw, tuple(c * up_component for c in up_axis)))

    forward_axis = _cross(right_axis, up_axis)

    print(
        f"[gate_latch_pocket_oriented_movie] camera basis: right={right_axis}, "
        f"up={up_axis}, forward={forward_axis}"
    )

    rot = [
        right_axis[0], right_axis[1], right_axis[2], 0,
        up_axis[0], up_axis[1], up_axis[2], 0,
        forward_axis[0], forward_axis[1], forward_axis[2], 0,
        0, 0, 0, 1,
    ]
    cmd.transform_selection("mol", rot, homogenous=1, state=0)

    return add_time_label_anchor()


def add_time_label_anchor():
    """Places a fixed-in-model-space pseudoatom above the pocket, post-reorientation.

    Since orient_gate_left_latch_right() already rotated the model so "up"
    is world +Y, the anchor just needs to sit above the (rotated)
    landmark_sel region along +Y -- no camera-vector bookkeeping needed,
    unlike gate_latch_movie_hero.py's version of this function.

    Only carries the time text now ("N ns") -- the distance readout moved
    to a separate, per-frame-repositioned label next to the dashed line
    itself (see render_frames()), by request, so the number visibly tracks
    the line as it shortens. Margin bumped from the original +1.5 to +6.0:
    the +1.5 version got the top of the label glyphs clipped by the frame
    edge (confirmed by inspecting a rendered frame) -- ZOOM_BUFFER=3.5
    leaves very little headroom at this tight a zoom, and label text has
    real pixel height that a margin defined purely in landmark-radius terms
    doesn't account for.

    Returns:
        str: Name of the created pseudoatom object ("time_label_anchor").
    """
    anchor = "time_label_anchor"
    landmark_model = cmd.get_model("landmark_sel", state=1)
    center = _centroid(landmark_model)
    radius = max(
        (
            (a.coord[0] - center[0]) ** 2
            + (a.coord[1] - center[1]) ** 2
            + (a.coord[2] - center[2]) ** 2
        )
        ** 0.5
        for a in landmark_model.atom
    )
    label_pos = (center[0], center[1] + radius + 3.5, center[2])

    cmd.pseudoatom(anchor, pos=list(label_pos))
    cmd.hide("everything", anchor)
    cmd.set("label_size", 30, anchor)
    cmd.set("label_color", "black", anchor)
    cmd.set("label_font_id", 7, anchor)
    cmd.set("label_outline_color", "white", anchor)
    return anchor


def set_fixed_camera(anchor):
    """Zooms once on the reoriented landmark region; no further camera moves.

    Args:
        anchor (str): Name of the time-label pseudoatom, included in the
            zoom target so it stays on screen (see add_time_label_anchor()).
    """
    cmd.frame(1)
    # Zoom on landmark_sel alone (not the anchor) -- the anchor is a single
    # point well off to one side of the landmark region's mass, and
    # including it in the zoom target pulled the whole frame outward/
    # off-center to fit it, leaving too much empty space (found via
    # smoke-test render). It's placed close enough (see
    # add_time_label_anchor()) to fall inside this buffer anyway.
    cmd.zoom("landmark_sel", buffer=ZOOM_BUFFER)


def render_frames(n_states, frame_dir, frame_start, frame_stop):
    """Renders one OpenGL (non-ray-traced) PNG per trajectory state in [frame_start, frame_stop].

    Also redraws the gate-Ca88-to-latch-Ca116 distance dash every frame (see
    module docstring "Distance measurement") -- recreated from scratch each
    time (delete + cmd.distance(..., state=s)) rather than relying on any
    auto-updating behavior, so both the dash geometry and its numeric
    readout are always computed fresh for the frame actually being
    rendered. By request, the distance number now lives in its own label
    (DIST_LABEL_OBJ_NAME) parked just above the dashed line's own midpoint
    and recreated every frame at that frame's midpoint, rather than folded
    into the fixed top-of-frame time HUD -- so the number visibly sits next
    to, and moves with, the line it's measuring as the line shortens.
    time_label_anchor now only ever shows "N ns".

    frame_start/frame_stop (both 1-indexed, inclusive) let one process
    render only a sub-range of the loaded states -- see module docstring
    "Chunked rendering". PNG filenames use the global state index `s`
    (not a chunk-local counter), so multiple chunk invocations interleave
    correctly into one frame_dir for the final stitch.

    Args:
        n_states (int): Number of loaded trajectory states (all states are
            loaded regardless of frame_start/frame_stop; only rendering is
            restricted).
        frame_dir (str): Directory to write frame_00000.png etc. into.
        frame_start (int): First state (1-indexed) to render.
        frame_stop (int): Last state (1-indexed, inclusive) to render.
    """
    cmd.viewport(VIEWPORT_W, VIEWPORT_H)

    for s in range(frame_start, frame_stop + 1):
        cmd.frame(s)
        sim_ns = T0_NS + (s - 1) * STRIDE * TIME_SPACING_NS

        gate_c = cmd.get_atom_coords("gate_ca88_sel", state=s)
        latch_c = cmd.get_atom_coords("latch_ca116_sel", state=s)
        gate_latch_dist_nm = (
            sum((gate_c[i] - latch_c[i]) ** 2 for i in range(3)) ** 0.5
        ) / 10.0

        cmd.delete(DIST_OBJ_NAME)
        cmd.distance(DIST_OBJ_NAME, "gate_ca88_sel", "latch_ca116_sel", state=s)
        cmd.hide("labels", DIST_OBJ_NAME)
        cmd.color(DIST_DASH_COLOR, DIST_OBJ_NAME)
        cmd.set("dash_width", DIST_DASH_WIDTH, DIST_OBJ_NAME)
        cmd.set("dash_radius", DIST_DASH_RADIUS, DIST_OBJ_NAME)
        cmd.set("dash_gap", DIST_DASH_GAP, DIST_OBJ_NAME)
        cmd.set("dash_length", DIST_DASH_LENGTH, DIST_OBJ_NAME)

        # Distance number, repositioned every frame at the dash's own
        # midpoint (model is pre-rotated so "up" is always world +Y here --
        # see orient_gate_left_latch_right() -- so a plain +Y nudge lifts
        # the label clear of the line itself without any camera-vector
        # bookkeeping).
        mid = tuple((gate_c[i] + latch_c[i]) / 2.0 for i in range(3))
        dist_label_pos = (mid[0], mid[1] + DIST_LABEL_Y_OFFSET, mid[2])
        cmd.delete(DIST_LABEL_OBJ_NAME)
        cmd.pseudoatom(DIST_LABEL_OBJ_NAME, pos=list(dist_label_pos))
        cmd.hide("everything", DIST_LABEL_OBJ_NAME)
        cmd.set("label_size", DIST_LABEL_SIZE, DIST_LABEL_OBJ_NAME)
        cmd.set("label_color", "black", DIST_LABEL_OBJ_NAME)
        cmd.set("label_font_id", 7, DIST_LABEL_OBJ_NAME)
        cmd.set("label_outline_color", "white", DIST_LABEL_OBJ_NAME)
        cmd.label(DIST_LABEL_OBJ_NAME, '"%.2f nm"' % gate_latch_dist_nm)

        cmd.label("time_label_anchor", '"%d ns"' % int(round(sim_ns)))

        png_path = os.path.join(frame_dir, f"frame_{s - 1:05d}.png")
        cmd.png(png_path, ray=0, quiet=1)

        if s == frame_start or s % 25 == 0 or s == frame_stop:
            print(
                f"[gate_latch_pocket_oriented_movie] rendered frame {s}/{n_states} "
                f"(t={sim_ns:.1f} ns, Ca88-Ca116={gate_latch_dist_nm:.3f} nm)"
            )


def stitch_movie(frame_dir, n_states):
    """Stitches the rendered PNG frames into an MP4 with ffmpeg.

    Args:
        frame_dir (str): Directory of frame_%05d.png files (see render_frames).
        n_states (int): Number of frames, used to report movie duration.

    Raises:
        RuntimeError: ffmpeg exits non-zero; frames are left in `frame_dir`
            for debugging.
    """
    os.makedirs(OUT_DIR, exist_ok=True)
    ffmpeg_cmd = [
        FFMPEG_BIN,
        "-y",
        "-framerate", str(FPS),
        "-i", os.path.join(frame_dir, "frame_%05d.png"),
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-crf", "18",
        "-movflags", "+faststart",
        OUT_MP4,
    ]
    print("[gate_latch_pocket_oriented_movie] running:", " ".join(ffmpeg_cmd))
    result = subprocess.run(ffmpeg_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(
            f"[gate_latch_pocket_oriented_movie] ffmpeg failed (exit {result.returncode}); "
            f"frames left in {frame_dir} for debugging."
        )

    size_bytes = os.path.getsize(OUT_MP4)
    duration_s = n_states / FPS
    print(
        f"[gate_latch_pocket_oriented_movie] wrote {OUT_MP4} "
        f"({size_bytes / 1e6:.1f} MB, {n_states} frames, "
        f"~{duration_s:.1f} s @ {FPS} fps)"
    )

    if not KEEP_FRAMES:
        shutil.rmtree(frame_dir, ignore_errors=True)
        print(f"[gate_latch_pocket_oriented_movie] cleaned up temp frame dir {frame_dir}")
    else:
        print(f"[gate_latch_pocket_oriented_movie] kept temp frame dir {frame_dir} (GATE_LATCH_KEEP_FRAMES=1)")


def main():
    """Runs load/fit/style/orient/render for this invocation's frame range, then stitches.

    See module docstring "Chunked rendering". Three modes, all driven by
    env vars:
      - STITCH_ONLY: skip loading/rendering entirely, just ffmpeg-stitch
        whatever PNGs already exist in FRAME_DIR (from prior chunk runs).
      - Chunked (FRAME_START/FRAME_STOP set narrower than the full range,
        typically with SKIP_STITCH=1): load everything (cheap -- states are
        just coordinate arrays) but only render+PNG that sub-range, then
        skip stitching so the process can exit and free its memory before
        the next chunk starts.
      - Default (no env vars): identical to the original single-process
        behavior -- load, render every frame, stitch, done.
    """
    if STITCH_ONLY:
        # Count actual PNGs rather than trust N_OUTPUT_FRAMES for this
        # report -- that config value (computed from RAW_N_FRAMES_TOTAL/
        # STRIDE before any trajectory is actually loaded) is off by one
        # from what load_traj really produces for this xtc (401 vs the 400
        # states actually loaded; harmless everywhere else since
        # render_frames() always used the real n_states, not this config
        # value, for its loop bound -- but STITCH_ONLY never loads the
        # trajectory, so there's no real n_states available here).
        n_frames_found = len(glob.glob(os.path.join(FRAME_DIR, "frame_*.png")))
        print(
            f"[gate_latch_pocket_oriented_movie] STITCH_ONLY mode: stitching "
            f"{n_frames_found} existing frames in {FRAME_DIR}"
        )
        stitch_movie(FRAME_DIR, n_frames_found)
        return

    print(
        f"[gate_latch_pocket_oriented_movie] mode={'SMOKE TEST' if IS_SMOKE_TEST else 'FULL RUN'}, "
        f"stride={STRIDE}, n_output_frames={N_OUTPUT_FRAMES}, "
        f"frame_range=[{FRAME_START},{FRAME_STOP}], skip_stitch={SKIP_STITCH}, "
        f"fps={FPS}, out={OUT_MP4}"
    )

    setup_global_render_settings()
    n_states = load_topology_and_trajectory()
    define_selections()
    remove_rigid_body_drift()
    style_scene()
    anchor = orient_gate_left_latch_right()
    set_fixed_camera(anchor)

    os.makedirs(FRAME_DIR, exist_ok=True)
    print(
        f"[gate_latch_pocket_oriented_movie] rendering states "
        f"{FRAME_START}-{FRAME_STOP} of {n_states} into {FRAME_DIR}"
    )
    render_frames(n_states, FRAME_DIR, FRAME_START, FRAME_STOP)

    if SKIP_STITCH:
        print(
            "[gate_latch_pocket_oriented_movie] SKIP_STITCH set -- leaving frames "
            f"in {FRAME_DIR} for a later GATE_LATCH_STITCH_ONLY=1 run"
        )
    else:
        stitch_movie(FRAME_DIR, n_states)


main()
