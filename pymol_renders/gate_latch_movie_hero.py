#!/usr/bin/env python
"""Renders a movie of gate/latch closure dynamics for pair_3101_binder.

PYR1+LCA complex, restyled as a richer, more "hero shot"-like companion
piece to pymol_renders/pyr1_lca_hero_complex.py (which produced
pymol_renders/output/pyr1_lca_biosensor.png): a cartoon+translucent-surface
hybrid, the same blue/orange/pink landmark palette, and a cyan
ball-and-stick ligand, instead of the original gate_latch_movie.py's flat
gray cartoon + yellow sticks. It is a new, independent file --
gate_latch_movie.py and its output MP4 are left untouched (see that
script's docstring for the un-restyled version's own design history). Where
this script's own choices depart from the original (surface layer,
palette, ligand treatment, and the transparency-vs-occlusion tradeoff that
follows from adding a surface to a *close-up, moving* shot rather than a
static whole-fold one), that reasoning is recorded below rather than
repeating gate_latch_movie.py's still-applicable design decisions
(rigid-body-drift removal, fixed camera, on-screen time label, OpenGL-speed
rendering) verbatim.

This script targets the t=0-30 ns window
(gate_closure_pair3101_0_30ns_PL.xtc) of pair_3101_binder, chosen after
scanning gate_latch_timeseries.xvg across ~38 binder sequences for a
cleaner, more visually obvious closure event than seq14_binder's messy
multi-modal one (seq14_binder had a closing transition immediately
followed by a wide partial-reopening excursion -- see this file's git
history for that version). Smoothed (10-point rolling average) gate-latch
distance for pair_3101_binder: open baseline ~0.95-1.06 nm over t=0-2 ns
(peak 1.06 nm at t=0), a sharp, essentially monotonic closing transition
down to ~0.73 nm over t=2-4 ns, then a stable closed state in a tight
~0.70-0.78 nm band (std ~0.02-0.03 nm) from t=4 ns through t=33 ns -- a
real settled closed state, not breathing between open/closed like
seq14_binder. A second reopening starts around t~34 ns, outside this
window: a single clean closure transition, better material for "clearly
shows the gate closure event."

Mechanistic story being illustrated: in binders, the gate loop (resi
84-90) and latch (resi 114-118) close down over the bound ligand (LIG).
This script renders every Nth frame as a fast (non-ray-traced) PNG with a
fixed camera, then stitches the PNGs into an MP4 with ffmpeg, so any
apparent motion on screen is real conformational change in the loops, not
camera movement or rigid-body tumbling of the whole protein.

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/gate_latch_movie_hero.py

Smoke test (renders only a handful of frames, fast, writes to a
"_smoketest" suffixed mp4 so it never clobbers the real deliverable):

    GATE_LATCH_MAX_FRAMES=25 /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/gate_latch_movie_hero.py

Other environment variable overrides (all optional):
    GATE_LATCH_STRIDE: raw-trajectory-frame stride between rendered output
        frames (default 2; the source trajectory has 801 raw frames at
        37.5 ps spacing spanning t=0-30 ns, so stride 2 gives ~401 output
        frames, ~40.1 s at 10 fps -- 0.75 ns of sim time per second of
        video, the same slow pace that worked for the prior seq14_binder
        0-25ns version, landing in the requested ~30-40 s total range).
    GATE_LATCH_MAX_FRAMES: if set, caps the number of output frames
        loaded/rendered (smoke-test mode; also adds a "_smoketest" suffix
        to the output filename).
    GATE_LATCH_FPS: output movie frame rate (default 10).
    GATE_LATCH_KEEP_FRAMES: "1" to keep the temp PNG frame directory after
        a successful ffmpeg run (default "0", deletes).
    FFMPEG_BIN: path to the ffmpeg binary (default /opt/homebrew/bin/ffmpeg).

Design decisions worth knowing about if you're re-running or adapting this:

  1. Removing rigid-body drift before rendering. The trajectory is a
     pre-trimmed, PBC-handled protein+ligand-only .xtc, but GROMACS PBC
     removal does not remove the protein's own rigid-body translation/
     rotation. All loaded states are fit to state 1 (t=0 ns for this
     window) with cmd.intra_fit() before any camera/orient calls, so the
     movie shows only internal conformational change, not tumbling. The
     fit selection ("core_fit_sel") is protein CA atoms excluding the gate
     (84-90), latch (114-118), and the Lb7a5/C-terminal recoil region
     (148-166; see CLAUDE.md "Key domain conventions") -- fitting on the
     stable scaffold while leaving the flexible/functional loops we want
     to see moving out of the least-squares fit, so their motion isn't
     partially averaged into the global superposition.

  2. Fixed camera. cmd.orient() + cmd.zoom() are called exactly once, on
     state 1, on the gate+latch+ligand region. No camera-moving command
     (orient/zoom/turn/move) is called again in the per-frame render loop
     -- only cmd.frame(state) changes, swapping which already-fitted
     coordinate set is displayed under the same camera.

  3. On-screen running time label. Rather than PyMOL view-matrix decoding
     (fragile, version-dependent) to place a fixed 2D HUD label, this
     script reuses the atomic-coordinate-only placement trick from
     pymol_renders/medoid_comparison.py: a label anchor pseudoatom is
     placed at a fixed point in the (now drift-corrected) model coordinate
     frame, offset outward from the ligand along the protein-core ->
     ligand vector. Because the coordinate frame is drift-corrected (step
     1) and the camera never moves (step 2), that fixed model-space point
     stays fixed on screen for the whole movie. The label text is
     re-issued every frame with the real simulation time in ns, computed
     from frame index * stride * 37.5 ps + T0_NS (not from PyMOL's
     movie-frame counter or the ray/opengl clock). T0_NS = 0.0 for this
     0-30 ns window (it was 40.0 for the original 40-500 ns full-trajectory
     version of this script).

  4. Speed. All output frames are captured with PyMOL's OpenGL renderer
     (cmd.png(..., ray=0) against a fixed cmd.viewport()), never cmd.ray()
     -- ray-tracing every frame individually would take far too long, even
     at this window's smaller frame count (~401 frames vs. ~1023 for the
     original full-trajectory version).

  5. Restyle vs. pyr1_lca_biosensor.png -- camera and transparency levels
     went through one extra iteration to get right. The first attempt kept
     gate_latch_movie.py's original tight camera (cmd.orient("landmark_sel")
     + buffer=8 zoom on just gate+latch+ligand) and added the surface layer
     on top. Checked against a smoke-test render, that combination was
     unusable: at that close a distance, the translucent surface's own
     local curvature filled the entire frame like a foggy wall, hiding the
     cartoon, the gate/latch loops, and most of the ligand behind it --
     exactly the opposite of "more complex and aesthetically pleasing." The
     fix was to stop trying to compensate with transparency levels alone
     and instead adopt pyr1_lca_hero_complex.py's whole-protein
     pocket-opening-axis camera (set_pocket_opening_view(), ported into
     this script) with a zoom target of the whole protein+ligand, not just
     the landmark loops (see set_fixed_camera()'s own "v2" docstring note).
     With that wider framing, CARTOON_TRANSPARENCY and SURFACE_TRANSPARENCY
     could go back to matching pyr1_lca_hero_complex.py's values exactly,
     since the same "wide enough for the surface to read as glass"
     reasoning applies. The tradeoff: the movie now shows gate/latch
     closure in the context of the whole fold rather than as a tight
     mechanistic close-up -- a deliberate choice, since the ask was to
     match pyr1_lca_biosensor.png's structure. Surface geometry is
     recomputed every frame (cmd.show("surface", ...) plus the per-state
     coordinate swap in render_frames()); watch the smoke test's
     wall-clock time before committing to a full ~401-frame run, since a
     translucent surface is more expensive per frame than the original
     cartoon-only scene.

     v3 (current): by request, after watching the v2 render -- the
     gate/latch closure itself was hard to follow at v2's whole-fold
     framing (the loops are a small feature of a large view). Fixed by
     going back to a tight zoom on landmark_sel (gate+latch+ligand only,
     ZOOM_BUFFER=6, tighter than the original script's buffer=8) and
     dropping the surface layer entirely (see style_scene()) -- the same
     "surface fills the frame like fog at close range" problem from the
     first v1 attempt would recur otherwise, and this time there's no
     wider-camera escape hatch, since the whole point is to be close. The
     pocket-opening-axis rotation (set_pocket_opening_view()) is kept even
     without the surface, since it's still a better angle than a plain
     orient(); only the zoom target and buffer changed. CARTOON_TRANSPARENCY
     was also raised (0.0 -> 0.25) so the base cartoon doesn't compete with
     the now more-prominent opaque gate/latch loops at this closer
     distance. The HUD label anchor's bounding reference also had to move
     from the whole protein to landmark_sel, since a whole-protein radius
     placed the label off-screen once the camera was cropped tight (see
     add_time_label_anchor()).
"""

import os
import shutil
import subprocess
import tempfile

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

# t=0-30 ns closure-event window for pair_3101_binder, extracted using the
# same Protein_LIG index group convention as the earlier seq14_binder
# windows. Atom count (2956, confirmed matching this sequence's own
# medoid_PL.pdb via grep count of ATOM records) differs from seq14_binder's
# 2961 -- pair_3101_binder has different mutations, so atom counts don't
# carry over across sequences. Lives in this repo's scratch dir, not
# OneDrive, since it's a derived intermediate.
XTC = os.path.join(
    REPO_ROOT,
    "pymol_renders",
    "scratch_traj",
    "gate_closure_pair3101_0_30ns_PL.xtc",
)

# 801 raw frames at 37.5 ps spacing, spanning simulation time 0-30 ns. Raw
# frame index 0 (1-indexed frame 1 in PyMOL's load_traj convention) = t = 0.0
# ns (unlike the old 40-500 ns full-trajectory version of this script).
RAW_N_FRAMES_TOTAL = 801
TIME_SPACING_NS = 0.0375  # 37.5 ps
T0_NS = 0.0

PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
LIGAND_RESN = "LIG"

# Structural landmarks (see CLAUDE.md "Key domain conventions")
GATE_RESI = "84-90"
LATCH_RESI = "114-118"
# Excluded from the rigid-body-drift-removal fit selection (see docstring
# point 1): gate, latch, Lb7a5 loop, and the C-terminal recoil helix.
FLEXIBLE_EXCLUDE_RESI = "84-90+114-118+148-166"

# Landmark colors -- matches pymol_renders/pyr1_lca_hero_complex.py /
# pyr1_lca_biosensor.png (gate orange is CLAUDE.md's standard landmark
# color; latch pink is a by-request departure from CLAUDE.md's usual
# landmark purple, #785EF0, first made for that still image -- see this
# script's docstring "Restyle" note).
GATE_COLOR_HEX = "#FE6100"     # orange, CLAUDE.md landmark color
LATCH_COLOR_HEX = "#DC267F"    # pink, matches pyr1_lca_biosensor.png
BASE_PROTEIN_HEX = "#5B7FBE"   # medium blue, matches pyr1_lca_biosensor.png
                                # (replaces the original script's flat
                                # gray80)
LIGAND_CARBON_HEX = "#00E5FF"  # bright cyan, matches
                                # pyr1_lca_biosensor.png (replaces the
                                # original script's util.cbay yellow)

# Base (non-gate/latch) protein cartoon transparency. v3 (current): tuned
# for the tight gate/latch-region zoom (see this script's docstring "v3"
# note) -- translucent enough that the base cartoon doesn't visually
# compete with or partially hide the opaque gate/latch loops from a close
# camera distance, but present enough to still read as "the rest of the
# fold" rather than disappearing. Gate/latch themselves are always forced
# fully opaque regardless (see style_scene()). Higher than the original
# gate_latch_movie.py's 0.15 since this base color is a more visually
# present blue rather than a faint gray80.
CARTOON_TRANSPARENCY = 0.25

# Zoom buffer (Angstroms) around landmark_sel (gate+latch+ligand) for the
# tight v3 camera (see set_fixed_camera()'s "v3" note). Tighter than the
# original gate_latch_movie.py's buffer=8 -- since this restyle's whole
# point is a closer, easier-to-follow view of the closure motion, not just
# a like-for-like recreation of the original framing.
ZOOM_BUFFER = 6

# Ligand geometry (ball-and-stick, matches pyr1_lca_hero_complex.py)
LIGAND_STICK_RADIUS = 0.3
LIGAND_SPHERE_SCALE = 0.26

VIEWPORT_W = 1600
VIEWPORT_H = 900

FFMPEG_BIN = os.environ.get("FFMPEG_BIN", "/opt/homebrew/bin/ffmpeg")

# --------------------------------------------------------------------------
# Env-var-driven run parameters (see module docstring for smoke-test usage)
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

# 1-indexed raw-trajectory frame number to stop load_traj at.
RAW_STOP = min(RAW_N_FRAMES_TOTAL, (N_OUTPUT_FRAMES - 1) * STRIDE + 1)

OUT_NAME = "gate_latch_closure_pair3101_binder_0-30ns_hero"
if IS_SMOKE_TEST:
    OUT_NAME += "_smoketest"
OUT_MP4 = os.path.join(OUT_DIR, OUT_NAME + ".mp4")


# --------------------------------------------------------------------------
# Small color / vector helpers (pure python, no numpy dependency) -- same
# pattern as pymol_renders/medoid_comparison.py
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
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 1)
    cmd.set("orthoscopic", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("ray_trace_fog", 0)
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.2)
    cmd.set("ambient", 0.4)
    # Many-state memory/perf: build cartoon geometry per-frame instead of
    # caching every state's geometry simultaneously.
    cmd.set("defer_builds_mode", 3)


def load_topology_and_trajectory():
    """Loads the topology PDB and trajectory, strips hydrogens.

    Returns:
        int: Number of trajectory states loaded.
    """
    cmd.load(PDB, "mol")
    n_atoms_pdb = cmd.count_atoms("mol")
    print(f"[gate_latch_movie_hero] loaded topology: {n_atoms_pdb} atoms from {PDB}")

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
        f"[gate_latch_movie_hero] loaded {n_states} states from trajectory "
        f"(stride={STRIDE}, raw stop frame={RAW_STOP} of {RAW_N_FRAMES_TOTAL})"
    )

    # Strip hydrogens for lighter/faster cartoon+stick rendering across
    # all loaded states (removal applies across all already-loaded states).
    cmd.remove("mol and hydro")
    print(f"[gate_latch_movie_hero] after removing hydrogens: {cmd.count_atoms('mol')} atoms")

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

    cmd.select("protein_sel", protein_sel)
    cmd.select("gate_sel", gate_sel)
    cmd.select("latch_sel", latch_sel)
    cmd.select("ligand_sel", ligand_sel)
    cmd.select("landmark_sel", "gate_sel or latch_sel or ligand_sel")
    cmd.select("core_fit_sel", core_fit_sel)

    n_core = cmd.count_atoms("core_fit_sel")
    n_gate = cmd.count_atoms("gate_sel")
    n_latch = cmd.count_atoms("latch_sel")
    n_lig = cmd.count_atoms("ligand_sel")
    print(
        f"[gate_latch_movie_hero] selections: core_fit={n_core} atoms, "
        f"gate={n_gate} atoms, latch={n_latch} atoms, ligand={n_lig} atoms"
    )
    if n_core == 0 or n_gate == 0 or n_latch == 0 or n_lig == 0:
        raise RuntimeError(
            "[gate_latch_movie_hero] one or more required selections is empty; "
            "check chain/resi/resn conventions before rendering."
        )


def remove_rigid_body_drift():
    """Fits every loaded state onto state 1 using stable-core CA atoms only.

    Gate/latch/Lb7a5/recoil loops are excluded from the fit selection, so
    apparent motion in the movie is real conformational change, not
    whole-protein tumbling. intra_fit applies the resulting rigid transform
    to the whole object (all atoms, all selections) per state, so the
    ligand and gate/latch loops are carried along correctly relative to the
    rest of the protein.
    """
    rms_list = cmd.intra_fit("core_fit_sel", state=1)
    print(
        f"[gate_latch_movie_hero] intra_fit on core_fit_sel done; "
        f"max per-state RMSD to state 1 = {max(rms_list):.3f} A"
    )


def style_scene():
    """Applies cartoon/stick styling: blue protein, colored gate/latch, ligand.

    Matches pyr1_lca_hero_complex.py's blue/orange/pink/cyan palette, but
    (as of v3) without that script's translucent molecular surface -- see
    this script's docstring "v3" note for why the surface doesn't survive
    the switch to a tight, closure-motion-focused zoom.
    """
    cmd.hide("everything", "mol")

    cmd.set_color("base_color", list(hex_to_rgb01(BASE_PROTEIN_HEX)))
    cmd.set_color("gate_color", list(hex_to_rgb01(GATE_COLOR_HEX)))
    cmd.set_color("latch_color", list(hex_to_rgb01(LATCH_COLOR_HEX)))

    cmd.color("base_color", "protein_sel")
    cmd.color("gate_color", "gate_sel")
    cmd.color("latch_color", "latch_sel")

    cmd.show("cartoon", "protein_sel")
    cmd.set("cartoon_transparency", CARTOON_TRANSPARENCY, "protein_sel")
    # Keep gate/latch fully opaque even though the rest of the protein is
    # translucent, so the closing loops stay the visual focal point.
    cmd.set("cartoon_transparency", 0.0, "gate_sel or latch_sel")
    cmd.set("cartoon_loop_radius", 0.4, "gate_sel or latch_sel")
    cmd.set("cartoon_tube_radius", 0.4, "gate_sel or latch_sel")

    # No molecular surface in this version -- see this script's docstring
    # "v3" note. A whole-protein surface read as glass under the wide
    # pocket-opening framing, but once the camera is zoomed tight on just
    # the gate/latch/ligand region (for closure-motion clarity), that same
    # surface's local curvature fills the frame and hides the loops it's
    # supposed to be showing. Dropping it was the direct fix.

    cmd.show("sticks", "(gate_sel or latch_sel) and not name C+N+O")
    cmd.set("stick_radius", 0.18, "gate_sel or latch_sel")
    cmd.color("gate_color", "gate_sel and elem C")
    cmd.color("latch_color", "latch_sel and elem C")

    # Ligand: ball-and-stick with cyan carbons, matching
    # pyr1_lca_hero_complex.py (replaces the original script's plain
    # yellow-carbon sticks, util.cbay).
    cmd.set_color("lig_carbon_color", list(hex_to_rgb01(LIGAND_CARBON_HEX)))
    cmd.show("sticks", "ligand_sel")
    cmd.show("spheres", "ligand_sel")
    cmd.set("stick_radius", LIGAND_STICK_RADIUS, "ligand_sel")
    cmd.set("sphere_scale", LIGAND_SPHERE_SCALE, "ligand_sel")
    cmd.color("lig_carbon_color", "ligand_sel")
    util.cnc("ligand_sel")

    cmd.deselect()


def set_pocket_opening_view():
    """Builds a camera rotation looking into the pocket from the gate/latch end.

    Ported from pyr1_lca_hero_complex.py's set_pocket_opening_view() (see
    that script's module docstring "Camera" section for the full gluLookAt
    derivation): outward_dir = normalize(gate_latch_centroid -
    protein_core_centroid) at state 1, a camera basis built from it, and
    only the rotation block written via cmd.set_view() (translation/scale
    left for cmd.zoom() in set_fixed_camera()). Computed on the
    already-drift-corrected coordinate frame (remove_rigid_body_drift() has
    already run), so the resulting fixed rotation stays valid for the whole
    movie exactly as the original tight-zoom camera did.

    Returns:
        tuple[bool, tuple[float, float, float] | None]: (True, up_vector)
        if the deterministic view was set, where up_vector is the model-
        space "up" direction of the resulting camera (used by
        add_time_label_anchor() to place the HUD label above the molecule
        on screen); (False, None) if it fell back (missing atoms in the
        core or gate+latch CA selection -- shouldn't happen given
        define_selections()'s atom-count checks).
    """
    protein_ca_sel = "protein_sel and name CA"
    gate_latch_ca_sel = f"({protein_ca_sel}) and (resi {GATE_RESI} or resi {LATCH_RESI})"

    protein_model = cmd.get_model(protein_ca_sel, state=1)
    gate_latch_model = cmd.get_model(gate_latch_ca_sel, state=1)

    if not protein_model.atom or not gate_latch_model.atom:
        print(
            "[gate_latch_movie_hero] WARNING: could not compute "
            "deterministic pocket-opening camera; falling back to "
            "cmd.orient()"
        )
        return False, None

    protein_core_centroid = _centroid(protein_model)
    gate_latch_centroid = _centroid(gate_latch_model)
    outward_dir = _norm(_sub(gate_latch_centroid, protein_core_centroid))

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
    return True, up


def add_time_label_anchor(up_vector):
    """Places a fixed-in-model-space pseudoatom to anchor the time HUD label.

    v2 (see this script's docstring "Restyle" note): positioned above the
    whole protein+ligand along the camera's own "up" vector, so the label
    sits above the molecule's silhouette on screen regardless of viewing
    angle. v1 (ported unchanged from gate_latch_movie.py) placed it just
    outside the ligand along the stable-core -> ligand direction -- that
    made sense for the original tight landmark-only camera, but under the
    new whole-protein pocket-opening camera (set_pocket_opening_view()) it
    placed the label right on top of the latch loop's on-screen position,
    overlapping the sticks. No camera/view matrix decoding is needed
    either way: since the coordinate frame is drift-corrected
    (remove_rigid_body_drift) and the camera never moves, a fixed
    model-space point stays fixed on screen for the whole movie.

    Args:
        up_vector: Model-space "up" direction of the fixed camera, from
            set_pocket_opening_view(). If None (that view's deterministic
            derivation fell back), falls back to v1's ligand-outward
            placement instead, since there's no camera-relative direction
            to place it along.

    Returns:
        str: Name of the created pseudoatom object ("time_label_anchor").
    """
    anchor = "time_label_anchor"

    if up_vector is not None:
        # v3: bounding reference is landmark_sel (gate+latch+ligand), not
        # the whole protein -- now that the camera is zoomed tight on just
        # that region (see set_fixed_camera()'s "v3" note), a whole-protein
        # radius would place the label way outside the cropped frame.
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
        # Pushed well clear of the landmark region's silhouette along the
        # camera's own "up" axis, so the label never overlaps the
        # cartoon/gate/latch regardless of how the loops move.
        label_pos = _add(center, up_vector, radius + 8.0)
    else:
        core_model = cmd.get_model("core_fit_sel", state=1)
        ligand_model = cmd.get_model("ligand_sel", state=1)
        core_c = _centroid(core_model)
        ligand_c = _centroid(ligand_model)
        outward = _norm(_sub(ligand_c, core_c))
        ligand_radius = max(
            (
                (a.coord[0] - ligand_c[0]) ** 2
                + (a.coord[1] - ligand_c[1]) ** 2
                + (a.coord[2] - ligand_c[2]) ** 2
            )
            ** 0.5
            for a in ligand_model.atom
        )
        label_pos = _add(ligand_c, outward, ligand_radius + 14.0)

    cmd.pseudoatom(anchor, pos=list(label_pos))
    cmd.hide("everything", anchor)
    cmd.set("label_size", 30, anchor)
    cmd.set("label_color", "black", anchor)
    cmd.set("label_font_id", 7, anchor)
    cmd.set("label_outline_color", "white", anchor)
    return anchor


def set_fixed_camera():
    """Orients and zooms the camera exactly once, on state 1.

    No camera-moving command is called anywhere else in this script (see
    module docstring, design decision 2).

    v1 of this restyle used cmd.orient("landmark_sel") + a tight buffer=8
    zoom on just landmark_sel (gate+latch+ligand), unchanged from the
    original gate_latch_movie.py. Seen against a full-protein translucent
    surface, that tight framing put the camera close enough that the
    surface's own local curvature filled the whole viewport like a foggy
    wall, hiding the cartoon, the gate/latch loops, and most of the ligand
    behind it. v2 fixed that by reusing pyr1_lca_hero_complex.py's
    pocket-opening-axis camera (set_pocket_opening_view()) and zooming on
    the whole protein+ligand instead -- matching pyr1_lca_biosensor.png's
    structure, but by request (the gate/latch closure itself read as too
    small/hard to follow at that wider framing) v3 (current) drops the
    surface entirely (see style_scene()) and goes back to a tight zoom, so
    the loop motion is the dominant thing on screen. The pocket-opening
    rotation is kept even without the surface, since it's still a better
    camera angle than a plain cmd.orient() (looking into the pocket from
    the gate/latch side rather than face-on into the beta-sheet core, see
    pyr1_lca_hero_complex.py's "Camera" docstring section) -- only the zoom
    *target* reverts to landmark_sel-only, tighter than v1's buffer=8 (see
    ZOOM_BUFFER).

    Creates the HUD label anchor itself (add_time_label_anchor()) after
    the camera rotation is set, since v2/v3's anchor placement needs the
    camera's "up" vector (see that function's docstring).

    Returns:
        str: Name of the created time-label pseudoatom (from
        add_time_label_anchor()), for render_frames() to re-label per
        frame.
    """
    cmd.frame(1)
    view_was_set, up_vector = set_pocket_opening_view()
    if not view_was_set:
        cmd.orient("landmark_sel")

    anchor = add_time_label_anchor(up_vector)
    cmd.zoom(f"landmark_sel or {anchor}", buffer=ZOOM_BUFFER)
    return anchor


def render_frames(n_states, frame_dir):
    """Renders one OpenGL (non-ray-traced) PNG per trajectory state.

    Args:
        n_states (int): Number of loaded trajectory states to render.
        frame_dir (str): Directory to write frame_00000.png etc. into.
    """
    cmd.viewport(VIEWPORT_W, VIEWPORT_H)

    for s in range(1, n_states + 1):
        cmd.frame(s)
        sim_ns = T0_NS + (s - 1) * STRIDE * TIME_SPACING_NS
        cmd.label("time_label_anchor", '"%d ns"' % int(round(sim_ns)))

        png_path = os.path.join(frame_dir, f"frame_{s - 1:05d}.png")
        cmd.png(png_path, ray=0, quiet=1)

        if s == 1 or s % 25 == 0 or s == n_states:
            print(f"[gate_latch_movie_hero] rendered frame {s}/{n_states} (t={sim_ns:.1f} ns)")


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
    print("[gate_latch_movie_hero] running:", " ".join(ffmpeg_cmd))
    result = subprocess.run(ffmpeg_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(
            f"[gate_latch_movie_hero] ffmpeg failed (exit {result.returncode}); "
            f"frames left in {frame_dir} for debugging."
        )

    size_bytes = os.path.getsize(OUT_MP4)
    duration_s = n_states / FPS
    print(
        f"[gate_latch_movie_hero] wrote {OUT_MP4} "
        f"({size_bytes / 1e6:.1f} MB, {n_states} frames, "
        f"~{duration_s:.1f} s @ {FPS} fps)"
    )

    if not KEEP_FRAMES:
        shutil.rmtree(frame_dir, ignore_errors=True)
        print(f"[gate_latch_movie_hero] cleaned up temp frame dir {frame_dir}")
    else:
        print(f"[gate_latch_movie_hero] kept temp frame dir {frame_dir} (GATE_LATCH_KEEP_FRAMES=1)")


def main():
    """Runs the full pipeline: load, fit, style, camera, render, stitch."""
    print(
        f"[gate_latch_movie_hero] mode={'SMOKE TEST' if IS_SMOKE_TEST else 'FULL RUN'}, "
        f"stride={STRIDE}, n_output_frames={N_OUTPUT_FRAMES}, fps={FPS}, "
        f"out={OUT_MP4}"
    )

    setup_global_render_settings()
    n_states = load_topology_and_trajectory()
    define_selections()
    remove_rigid_body_drift()
    style_scene()
    set_fixed_camera()

    frame_dir = tempfile.mkdtemp(prefix="gate_latch_frames_")
    print(f"[gate_latch_movie_hero] rendering {n_states} frames into {frame_dir}")
    render_frames(n_states, frame_dir)
    stitch_movie(frame_dir, n_states)


main()
