#!/usr/bin/env python
"""Renders a movie of gate/latch closure dynamics for pair_3101_binder.

PYR1+LCA complex, for a research talk.

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

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/gate_latch_movie.py

Smoke test (renders only a handful of frames, fast, writes to a
"_smoketest" suffixed mp4 so it never clobbers the real deliverable):

    GATE_LATCH_MAX_FRAMES=25 /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/gate_latch_movie.py

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

# Consistent landmark colors (see CLAUDE.md IBM colorblind-safe palette /
# pymol_renders/medoid_comparison.py convention)
GATE_COLOR_HEX = "#FE6100"   # orange
LATCH_COLOR_HEX = "#785EF0"  # purple
PROTEIN_GRAY = "gray80"

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

OUT_NAME = "gate_latch_closure_pair3101_binder_0-30ns"
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
    print(f"[gate_latch_movie] loaded topology: {n_atoms_pdb} atoms from {PDB}")

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
        f"[gate_latch_movie] loaded {n_states} states from trajectory "
        f"(stride={STRIDE}, raw stop frame={RAW_STOP} of {RAW_N_FRAMES_TOTAL})"
    )

    # Strip hydrogens for lighter/faster cartoon+stick rendering across
    # all loaded states (removal applies across all already-loaded states).
    cmd.remove("mol and hydro")
    print(f"[gate_latch_movie] after removing hydrogens: {cmd.count_atoms('mol')} atoms")

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
        f"[gate_latch_movie] selections: core_fit={n_core} atoms, "
        f"gate={n_gate} atoms, latch={n_latch} atoms, ligand={n_lig} atoms"
    )
    if n_core == 0 or n_gate == 0 or n_latch == 0 or n_lig == 0:
        raise RuntimeError(
            "[gate_latch_movie] one or more required selections is empty; "
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
        f"[gate_latch_movie] intra_fit on core_fit_sel done; "
        f"max per-state RMSD to state 1 = {max(rms_list):.3f} A"
    )


def style_scene():
    """Applies cartoon/stick styling: gray protein, colored gate/latch, ligand."""
    cmd.hide("everything", "mol")

    cmd.show("cartoon", "protein_sel")
    cmd.color(PROTEIN_GRAY, "protein_sel")
    cmd.set("cartoon_transparency", 0.15, "protein_sel")

    cmd.set_color("gate_color", list(hex_to_rgb01(GATE_COLOR_HEX)))
    cmd.set_color("latch_color", list(hex_to_rgb01(LATCH_COLOR_HEX)))

    cmd.color("gate_color", "gate_sel")
    cmd.color("latch_color", "latch_sel")
    # Keep gate/latch fully opaque even though the rest of the protein is
    # translucent, so the closing loops stay the visual focal point.
    cmd.set("cartoon_transparency", 0.0, "gate_sel or latch_sel")
    cmd.set("cartoon_loop_radius", 0.4, "gate_sel or latch_sel")
    cmd.set("cartoon_tube_radius", 0.4, "gate_sel or latch_sel")

    cmd.show("sticks", "(gate_sel or latch_sel) and not name C+N+O")
    cmd.set("stick_radius", 0.18, "gate_sel or latch_sel")
    cmd.color("gate_color", "gate_sel and elem C")
    cmd.color("latch_color", "latch_sel and elem C")

    cmd.show("sticks", "ligand_sel")
    cmd.set("stick_radius", 0.22, "ligand_sel")
    util.cbay("ligand_sel")

    cmd.deselect()


def add_time_label_anchor():
    """Places a fixed-in-model-space pseudoatom to anchor the time HUD label.

    Positioned just outside the ligand along the (stable-core -> ligand)
    direction at state 1, mirroring the atomic-coordinate-only placement
    trick in pymol_renders/medoid_comparison.py::add_pocket_label. No
    camera/view matrix decoding is needed: since the coordinate frame is
    drift-corrected (remove_rigid_body_drift) and the camera never moves, a
    fixed model-space point stays fixed on screen for the whole movie.

    Returns:
        str: Name of the created pseudoatom object ("time_label_anchor").
    """
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
    # Pushed well clear of the ligand/pocket mouth so the label doesn't
    # visually collide with the gate/latch closure action itself.
    label_pos = _add(ligand_c, outward, ligand_radius + 14.0)

    anchor = "time_label_anchor"
    cmd.pseudoatom(anchor, pos=list(label_pos))
    cmd.hide("everything", anchor)
    cmd.set("label_size", 30, anchor)
    cmd.set("label_color", "black", anchor)
    cmd.set("label_font_id", 7, anchor)
    cmd.set("label_outline_color", "white", anchor)
    return anchor


def set_fixed_camera(anchor):
    """Orients and zooms the camera exactly once, on state 1.

    No camera-moving command is called anywhere else in this script (see
    module docstring, design decision 2).

    buffer is larger than the medoid_comparison.py reference (5 A) as a
    safety margin. pair_3101_binder state 1 (t=0) sits at or near the
    open-baseline peak (~0.95-1.06 nm gate-latch distance, per the smoothed
    timeseries), and the rest of the window (closing to ~0.73 nm, then a
    stable closed band of ~0.70-0.78 nm through t=30 ns) stays at or inside
    that span, so this sequence shouldn't need extra headroom beyond the
    t=0 framing snapshot -- unlike the earlier seq14_binder window. The
    larger buffer is kept anyway since it's harmless (just slightly more
    zoomed out); check the smoke test and tighten back toward 5 for a
    closer framing if needed.

    Args:
        anchor (str): Pseudoatom name from add_time_label_anchor, included
            in the zoom selection so the HUD label stays in frame.
    """
    cmd.frame(1)
    cmd.orient("landmark_sel")
    cmd.zoom(f"landmark_sel or {anchor}", buffer=8)


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
            print(f"[gate_latch_movie] rendered frame {s}/{n_states} (t={sim_ns:.1f} ns)")


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
    print("[gate_latch_movie] running:", " ".join(ffmpeg_cmd))
    result = subprocess.run(ffmpeg_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(
            f"[gate_latch_movie] ffmpeg failed (exit {result.returncode}); "
            f"frames left in {frame_dir} for debugging."
        )

    size_bytes = os.path.getsize(OUT_MP4)
    duration_s = n_states / FPS
    print(
        f"[gate_latch_movie] wrote {OUT_MP4} "
        f"({size_bytes / 1e6:.1f} MB, {n_states} frames, "
        f"~{duration_s:.1f} s @ {FPS} fps)"
    )

    if not KEEP_FRAMES:
        shutil.rmtree(frame_dir, ignore_errors=True)
        print(f"[gate_latch_movie] cleaned up temp frame dir {frame_dir}")
    else:
        print(f"[gate_latch_movie] kept temp frame dir {frame_dir} (GATE_LATCH_KEEP_FRAMES=1)")


def main():
    """Runs the full pipeline: load, fit, style, camera, render, stitch."""
    print(
        f"[gate_latch_movie] mode={'SMOKE TEST' if IS_SMOKE_TEST else 'FULL RUN'}, "
        f"stride={STRIDE}, n_output_frames={N_OUTPUT_FRAMES}, fps={FPS}, "
        f"out={OUT_MP4}"
    )

    setup_global_render_settings()
    n_states = load_topology_and_trajectory()
    define_selections()
    remove_rigid_body_drift()
    style_scene()
    anchor = add_time_label_anchor()
    set_fixed_camera(anchor)

    frame_dir = tempfile.mkdtemp(prefix="gate_latch_frames_")
    print(f"[gate_latch_movie] rendering {n_states} frames into {frame_dir}")
    render_frames(n_states, frame_dir)
    stitch_movie(frame_dir, n_states)


main()
