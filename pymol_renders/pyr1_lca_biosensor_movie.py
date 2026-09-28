#!/usr/bin/env python
"""Renders a slow 360-degree turntable movie of the pyr1_lca_biosensor.png structure.

Same PYR1 (lca_001_binder) + LCA complex, same coordinate-matched input
PDBs, and the exact same load/style/starting-camera code as
pymol_renders/pyr1_lca_hero_complex.py (which produced
pymol_renders/output/pyr1_lca_biosensor.png): cartoon+translucent-surface
hybrid, the blue/orange/pink landmark palette, cyan ball-and-stick ligand,
white background, and the pocket-opening-axis camera. Where this script
differs from that one is purely "movie instead of single ray-traced still":
it spins the camera slowly about the vertical screen axis for one full
360-degree turn so every side of the molecule is visible, and renders with
PyMOL's fast OpenGL path instead of cmd.ray() -- see "Render mode" below for
why.

There is no trajectory here (unlike pymol_renders/gate_latch_movie_hero.py,
which animates MD simulation time for pair_3101_binder): this is one static
structure, and "movie" means a turntable rotation of that one frame, not
simulation dynamics. Camera scale/framing is set once, exactly as in
pyr1_lca_hero_complex.py, and the render loop calls only cmd.turn() between
frames, never cmd.zoom()/cmd.orient() again, so the molecule doesn't drift
or resize during the spin.

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/pyr1_lca_biosensor_movie.py

Smoke test (renders only a handful of frames, fast, writes to a
"_smoketest" suffixed mp4 so it never clobbers the real deliverable):

    PYR1_MOVIE_MAX_FRAMES=20 /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/pyr1_lca_biosensor_movie.py

Other environment variable overrides (all optional):
    PYR1_MOVIE_DURATION_S: seconds for one full 360-degree rotation
        (default 30 -- a slow, easy-to-watch turntable pace of 12
        degrees/second; lower this for a faster spin).
    PYR1_MOVIE_FPS: output movie frame rate (default 24; also sets how many
        frames make up the rotation, since n_frames = duration * fps).
    PYR1_MOVIE_MAX_FRAMES: if set, caps the number of output frames
        rendered (smoke-test mode; also adds a "_smoketest" suffix to the
        output filename). The rotation still completes a full 360 degrees
        over the *full* frame count, so a smoke test only shows the first
        slice of the turn, not a fast preview of the whole thing.
    PYR1_MOVIE_KEEP_FRAMES: "1" to keep the temp PNG frame directory after
        a successful ffmpeg run (default "0", deletes).
    FFMPEG_BIN: path to the ffmpeg binary (default /opt/homebrew/bin/ffmpeg).

Render mode -- OpenGL (cmd.png(..., ray=0)), not ray-traced. A timing check
against this exact scene (cartoon+surface hybrid, ray_shadows on,
light_count 3, antialias 2, matching pyr1_lca_hero_complex.py's full
lighting setup) measured ~59 s/frame to ray-trace at 1600x1200 on this
machine -- multiplied out over hundreds of frames for a smooth 360-degree
turn, that's hours, not minutes. The same scene renders via OpenGL in ~2
s/frame. This mirrors gate_latch_movie_hero.py's own "Speed" design
decision (see that script's docstring, point 4): every existing movie in
this repo uses cmd.png(ray=0) against a fixed viewport for exactly this
reason, reserving cmd.ray() for single showcase stills. The tradeoff is
that OpenGL's transparency/lighting approximation reads slightly flatter
than the ray-traced still -- acceptable for a rotating movie watched in
motion, where per-frame photorealism matters less than smooth playback.

Axis of rotation -- "the z (vertical) axis" as the request framed it means
rotating the view around vertical while both top and bottom of the molecule
stay fixed on screen (a standard camera turntable), so every side comes
into view over one full turn. In PyMOL's own cmd.turn(axis, angle) API the
screen-space vertical axis is confusingly named "y" (turn's axes are
screen-relative: "x" pitches the view, "y" yaws it around vertical, "z"
rolls/spins it in-plane like a clock hand, which would NOT reveal new
sides). This script calls cmd.turn("y", ...) each frame to get the intended
vertical-axis turntable motion; see TURN_AXIS below.
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

PROTEIN_PDB = os.path.join(
    REPO_ROOT, "binders/lca_001_binder/protein_lca001_fixed_H.pdb"
)
LIGAND_PDB = os.path.join(
    REPO_ROOT, "binders/lca_001_binder/ligand_lca001.pdb"
)

PROTEIN_OBJ = "pyr1"
LIGAND_OBJ = "lca_ligand"

PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
LIGAND_RESN = "LIG"

# Structural landmarks (see CLAUDE.md "Key domain conventions")
GATE_RESI = "84-90"
LATCH_RESI = "114-118"

# Colors -- identical to pyr1_lca_hero_complex.py / pyr1_lca_biosensor.png
BASE_CARTOON_HEX = "#5B7FBE"
GATE_HEX = "#FE6100"
LATCH_HEX = "#DC267F"
LIGAND_CARBON_HEX = "#00E5FF"
SURFACE_HEX = "#B9CBE8"
BG_HEX = "#FFFFFF"

# Surface / cartoon hybrid -- same values as pyr1_lca_hero_complex.py
SURFACE_TRANSPARENCY = 0.55
CARTOON_TRANSPARENCY = 0.0

# Ligand geometry -- same values as pyr1_lca_hero_complex.py
LIGAND_STICK_RADIUS = 0.3
LIGAND_SPHERE_SCALE = 0.26

# Movie frame size. Smaller than pyr1_lca_hero_complex.py's 2400x1800 still
# (this is a fast OpenGL turntable, not a poster-quality ray-traced render),
# but still comfortably above 1080p-equivalent detail for a rotating shot.
VIEWPORT_W = 1600
VIEWPORT_H = 1200

# Camera framing -- same buffer as pyr1_lca_hero_complex.py, so frame 0 of
# this movie matches pyr1_lca_biosensor.png's framing before the turn starts.
ZOOM_BUFFER = 1.5

# Screen-space turn axis for a vertical-axis turntable (see module docstring
# "Axis of rotation" -- PyMOL's cmd.turn() names this axis "y", not "z").
TURN_AXIS = "y"

FFMPEG_BIN = os.environ.get("FFMPEG_BIN", "/opt/homebrew/bin/ffmpeg")

# --------------------------------------------------------------------------
# Env-var-driven run parameters (see module docstring for smoke-test usage)
# --------------------------------------------------------------------------
DURATION_S = float(os.environ.get("PYR1_MOVIE_DURATION_S", "30"))
FPS = int(os.environ.get("PYR1_MOVIE_FPS", "24"))
KEEP_FRAMES = os.environ.get("PYR1_MOVIE_KEEP_FRAMES", "0") == "1"
_max_frames_env = os.environ.get("PYR1_MOVIE_MAX_FRAMES", "").strip()

FULL_N_FRAMES = max(2, round(DURATION_S * FPS))
DEG_PER_FRAME = 360.0 / FULL_N_FRAMES

if _max_frames_env:
    IS_SMOKE_TEST = True
    N_RENDERED_FRAMES = min(int(_max_frames_env), FULL_N_FRAMES)
else:
    IS_SMOKE_TEST = False
    N_RENDERED_FRAMES = FULL_N_FRAMES

OUT_NAME = "pyr1_lca_biosensor_turntable"
if IS_SMOKE_TEST:
    OUT_NAME += "_smoketest"
OUT_MP4 = os.path.join(OUT_DIR, OUT_NAME + ".mp4")


# --------------------------------------------------------------------------
# Small color / vector helpers (pure python, no numpy dependency) --
# duplicated from pyr1_lca_hero_complex.py rather than imported, matching
# that script's own "stand alone" convention (see its module docstring).
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
# Scene construction -- identical to pyr1_lca_hero_complex.py except the
# global render settings drop the ray-trace-only knobs (ray_shadows,
# ray_trace_mode, ray_interior_color, spec_direct*) that have no effect in
# OpenGL mode, and antialias is reduced from 2 to 1 to keep per-frame time
# down (see module docstring "Render mode").
# --------------------------------------------------------------------------
def setup_global_render_settings():
    """Resets the PyMOL session and applies global render, lighting, and background settings."""
    cmd.reinitialize()

    cmd.set_color("bg_flat", list(hex_to_rgb01(BG_HEX)))
    cmd.bg_color("bg_flat")
    cmd.set("bg_gradient", 0)
    cmd.set("ray_opaque_background", 1)

    cmd.set("antialias", 1)
    cmd.set("orthoscopic", 0)
    cmd.set("field_of_view", 28)

    cmd.set("light_count", 3)
    cmd.set("ambient", 0.28)
    cmd.set("direct", 0.65)
    cmd.set("reflect", 0.45)
    cmd.set("specular", 0.4)
    cmd.set("shininess", 60)
    cmd.set("spec_power", 150)

    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.5)

    cmd.set("transparency_mode", 2)
    cmd.set("two_sided_lighting", 1)

    cmd.set("surface_quality", 1)  # lower than the still's 2 -- built once
                                    # (not per frame, since the turn only
                                    # moves the camera, not the coordinates)
                                    # but no need to pay for extra mesh
                                    # detail that won't be visible at movie
                                    # resolution
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)


def load_structures():
    """Loads the protein and ligand PDBs as two separate objects, one session.

    Raises:
        SystemExit: No protein polymer atoms, or no ligand atoms, loaded.
    """
    cmd.load(PROTEIN_PDB, PROTEIN_OBJ)
    cmd.load(LIGAND_PDB, LIGAND_OBJ)

    cmd.remove(f"{PROTEIN_OBJ} and hydro")
    cmd.remove(f"{PROTEIN_OBJ} and solvent")

    n_protein = cmd.count_atoms(f"{PROTEIN_OBJ} and polymer")
    n_ligand = cmd.count_atoms(f"{LIGAND_OBJ} and resn {LIGAND_RESN}")
    print(
        f"[pyr1_lca_biosensor_movie] loaded protein: {n_protein} atoms "
        f"(chain {PROTEIN_CHAIN} expected); ligand: {n_ligand} atoms "
        f"(resn {LIGAND_RESN})"
    )
    if n_protein == 0:
        raise SystemExit(
            f"[pyr1_lca_biosensor_movie] ABORT: no protein polymer atoms "
            f"loaded from {PROTEIN_PDB}"
        )
    if n_ligand == 0:
        raise SystemExit(
            f"[pyr1_lca_biosensor_movie] ABORT: no resn {LIGAND_RESN} atoms "
            f"loaded from {LIGAND_PDB}"
        )


def style_scene():
    """Styles the protein as a cartoon+surface hybrid and the ligand as ball-and-stick.

    Identical styling to pyr1_lca_hero_complex.py -- see that script's
    module docstring for the full rationale (shared per-atom coloring
    between cartoon and surface, decoupled surface_color for a clean glass
    look, cyan ball-and-stick ligand).
    """
    protein_sel = f"{PROTEIN_OBJ} and chain {PROTEIN_CHAIN} and polymer"
    gate_sel = f"{protein_sel} and resi {GATE_RESI}"
    latch_sel = f"{protein_sel} and resi {LATCH_RESI}"
    ligand_sel = f"{LIGAND_OBJ} and resn {LIGAND_RESN}"

    cmd.select("gate", gate_sel)
    cmd.select("latch", latch_sel)
    cmd.select("lig", ligand_sel)

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

    cmd.hide("everything", LIGAND_OBJ)
    cmd.show("sticks", ligand_sel)
    cmd.show("spheres", ligand_sel)
    cmd.set("stick_radius", LIGAND_STICK_RADIUS, ligand_sel)
    cmd.set("sphere_scale", LIGAND_SPHERE_SCALE, ligand_sel)
    cmd.set_color("lig_carbon_color", list(hex_to_rgb01(LIGAND_CARBON_HEX)))
    cmd.color("lig_carbon_color", ligand_sel)
    util.cnc(ligand_sel)

    cmd.deselect()


def set_pocket_opening_view():
    """Builds a camera rotation looking into the pocket from the gate/latch end.

    Ported unchanged from pyr1_lca_hero_complex.py (see that script's
    module docstring "Camera" section for the full derivation). Only sets
    this movie's *starting* orientation, frame 0 -- frame_camera() zooms
    once after this, and the per-frame render loop then turns away from it.

    Returns:
        bool: True if the deterministic view was set; False if it fell back.
    """
    protein_ca_sel = (
        f"{PROTEIN_OBJ} and chain {PROTEIN_CHAIN} and polymer and name CA"
    )
    gate_latch_ca_sel = f"({protein_ca_sel}) and (resi {GATE_RESI} or resi {LATCH_RESI})"

    protein_model = cmd.get_model(protein_ca_sel)
    gate_latch_model = cmd.get_model(gate_latch_ca_sel)

    if not protein_model.atom or not gate_latch_model.atom:
        print(
            "[pyr1_lca_biosensor_movie] WARNING: could not compute "
            "deterministic pocket-opening camera (missing CA atoms in "
            "protein or gate+latch selection); falling back to cmd.orient()"
        )
        return False

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
    return True


def frame_camera():
    """Sets the pocket-opening starting camera, then zooms on the whole complex.

    Called exactly once, before the render loop. No camera-moving command
    other than cmd.turn(TURN_AXIS, ...) runs after this, so scale/framing
    stay fixed for the whole turntable (only the viewing angle changes).
    """
    got_deterministic_view = set_pocket_opening_view()
    if not got_deterministic_view:
        cmd.orient("lig or gate or latch")

    cmd.zoom(f"{PROTEIN_OBJ} or {LIGAND_OBJ}", buffer=ZOOM_BUFFER)


def render_frames(frame_dir):
    """Renders one OpenGL (non-ray-traced) PNG per turntable step.

    Turns the camera TURN_AXIS ("y", PyMOL's name for the screen-vertical
    axis -- see module docstring "Axis of rotation") by DEG_PER_FRAME
    before each frame, so frame 0 is the pocket-opening starting view set
    by frame_camera() and the sequence sweeps through a full 360 degrees
    over FULL_N_FRAMES (only N_RENDERED_FRAMES of which are actually
    written out in smoke-test mode).

    Args:
        frame_dir (str): Directory to write frame_00000.png etc. into.
    """
    cmd.viewport(VIEWPORT_W, VIEWPORT_H)

    for i in range(N_RENDERED_FRAMES):
        if i > 0:
            cmd.turn(TURN_AXIS, DEG_PER_FRAME)

        png_path = os.path.join(frame_dir, f"frame_{i:05d}.png")
        cmd.png(png_path, ray=0, quiet=1)

        if i == 0 or (i + 1) % 25 == 0 or i == N_RENDERED_FRAMES - 1:
            deg = i * DEG_PER_FRAME
            print(
                f"[pyr1_lca_biosensor_movie] rendered frame {i + 1}/"
                f"{N_RENDERED_FRAMES} ({deg:.1f} deg of 360)"
            )


def stitch_movie(frame_dir):
    """Stitches the rendered PNG frames into an MP4 with ffmpeg.

    Args:
        frame_dir (str): Directory of frame_%05d.png files (see render_frames).

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
    print("[pyr1_lca_biosensor_movie] running:", " ".join(ffmpeg_cmd))
    result = subprocess.run(ffmpeg_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(
            f"[pyr1_lca_biosensor_movie] ffmpeg failed (exit {result.returncode}); "
            f"frames left in {frame_dir} for debugging."
        )

    size_bytes = os.path.getsize(OUT_MP4)
    duration_s = N_RENDERED_FRAMES / FPS
    print(
        f"[pyr1_lca_biosensor_movie] wrote {OUT_MP4} "
        f"({size_bytes / 1e6:.1f} MB, {N_RENDERED_FRAMES} frames, "
        f"~{duration_s:.1f} s @ {FPS} fps)"
    )

    if not KEEP_FRAMES:
        shutil.rmtree(frame_dir, ignore_errors=True)
        print(f"[pyr1_lca_biosensor_movie] cleaned up temp frame dir {frame_dir}")
    else:
        print(f"[pyr1_lca_biosensor_movie] kept temp frame dir {frame_dir} (PYR1_MOVIE_KEEP_FRAMES=1)")


def main():
    """Runs the full pipeline: load, style, camera, render, stitch."""
    print(
        f"[pyr1_lca_biosensor_movie] mode={'SMOKE TEST' if IS_SMOKE_TEST else 'FULL RUN'}, "
        f"duration={DURATION_S}s, fps={FPS}, full_n_frames={FULL_N_FRAMES}, "
        f"deg_per_frame={DEG_PER_FRAME:.3f}, n_rendered_frames={N_RENDERED_FRAMES}, "
        f"out={OUT_MP4}"
    )

    os.makedirs(OUT_DIR, exist_ok=True)
    setup_global_render_settings()

    load_structures()
    style_scene()
    frame_camera()

    frame_dir = tempfile.mkdtemp(prefix="pyr1_biosensor_turntable_frames_")
    print(f"[pyr1_lca_biosensor_movie] rendering {N_RENDERED_FRAMES} frames into {frame_dir}")
    render_frames(frame_dir)
    stitch_movie(frame_dir)


main()
