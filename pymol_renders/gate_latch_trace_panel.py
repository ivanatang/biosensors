#!/usr/bin/env python
"""Adds a synced distance-trace strip under the gate_latch_closure_ghost_movie.py
frames and stitches the final movie.

Second half of the pipeline described in gate_latch_closure_ghost_movie.py;
run it after every 3D chunk has rendered, with the same GATE_MOVIE key (see
gate_closure_movie_config.py). Needs matplotlib + Pillow, which PyMOL's
bundled python lacks:

    GATE_MOVIE=pair3087_lig3oh ~/miniforge3/bin/python pymol_renders/gate_latch_trace_panel.py

The trace is the movie's trace_xvg (gmx distance on the same two atoms the 3D
dash measures) at full 37.5 ps resolution, so the cursor value and the
in-scene dash label agree.

Each 1600x1000 3D frame is stacked over a 1600x400 strip into a 1600x1400
frame, then ffmpeg-encoded. Set GATE_LATCH_KEEP_FRAMES=1 to keep the
intermediate frame directories.
"""

import glob
import os
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from PIL import Image  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate_closure_movie_config import REPO_ROOT, TIME_SPACING_NS, get_movie  # noqa: E402

MOVIE = get_movie()

OUT_DIR = os.path.join(REPO_ROOT, "pymol_renders", "output")
OUT_NAME = MOVIE["out_name"]
FRAME3D_DIR = os.path.join(OUT_DIR, f".{OUT_NAME}_frames3d")
COMP_DIR = os.path.join(OUT_DIR, f".{OUT_NAME}_framescomp")
OUT_MP4 = os.path.join(OUT_DIR, OUT_NAME + ".mp4")

XVG = MOVIE["trace_xvg"]
T0_NS = MOVIE["t0_ns"]
T_END_NS = MOVIE["t_end_ns"]
STRIDE = int(os.environ.get("GATE_LATCH_STRIDE", "2"))
DIST_LABEL = f"Cα88–{MOVIE['partner_label']}"

FPS = int(os.environ.get("GATE_LATCH_FPS", "10"))
KEEP_FRAMES = os.environ.get("GATE_LATCH_KEEP_FRAMES", "0") == "1"
FFMPEG_BIN = os.environ.get("FFMPEG_BIN", "/opt/homebrew/bin/ffmpeg")

STRIP_W_PX = 1600
STRIP_H_PX = 400
DPI = 300
ROLL_NS = 0.5  # smoothing window for the bold trace; raw trace drawn underneath

TRACE_COLOR = "#648FFF"  # GROUP_COLOR["Binder"]
GATE_COLOR = "#FE6100"
LATCH_COLOR = "#DC267F"
GHOST_COLOR = "#C4C4C4"


def load_xvg(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Reads a two-column GROMACS .xvg, skipping '#' and '@' header lines.

    Args:
        path: Path to the .xvg file.

    Returns:
        Time in ns and distance in nm, as two 1D arrays.
    """
    rows = [
        [float(x) for x in line.split()[:2]]
        for line in open(path)
        if line.strip() and line[0] not in "#@"
    ]
    data = np.array(rows)
    return data[:, 0] / 1000.0, data[:, 1]


def read_ghost_time() -> float | None:
    """Returns the ghost's sim time, written by the 3D script into FRAME3D_DIR.

    Returns None when the movie has no ghost (show_ghost False); falls back
    to the first loaded state's time if the file is missing.
    """
    if not MOVIE["show_ghost"]:
        return None
    try:
        with open(os.path.join(FRAME3D_DIR, "ghost_time_ns.txt")) as f:
            return float(f.read())
    except (OSError, ValueError):
        return MOVIE["open_ns"]


def build_strip(t: np.ndarray, d: np.ndarray, ghost_ns: float | None):
    """Draws the static trace strip and returns handles for the per-frame cursor.

    Args:
        t: Time (ns) within the movie window.
        d: Ca88-to-partner distance (nm) at each time.
        ghost_ns: Time of the open-state ghost, marked with a dashed line and
            a legend entry; None to omit both.

    Returns:
        (fig, cursor_line, cursor_dot, title) for update_cursor().
    """
    fig, ax = plt.subplots(
        figsize=(STRIP_W_PX / DPI, STRIP_H_PX / DPI), dpi=DPI, constrained_layout=True
    )
    n = max(1, int(round(ROLL_NS / TIME_SPACING_NS)))
    smooth = np.convolve(d, np.ones(n) / n, mode="same")
    ax.plot(t, d, color="0.78", lw=0.4)
    ax.plot(t[n:-n], smooth[n:-n], color=TRACE_COLOR, lw=1.2)
    ax.set_xlim(T0_NS, T_END_NS)
    ax.set_ylim(*MOVIE["ylim"])
    ax.set_xlabel("time (ns)", fontsize=6, labelpad=1)
    ax.set_ylabel(f"{DIST_LABEL} (nm)", fontsize=6)
    ax.tick_params(labelsize=5, length=2, pad=1)
    ax.grid(True, alpha=0.4)

    handles = [
        Patch(color=GATE_COLOR, label="gate (84–90)"),
        Patch(color=LATCH_COLOR, label="latch (114–118)"),
    ]
    if ghost_ns is not None:
        handles.append(Patch(color=GHOST_COLOR, label=f"open-state ghost ({ghost_ns:.1f} ns)"))
        ax.axvline(ghost_ns, color=GHOST_COLOR, lw=1.0, ls="--")
    ax.legend(handles=handles, loc="upper right", fontsize=5, frameon=False, handlelength=1.2)

    cursor_line = ax.axvline(T0_NS, color="0.2", lw=0.7)
    (cursor_dot,) = ax.plot([], [], "o", color="0.1", ms=3)
    title = ax.set_title("", loc="left", fontsize=7, fontweight="bold", pad=2)
    return fig, cursor_line, cursor_dot, title


def main():
    """Draws the strip per 3D frame, stacks the two, and encodes the mp4.

    Raises:
        RuntimeError: No 3D frames found, or ffmpeg fails.
    """
    frames = sorted(glob.glob(os.path.join(FRAME3D_DIR, "frame_*.png")))
    if not frames:
        raise RuntimeError(f"no 3D frames in {FRAME3D_DIR}; run gate_latch_closure_ghost_movie.py first")

    t_all, d_all = load_xvg(XVG)
    keep = (t_all >= T0_NS - 1e-6) & (t_all <= T_END_NS)
    t, d = t_all[keep], d_all[keep]

    fig, cursor_line, cursor_dot, title = build_strip(t, d, read_ghost_time())
    os.makedirs(COMP_DIR, exist_ok=True)

    for k, path in enumerate(frames):
        state = int(os.path.basename(path)[6:11]) + 1
        t_ns = T0_NS + (state * STRIDE - 1) * TIME_SPACING_NS
        i = int(np.argmin(np.abs(t - t_ns)))
        cursor_line.set_xdata([t_ns, t_ns])
        cursor_dot.set_data([t[i]], [d[i]])
        title.set_text(f"t = {t_ns:.1f} ns    {DIST_LABEL} = {d[i]:.2f} nm")

        fig.canvas.draw()
        strip = Image.frombuffer("RGBA", fig.canvas.get_width_height(), fig.canvas.buffer_rgba()).convert("RGB")
        top = Image.open(path).convert("RGB")
        if strip.width != top.width:
            strip = strip.resize((top.width, round(strip.height * top.width / strip.width)))
        comp = Image.new("RGB", (top.width, top.height + strip.height), "white")
        comp.paste(top, (0, 0))
        comp.paste(strip, (0, top.height))
        comp.save(os.path.join(COMP_DIR, f"frame_{k:05d}.png"))
        if k % 50 == 0:
            print(f"[gate_latch_trace_panel] composited {k + 1}/{len(frames)} (t={t_ns:.2f} ns)")

    cmd = [
        FFMPEG_BIN, "-y", "-framerate", str(FPS),
        "-i", os.path.join(COMP_DIR, "frame_%05d.png"),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
        "-movflags", "+faststart", OUT_MP4,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"ffmpeg failed (exit {result.returncode}); frames left in {COMP_DIR}")
    print(
        f"[gate_latch_trace_panel] wrote {OUT_MP4} ({os.path.getsize(OUT_MP4) / 1e6:.1f} MB, "
        f"{len(frames)} frames, {len(frames) / FPS:.1f} s @ {FPS} fps)"
    )

    if not KEEP_FRAMES:
        shutil.rmtree(COMP_DIR, ignore_errors=True)
        shutil.rmtree(FRAME3D_DIR, ignore_errors=True)


main()
