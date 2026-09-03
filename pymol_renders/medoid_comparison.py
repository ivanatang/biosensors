#!/usr/bin/env python
"""Renders side-by-side binder vs nonbinder medoid comparison figures.

Publication-quality figures comparing a PYR1+LCA biosensor "binder" medoid
structure against a "nonbinder" medoid structure, for a research talk.

Mechanistic story being illustrated: in binders, the gate loop (resi
84-90) and latch (resi 114-118) close down over the bound ligand (LIG); in
nonbinders/false positives this gate-latch closure is disrupted even
though the ligand may still occupy the pocket.

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/biosensors/pymol_renders/medoid_comparison.py

Produces two ray-traced PNGs (one per structure, not a single composited
side-by-side image) in pymol_renders/output/:

    binder_pair_3059_medoid.png
    nonbinder_pair_0052_medoid.png

Shared camera (important): both structures load into the same PyMOL
session. The nonbinder is structurally aligned onto the binder first
(cmd.align, CA atoms of the stable scaffold, excluding the gate/latch resi
ranges from the fit so the alignment doesn't force the very loops we're
comparing to superimpose). One camera view (orientation + zoom) is then
derived from the binder's gate+latch+ligand region and applied identically
(cmd.set_view with the same matrix) to both renders, so a side-by-side
comparison of gate/latch closure reflects a real conformational
difference, not an artifact of two independently auto-oriented cameras.

Background is plain white in both panels (talk-slide friendly, no group
tint) -- panel identity (binder vs. nonbinder) is conveyed instead by (1)
filename (binder_... / nonbinder_...) and (2) an in-scene text label (e.g.
"Binder (pair_3059)") placed at one shared anchor point computed from both
structures' pocket geometry together (average ligand centroid, offset
outward from the average protein centroid by the larger of the two
ligands' radius plus a margin). Using one shared anchor for both panels
means the same "is this position clear of the structure" bounding logic
applies to both, rather than each panel getting its own independently
computed (and independently risky) position.

Gate (orange, #FE6100) and latch (purple, #785EF0) use one consistent
color each across both panels, since gate/latch are structural landmarks,
not themselves a binder/nonbinder label. Binder-vs-nonbinder is conveyed
only via in-scene label/filename, as above.

Gate/latch labels: this script only renders the colored structure; gate
and latch region labels are added afterward, in true pixel space, by
add_region_labels.py (run it on the PNGs this script produces). An earlier
version tried to draw the legend in-scene via CGO squares anchored on
cmd.get_view()'s decoded camera axes; that approach couldn't be verified
without executing PyMOL (no way to preview the render) and turned out not
to render at all, so it was replaced with the pixel-space post-processing
script instead, which is directly checkable.

To re-run on a different sequence pair, edit the STRUCTURES list below.
The alignment/shared-view logic assumes exactly two structures named
"binder" and "nonbinder" (via each dict's "group" key); adapt
align_structures()/compute_shared_view() to add more panels.
"""

import os

from pymol import cmd, util

# --------------------------------------------------------------------------
# Config - edit these to re-point the script at a different pair
# --------------------------------------------------------------------------
REPO_ROOT = "/Users/ivanatang/Developer/biosensors"
OUT_DIR = os.path.join(REPO_ROOT, "pymol_renders", "output")

STRUCTURES = [
    {
        "label": "Binder (pair_3059)",
        "group": "binder",
        "pdb": os.path.join(
            REPO_ROOT,
            "binders/pair_3059_binder/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_PL.pdb",
        ),
        "out_png": os.path.join(OUT_DIR, "binder_pair_3059_medoid.png"),
        "accent_hex": "#648FFF",  # project "Binder" color, CLAUDE.md palette
    },
    {
        "label": "Nonbinder (pair_0052)",
        "group": "nonbinder",
        "pdb": os.path.join(
            REPO_ROOT,
            "nonbinders/pair_0052_nb/prod_md_0p9_cutoff_3dt_64x1_16PME_642dd/medoid_PL.pdb",
        ),
        "out_png": os.path.join(OUT_DIR, "nonbinder_pair_0052_medoid.png"),
        "accent_hex": "#DC267F",  # project "False Positive" color, CLAUDE.md palette
    },
]

# The reference structure that the other structure(s) get aligned onto,
# and whose pocket region is used to derive the single shared camera view.
REFERENCE_GROUP = "binder"

# Structural landmarks (see CLAUDE.md "Key domain conventions")
GATE_RESI = "84-90"
LATCH_RESI = "114-118"
PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"
LIGAND_RESN = "LIG"

# Consistent landmark colors across BOTH panels (do not use these for
# binder-vs-nonbinder distinction, see module docstring)
GATE_COLOR_HEX = "#FE6100"   # orange
LATCH_COLOR_HEX = "#785EF0"  # purple
PROTEIN_GRAY = "gray85"

IMG_WIDTH = 2000
IMG_HEIGHT = 1600
PNG_DPI = 300

# Camera framing knobs
ZOOM_BUFFER = 5.0          # padding (Angstrom) around the framed selection
LABEL_MARGIN = 10.0        # Angstrom beyond the ligand radius for the shared label anchor


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


def _max_radius(model, center):
    """Returns the max distance from `center` to any atom in `model`.

    Args:
        model: PyMOL chempy Indexed/Storable model (has a .atom list).
        center (tuple[float, float, float]): Reference point.

    Returns:
        float: Maximum atom-to-center distance (Angstrom).
    """
    return max(
        (
            (a.coord[0] - center[0]) ** 2
            + (a.coord[1] - center[1]) ** 2
            + (a.coord[2] - center[2]) ** 2
        )
        ** 0.5
        for a in model.atom
    )


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


def load_and_style(struct):
    """Loads one structure and builds its cartoon/sticks styling.

    Loads into an object named after its group ("binder"/"nonbinder") and
    builds per-object-prefixed named selections so both structures can
    coexist in the same session without name collisions.

    Args:
        struct (dict): One entry of STRUCTURES.
    """
    obj = struct["group"]
    cmd.load(struct["pdb"], obj)
    cmd.remove(f"{obj} and hydro")

    protein_sel = f"{obj} and chain {PROTEIN_CHAIN} and polymer"
    gate_sel = f"{protein_sel} and resi {GATE_RESI}"
    latch_sel = f"{protein_sel} and resi {LATCH_RESI}"
    ligand_sel = f"{obj} and chain {LIGAND_CHAIN} and resn {LIGAND_RESN}"

    cmd.select(f"{obj}_gate", gate_sel)
    cmd.select(f"{obj}_latch", latch_sel)
    cmd.select(f"{obj}_ligand", ligand_sel)
    cmd.select(
        f"{obj}_landmark", f"{obj}_gate or {obj}_latch or {obj}_ligand"
    )

    # base protein cartoon
    cmd.hide("everything", obj)
    cmd.show("cartoon", protein_sel)
    cmd.color(PROTEIN_GRAY, protein_sel)

    # gate / latch: one consistent color each, thicker loop tube +
    # sidechain sticks so closure geometry (not just a flat tube) is
    # visible
    cmd.set_color("gate_color", list(hex_to_rgb01(GATE_COLOR_HEX)))
    cmd.set_color("latch_color", list(hex_to_rgb01(LATCH_COLOR_HEX)))

    cmd.color("gate_color", f"{obj}_gate")
    cmd.color("latch_color", f"{obj}_latch")

    cmd.set("cartoon_loop_radius", 0.4, f"{obj}_gate")
    cmd.set("cartoon_loop_radius", 0.4, f"{obj}_latch")
    cmd.set("cartoon_tube_radius", 0.4, f"{obj}_gate")
    cmd.set("cartoon_tube_radius", 0.4, f"{obj}_latch")

    cmd.show(
        "sticks", f"({gate_sel} or {latch_sel}) and not name C+N+O"
    )
    cmd.set("stick_radius", 0.18, f"{obj}_gate or {obj}_latch")
    cmd.color("gate_color", f"({gate_sel}) and elem C")
    cmd.color("latch_color", f"({latch_sel}) and elem C")

    # ligand: sticks, yellow-carbon scheme so it reads clearly against the
    # gray protein / orange gate / purple latch
    cmd.show("sticks", ligand_sel)
    cmd.set("stick_radius", 0.22, ligand_sel)
    util.cbay(ligand_sel)

    cmd.deselect()


def align_structures():
    """Aligns every non-reference structure onto REFERENCE_GROUP.

    Uses CA atoms of the stable scaffold, excluding the gate/latch resi
    ranges from the fit -- deliberate, since those are exactly the regions
    being compared, so forcing them into the fit would bias the "closure"
    comparison toward looking identical.
    """
    exclude = f"(resi {GATE_RESI} or resi {LATCH_RESI})"
    ref = REFERENCE_GROUP
    for struct in STRUCTURES:
        obj = struct["group"]
        if obj == ref:
            continue
        mobile_sel = (
            f"{obj} and chain {PROTEIN_CHAIN} and polymer and name CA and not {exclude}"
        )
        target_sel = (
            f"{ref} and chain {PROTEIN_CHAIN} and polymer and name CA and not {exclude}"
        )
        result = cmd.align(mobile_sel, target_sel, cycles=5)
        print(
            f"[medoid_comparison] aligned {obj!r} onto {ref!r} "
            f"(scaffold CA, gate/latch excluded from fit): "
            f"RMSD={result[0]:.3f} A over {result[1]} atoms "
            f"(raw RMSD before outlier rejection={result[3]:.3f} A over {result[4]} atoms)"
        )


def compute_shared_label_anchor():
    """Computes one label anchor position, shared by both panels.

    Computed from both structures' ligand + protein geometry together
    (post-alignment). Because it's the same 3D point used for both
    renders, and computed from the union of both structures, it's
    guaranteed to clear both ligands' immediate footprint by construction,
    not just one panel's.

    Returns:
        tuple[float, float, float] | None: The shared anchor position, or
        None if either selection matched no atoms.
    """
    protein_sel = " or ".join(
        f"({s['group']} and chain {PROTEIN_CHAIN} and polymer and name CA)"
        for s in STRUCTURES
    )
    ligand_sel = " or ".join(f"{s['group']}_ligand" for s in STRUCTURES)

    protein_model = cmd.get_model(protein_sel)
    ligand_model = cmd.get_model(ligand_sel)
    if not protein_model.atom or not ligand_model.atom:
        return None

    protein_c = _centroid(protein_model)
    ligand_c = _centroid(ligand_model)
    outward = _norm(_sub(ligand_c, protein_c))
    ligand_radius = _max_radius(ligand_model, ligand_c)

    return _add(ligand_c, outward, ligand_radius + LABEL_MARGIN)


def add_group_label(struct, position):
    """Adds a colored text label pseudoatom for one structure at a fixed position.

    Args:
        struct (dict): One entry of STRUCTURES.
        position (tuple[float, float, float]): Label anchor position, from
            compute_shared_label_anchor.

    Returns:
        str: Name of the created pseudoatom object.
    """
    obj = struct["group"]
    anchor = f"{obj}_label_anchor"
    cmd.pseudoatom(anchor, pos=list(position))
    cmd.hide("everything", anchor)

    color_name = f"{obj}_accent_color"
    cmd.set_color(color_name, list(hex_to_rgb01(struct["accent_hex"])))
    cmd.set("label_size", 30, anchor)
    cmd.set("label_color", color_name, anchor)
    cmd.set("label_font_id", 7, anchor)
    cmd.label(anchor, '"%s"' % struct["label"])
    return anchor


def compute_shared_view(anchor_names):
    """Computes the one camera view (orient + zoom) shared by both panels.

    Orients using the reference structure's gate+latch+ligand region (a
    stable, reproducible reference direction), then zooms tightly on the
    union of both structures' gate+latch+ligand atoms plus both label
    anchors, so the resulting camera keeps everything relevant in-frame
    (no ligand cropping) without excess dead canvas.

    Args:
        anchor_names (list[str]): Label pseudoatom names from
            add_group_label, included in the zoom selection.

    Returns:
        The PyMOL view matrix from cmd.get_view(), for cmd.set_view().
    """
    ref_landmark = f"{REFERENCE_GROUP}_landmark"
    cmd.orient(ref_landmark)

    combined_landmark = " or ".join(f"{s['group']}_landmark" for s in STRUCTURES)
    zoom_sel = combined_landmark + " or " + " or ".join(anchor_names)
    cmd.zoom(zoom_sel, buffer=ZOOM_BUFFER)

    return cmd.get_view()


def render_panel(struct, shared_view, group_to_obj_names):
    """Renders exactly one structure + its own label, using the shared camera.

    Every other structure/label object is disabled first.

    Args:
        struct (dict): One entry of STRUCTURES; the panel to render.
        shared_view: View matrix from compute_shared_view, applied
            identically to every panel.
        group_to_obj_names (dict): Maps group name to the tuple of PyMOL
            object/pseudoatom names to enable for that panel.
    """
    this_group = struct["group"]
    for group, names in group_to_obj_names.items():
        enable = group == this_group
        for name in names:
            (cmd.enable if enable else cmd.disable)(name)

    # Re-applied explicitly right before each render as a safety net:
    # enable/disable shouldn't perturb the view, but this guarantees both
    # panels are bit-for-bit the same camera.
    cmd.set_view(shared_view)

    cmd.bg_color("white")

    cmd.ray(IMG_WIDTH, IMG_HEIGHT)
    cmd.png(struct["out_png"], dpi=PNG_DPI, ray=0)
    print(f"[medoid_comparison] wrote {struct['out_png']}")


def main():
    """Runs the full pipeline: load, align, label, frame camera, render both panels."""
    os.makedirs(OUT_DIR, exist_ok=True)
    setup_global_render_settings()

    for struct in STRUCTURES:
        load_and_style(struct)

    align_structures()

    label_pos = compute_shared_label_anchor()
    anchor_by_group = {}
    for struct in STRUCTURES:
        anchor_by_group[struct["group"]] = add_group_label(struct, label_pos)

    shared_view = compute_shared_view(list(anchor_by_group.values()))

    group_to_obj_names = {
        s["group"]: (s["group"], anchor_by_group[s["group"]]) for s in STRUCTURES
    }

    for struct in STRUCTURES:
        render_panel(struct, shared_view, group_to_obj_names)


main()
