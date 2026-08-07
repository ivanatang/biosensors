#!/usr/bin/env python
"""
add_contact_residue_labels.py
--------------------------------
Labels the SPECIFIC gate/latch residues found to be nearest the bridging
water (e.g. "Ala89", "His115") -- not just the "Gate (84-90)"/"Latch
(114-118)" region labels add_water_network_labels.py already adds.

Since gate/latch residues are all colored the same orange/purple, there's
no per-residue color cluster to detect the way add_region_labels.py does
for whole regions. Instead, gate_latch_water_network.py renders a
throwaway MARKER pass with a small unique-cyan pseudoatom placed exactly
at each contact residue's atom position, alongside the real (clean,
marker-free) output pass. This script finds the cyan markers' pixel
positions in the marker pass, then writes plain TEXT-ONLY labels (no
dot/leader line, per the established preference in this project) onto
the real image at those positions.

Usage:
    python add_contact_residue_labels.py
"""
import os
import colorsys
from PIL import Image, ImageDraw, ImageFont

OUT_DIR = os.path.join(os.path.dirname(__file__), "output")
REAL_IMAGE = "pair_3085_binder_water_network.png"
MARKER_IMAGE = "pair_3085_binder_water_network_MARKERS_tmp.png"

# (label text, marker color to find, side to prefer if ambiguous)
MARKERS = [
    ("Ala89", (0, 255, 255)),   # gate contact residue, cyan marker
    ("His115", (0, 255, 255)),  # latch contact residue, cyan marker (same
                                  # color -- distinguished by which cluster
                                  # is which, see main())
]
HUE_TOL_DEG = 10
MIN_SAT = 0.5
MIN_VAL = 0.5
SAMPLE_STEP = 1

FONT_PATH = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_SIZE = 28
LABEL_GAP = 14
SAMPLE_MARGIN = 10
TEXT_COLOR = (30, 30, 30)


def rgb_to_hue_deg(r, g, b):
    h, s, v = colorsys.rgb_to_hsv(r / 255, g / 255, b / 255)
    return h * 360, s, v


def find_marker_clusters(img, target_rgb):
    """Cyan markers are small and there are two of them -- find ALL
    matching pixels, then split into up-to-two spatial clusters by
    simple x-position gap (gate marker and latch marker are far apart
    in this figure's layout), rather than assuming a single blob."""
    target_hue, _, _ = rgb_to_hue_deg(*target_rgb)
    w, h = img.size
    px = img.load()
    pts = []
    for y in range(0, h, SAMPLE_STEP):
        for x in range(0, w, SAMPLE_STEP):
            r, g, b = px[x, y][:3]
            hue, sat, val = rgb_to_hue_deg(r, g, b)
            if sat < MIN_SAT or val < MIN_VAL:
                continue
            dh = min(abs(hue - target_hue), 360 - abs(hue - target_hue))
            if dh <= HUE_TOL_DEG:
                pts.append((x, y))
    if not pts:
        return []

    pts.sort()
    clusters = [[pts[0]]]
    for p in pts[1:]:
        if p[0] - clusters[-1][-1][0] > 30:  # gap in x -> new cluster
            clusters.append([p])
        else:
            clusters[-1].append(p)

    out = []
    for c in clusters:
        xs = [p[0] for p in c]
        ys = [p[1] for p in c]
        out.append(dict(x_center=sum(xs) // len(xs), y_center=sum(ys) // len(ys), n=len(c)))
    return out


def main():
    marker_path = os.path.join(OUT_DIR, MARKER_IMAGE)
    real_path = os.path.join(OUT_DIR, REAL_IMAGE)

    marker_img = Image.open(marker_path).convert("RGB")
    clusters = find_marker_clusters(marker_img, (0, 255, 255))
    print(f"found {len(clusters)} cyan marker cluster(s): {clusters}")

    if len(clusters) != 2:
        raise RuntimeError(
            f"expected exactly 2 marker clusters (gate + latch), found "
            f"{len(clusters)}. Inspect {marker_path} manually."
        )

    # Left-most cluster = latch (purple side, left in current framing),
    # right-most = gate (orange side, right) -- matches the current
    # TURN_Y_DEG=-25 framing. Determine by x position rather than
    # hardcoding which index came first from the scan.
    clusters.sort(key=lambda c: c["x_center"])
    latch_marker, gate_marker = clusters[0], clusters[1]

    img = Image.open(real_path).convert("RGB")
    w, h = img.size
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

    for label, marker in [("Ala89", gate_marker), ("His115", latch_marker)]:
        bbox = draw.textbbox((0, 0), label, font=font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        # Place just below-right (gate) or below-left (latch) of the marker
        # point itself, offset enough to clear the sticks there, text-only.
        if marker is gate_marker:
            text_x = min(w - tw - SAMPLE_MARGIN, marker["x_center"] + LABEL_GAP)
        else:
            text_x = max(SAMPLE_MARGIN, marker["x_center"] - LABEL_GAP - tw)
        text_y = min(h - th - SAMPLE_MARGIN, marker["y_center"] + LABEL_GAP)
        draw.text((text_x, text_y), label, font=font, fill=TEXT_COLOR)
        print(f"placed '{label}' at ({text_x},{text_y}) near marker {marker}")

    img.save(real_path)
    print(f"[add_contact_residue_labels] labeled {real_path}")

    os.remove(marker_path)
    print(f"[add_contact_residue_labels] cleaned up temp marker pass {marker_path}")


if __name__ == "__main__":
    main()
