#!/usr/bin/env python
"""Labels the specific gate/latch residues nearest the bridging water.

Labels individual residues (e.g. "Ala89", "His115"), not just the whole
"Gate (84-90)"/"Latch (114-118)" regions add_water_network_labels.py
already labels.

Gate/latch residues are all colored the same orange/purple, so there's no
per-residue color cluster to detect the way add_region_labels.py does for
whole regions. Instead, gate_latch_water_network.py renders a throwaway
marker pass with a small unique-cyan pseudoatom placed exactly at each
contact residue's atom position, alongside the real (clean, marker-free)
output pass. This script finds the cyan markers' pixel positions in the
marker pass, then writes plain text-only labels (no dot/leader line, per
this project's established preference) onto the real image at those
positions.

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
    ("His115", (0, 255, 255)),  # latch contact residue, same cyan marker
                                 # color -- distinguished by cluster, see main()
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
    """Converts 0-255 RGB to (hue in degrees, saturation, value).

    Args:
        r (int): Red channel, 0-255.
        g (int): Green channel, 0-255.
        b (int): Blue channel, 0-255.

    Returns:
        tuple[float, float, float]: (hue_deg, saturation, value), each in
        their native HSV ranges (hue in [0, 360), sat/value in [0, 1]).
    """
    h, s, v = colorsys.rgb_to_hsv(r / 255, g / 255, b / 255)
    return h * 360, s, v


def find_marker_clusters(img, target_rgb):
    """Finds spatial clusters of pixels matching a target hue.

    Matches all pixels within HUE_TOL_DEG of `target_rgb`'s hue (subject to
    MIN_SAT/MIN_VAL), then splits them into clusters by x-position gaps
    rather than assuming a single blob, since the two cyan contact markers
    (gate and latch) sit far apart in this figure's layout.

    Args:
        img: PIL Image (RGB) to scan.
        target_rgb (tuple[int, int, int]): RGB color to match by hue.

    Returns:
        list[dict]: One dict per cluster (x_center, y_center, n), in the
        order clusters were encountered while scanning top-to-bottom.
    """
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
    """Labels Ala89/His115 on the real image using the marker pass, in place."""
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
