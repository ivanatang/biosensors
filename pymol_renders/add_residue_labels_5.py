#!/usr/bin/env python
"""Labels an arbitrary set of single-residue color highlights.

Generalizes add_region_labels.py's method (HSV hue clustering + plain text
labels placed in the nearest clear space, computed from the actual image)
to an arbitrary set of residues, here 5 on seq10_binder.

Usage:
    python add_residue_labels_5.py
"""
import os
import colorsys
from PIL import Image, ImageDraw, ImageFont

OUT_DIR = os.path.join(os.path.dirname(__file__), "output")
IMAGE = "seq10_binder_5residue_highlight.png"

RESIDUES = [
    # (label, RGB) -- label is the amino acid actually present at this
    # position in seq10_binder (confirmed against both the PDB and
    # feat_table_500ns.xlsx's Sequence column, 1-indexed match)
    ("Phe61",  (230, 159, 0)),   # #E69F00 orange
    ("Ile62",  (86, 180, 233)),  # #56B4E9 sky blue
    ("Asp81",  (0, 158, 115)),   # #009E73 bluish green
    ("Arg116", (178, 24, 43)),   # #B2182B crimson
    ("Val164", (204, 121, 167)), # #CC79A7 reddish purple
]
HUE_TOL_DEG = 12
MIN_SAT = 0.30
MIN_VAL = 0.20
SAMPLE_STEP = 2

FONT_PATH = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_SIZE = 32
LABEL_GAP = 16
SAMPLE_MARGIN = 10


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


def find_cluster(img, target_rgb):
    """Finds the pixel cluster matching a target hue.

    Args:
        img: PIL Image (RGB) to scan.
        target_rgb (tuple[int, int, int]): RGB color to match by hue.

    Returns:
        dict | None: Bounding box, center, and pixel count (x_min, x_max,
        y_min, y_max, x_center, y_center, n) of matching pixels, or None
        if no pixels matched.
    """
    target_hue, _, _ = rgb_to_hue_deg(*target_rgb)
    w, h = img.size
    px = img.load()
    xs, ys = [], []
    for y in range(0, h, SAMPLE_STEP):
        for x in range(0, w, SAMPLE_STEP):
            r, g, b = px[x, y][:3]
            hue, sat, val = rgb_to_hue_deg(r, g, b)
            if sat < MIN_SAT or val < MIN_VAL:
                continue
            dh = min(abs(hue - target_hue), 360 - abs(hue - target_hue))
            if dh <= HUE_TOL_DEG:
                xs.append(x)
                ys.append(y)
    if not xs:
        return None
    return dict(
        x_min=min(xs), x_max=max(xs), y_min=min(ys), y_max=max(ys),
        x_center=sum(xs) // len(xs), y_center=sum(ys) // len(ys),
        n=len(xs),
    )


def main():
    """Labels all 5 residues in RESIDUES onto IMAGE, in place."""
    path = os.path.join(OUT_DIR, IMAGE)
    img = Image.open(path).convert("RGB")
    w, h = img.size
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

    clusters = {}
    for label, rgb in RESIDUES:
        c = find_cluster(img, rgb)
        clusters[label] = c
        print(f"resi {label}: {c}")

    for label, rgb in RESIDUES:
        c = clusters[label]
        if c is None:
            print(f"  ! resi {label}: no matching pixels found, skipping label")
            continue
        text = label
        bbox = draw.textbbox((0, 0), text, font=font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]

        # Decide left/right placement based on which side of the image
        # midline the cluster sits on, so the label goes toward open
        # margin rather than back into the structure.
        if c["x_center"] <= w / 2:
            text_x = max(SAMPLE_MARGIN, c["x_min"] - LABEL_GAP - tw)
        else:
            text_x = min(w - tw - SAMPLE_MARGIN, c["x_max"] + LABEL_GAP)
        text_y = max(SAMPLE_MARGIN, min(h - th - SAMPLE_MARGIN, c["y_center"] - th // 2))

        draw.text((text_x, text_y), text, font=font, fill=tuple(rgb))

    img.save(path)
    print(f"\n[add_residue_labels_5] labeled {path}")


if __name__ == "__main__":
    main()
