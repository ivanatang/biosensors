#!/usr/bin/env python
"""Labels the gate and latch regions directly next to where they render.

Labels next to the actual rendered geometry instead of a corner legend.
Finds each region's stick-color pixel cluster by HSV hue matching (robust
to ray-traced shading/anti-aliasing, unlike exact RGB matching), then
places a plain text label (no pointer line or marker) in the nearest clear
side margin, vertically centered on the cluster. Pure PIL pixel-space work,
so no PyMOL camera-matrix guessing is needed.

Usage:
    python add_region_labels.py
"""
import os
import colorsys
from PIL import Image, ImageDraw, ImageFont

OUT_DIR = os.path.join(os.path.dirname(__file__), "output")
IMAGES = [
    "binder_pair_3059_medoid.png",
    "nonbinder_pair_0052_medoid.png",
]

GATE_RGB = (254, 97, 0)      # #FE6100
LATCH_RGB = (120, 94, 240)   # #785EF0
HUE_TOL_DEG = 18
MIN_SAT = 0.35
MIN_VAL = 0.25

FONT_PATH = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_SIZE = 34
TEXT_COLOR_GATE = (170, 60, 0)
TEXT_COLOR_LATCH = (75, 55, 165)
LABEL_GAP = 18     # gap between cluster edge and label text
SAMPLE_STEP = 3    # subsample pixels for speed


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
        dict: Bounding box and center (x_min, x_max, y_min, y_max,
        x_center, y_center) of matching pixels.

    Raises:
        RuntimeError: No pixels matched the target hue.
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
        raise RuntimeError(f"No pixels matched target hue {target_hue:.1f} deg")
    return dict(
        x_min=min(xs), x_max=max(xs), y_min=min(ys), y_max=max(ys),
        x_center=sum(xs) // len(xs), y_center=sum(ys) // len(ys),
    )


def label_image(path):
    """Labels the gate and latch regions on one PNG, in place.

    Args:
        path (str): Path to the PNG to modify and overwrite.
    """
    img = Image.open(path).convert("RGB")
    w, h = img.size
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

    gate = find_cluster(img, GATE_RGB)
    latch = find_cluster(img, LATCH_RGB)

    # Decided per-image from actual cluster centers, not assumed layout, in
    # case a future re-render flips gate/latch left-right.
    gate_is_left = gate["x_center"] <= latch["x_center"]

    def place(label_text, cluster, color, side):
        """Draws one label at the clear-margin side of its cluster."""
        bbox = draw.textbbox((0, 0), label_text, font=font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        anchor_y = cluster["y_center"]
        if side == "left":
            anchor_x = cluster["x_min"]
            text_x = max(10, anchor_x - LABEL_GAP - tw)
        else:
            anchor_x = cluster["x_max"]
            text_x = min(w - tw - 10, anchor_x + LABEL_GAP)
        text_y = anchor_y - th // 2
        text_y = max(10, min(h - th - 10, text_y))
        draw.text((text_x, text_y), label_text, font=font, fill=color)

    if gate_is_left:
        place("Gate (84-90)", gate, TEXT_COLOR_GATE, "left")
        place("Latch (114-118)", latch, TEXT_COLOR_LATCH, "right")
    else:
        place("Gate (84-90)", gate, TEXT_COLOR_GATE, "right")
        place("Latch (114-118)", latch, TEXT_COLOR_LATCH, "left")

    img.save(path)
    print(f"[add_region_labels] labeled {path}  gate={gate}  latch={latch}")


def main():
    """Labels every PNG in IMAGES."""
    for name in IMAGES:
        label_image(os.path.join(OUT_DIR, name))


if __name__ == "__main__":
    main()
