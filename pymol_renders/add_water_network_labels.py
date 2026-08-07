#!/usr/bin/env python
"""
add_water_network_labels.py
-----------------------------
Same pixel-space HSV-cluster-detection labeling approach as
add_region_labels.py / add_residue_labels_5.py, applied to the
gate-latch-ligand water network figure. Labels Gate and Latch only
(the explicit ask); the ligand's water-interacting oxygen and the
bridging water are already visually distinguished by color in the
render itself and don't need a text label to be legible.

Usage:
    python add_water_network_labels.py
"""
import os
import colorsys
from PIL import Image, ImageDraw, ImageFont

OUT_DIR = os.path.join(os.path.dirname(__file__), "output")
IMAGE = "pair_3085_binder_water_network.png"

REGIONS = [
    ("Gate (84-90)",    (254, 97, 0)),    # #FE6100 orange
    ("Latch (114-118)", (120, 94, 240)),  # #785EF0 purple
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
    h, s, v = colorsys.rgb_to_hsv(r / 255, g / 255, b / 255)
    return h * 360, s, v


def find_cluster(img, target_rgb):
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
    )


def main():
    path = os.path.join(OUT_DIR, IMAGE)
    img = Image.open(path).convert("RGB")
    w, h = img.size
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

    for label, rgb in REGIONS:
        c = find_cluster(img, rgb)
        if c is None:
            print(f"  ! {label}: no matching pixels found, skipping")
            continue
        print(f"{label}: {c}")

        bbox = draw.textbbox((0, 0), label, font=font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]

        if c["x_center"] <= w / 2:
            text_x = max(SAMPLE_MARGIN, c["x_min"] - LABEL_GAP - tw)
        else:
            text_x = min(w - tw - SAMPLE_MARGIN, c["x_max"] + LABEL_GAP)
        text_y = max(SAMPLE_MARGIN, min(h - th - SAMPLE_MARGIN, c["y_center"] - th // 2))

        draw.text((text_x, text_y), label, font=font, fill=rgb)

    img.save(path)
    print(f"\n[add_water_network_labels] labeled {path}")


if __name__ == "__main__":
    main()
