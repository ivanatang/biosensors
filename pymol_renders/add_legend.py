#!/usr/bin/env python
"""
add_legend.py
-------------
Post-process the rendered medoid comparison PNGs to add a gate/latch color
legend via true pixel-space compositing (PIL), rather than PyMOL in-scene
3D labels. Two prior attempts at camera-matrix-derived label/CGO placement
either overlapped the structure or failed to render at all (get_view() row
convention couldn't be verified without actually running PyMOL). This
sidesteps that entirely by drawing directly onto the finished PNG.

Usage:
    python add_legend.py
"""
import os
from PIL import Image, ImageDraw, ImageFont

OUT_DIR = os.path.join(os.path.dirname(__file__), "output")
IMAGES = [
    "binder_pair_3059_medoid.png",
    "nonbinder_pair_0052_medoid.png",
]

GATE_COLOR = (254, 97, 0)      # #FE6100
LATCH_COLOR = (120, 94, 240)   # #785EF0
TEXT_COLOR = (40, 40, 40)

FONT_PATH = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_SIZE = 34
SWATCH = 30
PAD = 18
LINE_GAP = 14
MARGIN = 40


def add_legend(path):
    img = Image.open(path).convert("RGBA")
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

    entries = [("Gate (84-90)", GATE_COLOR), ("Latch (114-118)", LATCH_COLOR)]

    # Measure text extents to size the legend box.
    text_widths = []
    text_heights = []
    for label, _ in entries:
        bbox = draw.textbbox((0, 0), label, font=font)
        text_widths.append(bbox[2] - bbox[0])
        text_heights.append(bbox[3] - bbox[1])
    row_h = max(SWATCH, max(text_heights)) + LINE_GAP
    box_w = PAD * 3 + SWATCH + max(text_widths)
    box_h = PAD * 2 + row_h * len(entries) - LINE_GAP

    # Bottom-left corner, above the bottom margin (empty background there
    # in both renders, confirmed by inspection).
    x0 = MARGIN
    y0 = img.height - MARGIN - box_h

    draw.rounded_rectangle(
        [x0, y0, x0 + box_w, y0 + box_h],
        radius=12,
        fill=(255, 255, 255, 235),
        outline=(120, 120, 120, 255),
        width=2,
    )

    y = y0 + PAD
    for label, color in entries:
        sw_y = y + (row_h - LINE_GAP - SWATCH) // 2
        draw.rounded_rectangle(
            [x0 + PAD, sw_y, x0 + PAD + SWATCH, sw_y + SWATCH],
            radius=6,
            fill=color,
        )
        text_y = y + (row_h - LINE_GAP - text_heights[0]) // 2 - 2
        draw.text((x0 + PAD * 2 + SWATCH, text_y), label, font=font, fill=TEXT_COLOR)
        y += row_h

    img.convert("RGB").save(path)
    print(f"[add_legend] legend added to {path} (box {box_w}x{box_h} at ({x0},{y0}))")


def main():
    for name in IMAGES:
        add_legend(os.path.join(OUT_DIR, name))


if __name__ == "__main__":
    main()
