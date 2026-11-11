"""Redraw the row and column labels of the histology composites at a legible print size.

Figure 3E/F: column labels (Healthy, DD, GG, GD) above each block and tissue labels, rotated and in two lines,
to the left of each block. Supplementary Figure 6A/6B: column labels ("DD, 1-year", ...) and rotated tissue labels.
Only the label margins are cleared and redrawn; micrographs, scale bars, arrows and the panel letters E and F are
left untouched.

usage (run from the repository root, after draw_histology_arrows.py):
  python3 relabel_histology_panels.py SCALEBAR_SPEC.json IN_DIR OUT_DIR
IN_DIR holds histology_panels_EF_arrows.png, SupplFig6A_arrows.png and SupplFig6B_arrows.png;
OUT_DIR receives histology_panels_EF.png, SupplFig6A.png and SupplFig6B.png.
"""
import json
import os
import sys
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None
spec_path, in_dir, out_dir = sys.argv[1:4]
fig3_job, s6a_job, s6b_job = json.load(open(spec_path, encoding="utf-8"))["jobs"]
os.makedirs(out_dir, exist_ok=True)
FONT_DIR = os.path.dirname(fig3_job["font"])          # Noto Sans directory given in the scale-bar spec
BOLD = os.path.join(FONT_DIR, "NotoSans-Bold.ttf")
REG = os.path.join(FONT_DIR, "NotoSans-Regular.ttf")
TISSUES = ["Foot muscle", "Digestive gland", "Gill", "Mantle"]


def centred_text(d, box, text, font):
    x0, y0, x1, y1 = box
    b = d.textbbox((0, 0), text, font=font)
    d.text(((x0 + x1) / 2 - (b[0] + b[2]) / 2, (y0 + y1) / 2 - (b[1] + b[3]) / 2), text, fill=(0, 0, 0), font=font)


def rotated_label(img, centre, text, font, spacing, max_w, max_h):
    probe = ImageDraw.Draw(Image.new("L", (10, 10)))
    b = probe.multiline_textbbox((0, 0), text, font=font, spacing=spacing, align="center")
    tw, th = int(round(b[2] - b[0])), int(round(b[3] - b[1]))
    assert th + 8 <= max_w and tw + 8 <= max_h, (text, tw, th, max_w, max_h)
    layer = Image.new("L", (tw + 20, th + 20), 0)
    ImageDraw.Draw(layer).multiline_text((10 - b[0], 10 - b[1]), text, fill=255, font=font, spacing=spacing, align="center")
    layer = layer.rotate(90, expand=True)
    x, y = round(centre[0] - layer.width / 2), round(centre[1] - layer.height / 2)
    img.paste(Image.new("RGB", layer.size, (0, 0, 0)), (x, y), layer)


def grid(tiles, nrow, ncol, first=0):
    cols = [(tiles[first + c][0], tiles[first + c][2]) for c in range(ncol)]
    rows = [(min(tiles[first + r * ncol + c][1] for c in range(ncol)), max(tiles[first + r * ncol + c][3] for c in range(ncol)))
            for r in range(nrow)]
    return cols, rows


# ---------------------------------------------------------------- Figure 3E/F
img = Image.open(os.path.join(in_dir, "histology_panels_EF_arrows.png")).convert("RGB")
d = ImageDraw.Draw(img)
col_font, row_font = ImageFont.truetype(BOLD, 58), ImageFont.truetype(BOLD, 50)
for first in (0, 16):                                    # E block, then F block
    cols, rows = grid(fig3_job["tiles"], 4, 4, first)
    top = rows[0][0]
    left = 0 if first == 0 else fig3_job["tiles"][15][2] + 1   # F labels sit between the two blocks
    d.rectangle([cols[0][0], 0, cols[-1][1], top - 1], fill=(255, 255, 255))
    d.rectangle([left, 140, cols[0][0] - 1, img.height - 1], fill=(255, 255, 255))   # letters E/F end above y=140
    for (x0, x1), name in zip(cols, ["Healthy", "DD", "GG", "GD"]):
        centred_text(d, (x0, 0, x1, top - 1), name, col_font)
    for (y0, y1), tis in zip(rows, TISSUES):
        rotated_label(img, ((left + cols[0][0]) / 2, (y0 + y1) / 2), tis.replace(" ", "\n"), row_font, 4,
                      cols[0][0] - left, y1 - y0)
img.save(os.path.join(out_dir, "histology_panels_EF.png"))
print("Figure 3E/F labels redrawn", img.size)

# ---------------------------------------------------------------- Supplementary Figure 6A/6B
for job, name, age in ((s6a_job, "SupplFig6A", "1-year"), (s6b_job, "SupplFig6B", "2-year")):
    img = Image.open(os.path.join(in_dir, f"{name}_arrows.png")).convert("RGB")
    d = ImageDraw.Draw(img)
    cols, rows = grid(job["tiles"], 4, 3)
    top = min(r[0] for r in rows)
    d.rectangle([cols[0][0], 0, cols[-1][1], top - 1], fill=(255, 255, 255))
    d.rectangle([0, top, cols[0][0] - 1, img.height - 1], fill=(255, 255, 255))
    col_font, row_font = ImageFont.truetype(REG, 62), ImageFont.truetype(REG, 60)
    for (x0, x1), g in zip(cols, ["DD", "GG", "GD"]):
        centred_text(d, (x0, top - 120, x1, top - 10), f"{g}, {age}", col_font)
    for (y0, y1), tis in zip(rows, TISSUES):
        rotated_label(img, (cols[0][0] - 60, (y0 + y1) / 2), tis, row_font, 0, 120, y1 - y0)
    img.save(os.path.join(out_dir, f"{name}.png"))
    print(name, "labels redrawn", img.size)
