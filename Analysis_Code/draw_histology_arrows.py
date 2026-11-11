"""Overlay lesion arrows on the scale-barred histology composites and build numbered review sheets.

Arrow tips are given as fractions of each Supplementary Figure 6 tile. The same micrographs appear at half size in
Figure 3E/F; each S6 tile is matched to its Figure 3 tile by a small shift search so the tips land on the same feature.

usage (run from the repository root, after redraw_histology_scalebars.py has written the *_scalebars.png images):
  python3 draw_histology_arrows.py SCALEBAR_SPEC.json ARROWS.json OUT_DIR [--review]
Writes OUT_DIR/histology_panels_EF_arrows.png, SupplFig6A_arrows.png, SupplFig6B_arrows.png and, with --review,
numbered review sheets (Healthy | DD | GG | GD per tissue) for the author.
"""
import json
import math
import os
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None
spec_path, arrows_path, out_dir = sys.argv[1:4]
review = "--review" in sys.argv
spec = json.load(open(spec_path, encoding="utf-8"))
_aj = json.load(open(arrows_path, encoding="utf-8"))
arrows = _aj["arrows"]
# arrow length as a fraction of the tile width, and white outline width (defaults = first-revision sizes)
STYLE = {"fig3_len_frac": None, "s6_len_frac": .115, "fig3_outline": 2, "s6_outline": 3}
STYLE.update(_aj.get("style", {}))
fig3_job, s6a_job, s6b_job = spec["jobs"]
os.makedirs(out_dir, exist_ok=True)
FONT_B = "/usr/share/fonts/truetype/noto/NotoSans-Bold.ttf"

DIRS = [(-1, -1), (1, -1), (-1, 1), (-1, 0), (0, -1), (1, 1)]


def arrow_poly(tip, tail, L):
    tx, ty = tip; ux, uy = (tip[0] - tail[0]) / L, (tip[1] - tail[1]) / L
    nx, ny = -uy, ux
    hl, hw, sw = .40 * L, .20 * L, .065 * L
    neck = (tx - ux * hl, ty - uy * hl)
    return [(tail[0] + nx * sw, tail[1] + ny * sw), (neck[0] + nx * sw, neck[1] + ny * sw), (neck[0] + nx * hw, neck[1] + ny * hw),
            (tx, ty), (neck[0] - nx * hw, neck[1] - ny * hw), (neck[0] - nx * sw, neck[1] - ny * sw), (tail[0] - nx * sw, tail[1] - ny * sw)]


def free(tile, tail):
    x0, y0, x1, y1 = tile; w, h = x1 - x0 + 1, y1 - y0 + 1
    fx, fy = (tail[0] - x0) / w, (tail[1] - y0) / h
    return .03 < fx < .97 and .03 < fy < .97 and not (fx > .55 and fy > .70)


def place(tile, frac, L, forced=None, also=None):
    """Tip at the fractional position; first direction whose tail stays inside the tile and away from the
    bottom-right scale bar. `also(direction)` can veto a direction (used to check the Figure 3 copy)."""
    x0, y0, x1, y1 = tile; w, h = x1 - x0 + 1, y1 - y0 + 1
    tip = (x0 + frac[0] * w, y0 + frac[1] * h)
    for dx, dy in ([forced] if forced else DIRS):
        n = math.hypot(dx, dy); tail = (tip[0] + dx / n * L, tip[1] + dy / n * L)
        if free(tile, tail) and (also is None or also((dx, dy))):
            return tip, tail, (dx, dy)
    raise ValueError(f"no free direction for arrow at {frac}")


def draw_arrow(d, tip, tail, L, outline):
    d.polygon(arrow_poly(tip, tail, L), fill=(0, 0, 0), outline=(255, 255, 255), width=outline)


def align(fig3, t3, s6, t6, search=10):
    """Shift (dx, dy) in Figure 3 pixels that best matches the S6 tile downsampled to the Figure 3 tile size."""
    x0, y0, x1, y1 = t3; w, h = x1 - x0 + 1, y1 - y0 + 1
    a = np.asarray(fig3.convert("L").crop((x0, y0, x1 + 1, y1 + 1))).astype(float)
    b = np.asarray(s6.convert("L").crop((t6[0], t6[1], t6[2] + 1, t6[3] + 1)).resize((w, h), Image.BILINEAR)).astype(float)
    m = search; best = (-2, 0, 0)
    core = b[m:h - m - int(h * .25), m:w - m]
    for dy in range(-m, m + 1):
        for dx in range(-m, m + 1):
            sub = a[m + dy:h - m - int(h * .25) + dy, m + dx:w - m + dx]
            c = np.corrcoef(sub.ravel(), core.ravel())[0, 1]
            if c > best[0]:
                best = (c, dx, dy)
    return best


fig3 = Image.open(fig3_job["dst"]).convert("RGB")
s6 = {"1-year": (Image.open(s6a_job["dst"]).convert("RGB"), s6a_job, 0), "2-year": (Image.open(s6b_job["dst"]).convert("RGB"), s6b_job, 16)}
f3_out = fig3.copy(); d3 = ImageDraw.Draw(f3_out)
s6_out = {k: v[0].copy() for k, v in s6.items()}
placed = []
for ar in arrows:
    img6, job6, base = s6[ar["age"]]
    t6 = job6["tiles"][ar["s6_tile"]]
    r, c = divmod(ar["s6_tile"], 3)
    i3 = base + r * 4 + c + 1
    t3 = fig3_job["tiles"][i3]
    corr, sx, sy = align(fig3, t3, img6, t6)
    w3, h3 = t3[2] - t3[0] + 1, t3[3] - t3[1] + 1
    w6, h6 = t6[2] - t6[0] + 1, t6[3] - t6[1] + 1
    k = w3 / w6
    L6 = STYLE["s6_len_frac"] * w6
    L3 = STYLE["fig3_len_frac"] * w3 if STYLE["fig3_len_frac"] else L6 * k
    # Figure 3: same fractional position, corrected by the matched shift
    tip3 = (t3[0] + ar["tip"][0] * w6 * k + sx, t3[1] + ar["tip"][1] * h6 * k + sy)

    def fig3_ok(dv):
        n = math.hypot(*dv)
        return free(t3, (tip3[0] + dv[0] / n * L3, tip3[1] + dv[1] / n * L3))

    tip6, tail6, dirv = place(t6, ar["tip"], L6, tuple(ar["dir"]) if ar.get("dir") else None, fig3_ok)
    draw_arrow(ImageDraw.Draw(s6_out[ar["age"]]), tip6, tail6, L6, STYLE["s6_outline"])
    n = math.hypot(*dirv); tail3 = (tip3[0] + dirv[0] / n * L3, tip3[1] + dirv[1] / n * L3)
    draw_arrow(d3, tip3, tail3, L3, STYLE["fig3_outline"])
    placed.append({"id": ar["id"], "fig3_tile": i3, "align_corr": round(float(corr), 3), "shift_px": [sx, sy],
                   "s6_tip": [round(v) for v in tip6], "s6_tail": [round(v) for v in tail6], "fig3_tip": [round(v) for v in tip3]})
    print(f"arrow {ar['id']:2d} {ar['age']} {ar['tissue']:15s} {ar['genotype']} S6 tile {ar['s6_tile']:2d} -> Fig3 tile {i3:2d}  match r={corr:.3f} shift=({sx},{sy})")

f3_out.save(os.path.join(out_dir, "histology_panels_EF_arrows.png"))
s6_out["1-year"].save(os.path.join(out_dir, "SupplFig6A_arrows.png"))
s6_out["2-year"].save(os.path.join(out_dir, "SupplFig6B_arrows.png"))
json.dump(placed, open(os.path.join(out_dir, "arrow_placement.json"), "w"), indent=1)

if review:
    font = ImageFont.truetype(FONT_B, 40); small = ImageFont.truetype(FONT_B, 34)
    tissues = ["Foot muscle", "Digestive gland", "Gill", "Mantle"]
    for age, (img6, job6, base) in s6.items():
        TW, TH = 831, 623; gap = 14; head = 70; lab = 230
        sheet = Image.new("RGB", (lab + 4 * (TW + gap), head + 4 * (TH + gap)), "white")
        d = ImageDraw.Draw(sheet)
        for ci, g in enumerate(["Healthy", "DD", "GG", "GD"]):
            d.text((lab + ci * (TW + gap) + TW // 2 - 60, 12), g, fill=(0, 0, 0), font=font)
        for r, tis in enumerate(tissues):
            y = head + r * (TH + gap)
            d.text((10, y + TH // 2 - 40), tis.replace(" ", "\n"), fill=(0, 0, 0), font=small)
            t3 = fig3_job["tiles"][base + r * 4]
            sheet.paste(f3_out.crop((t3[0], t3[1], t3[2] + 1, t3[3] + 1)).resize((TW, TH), Image.LANCZOS), (lab, y))
            for c in range(3):
                t6 = job6["tiles"][r * 3 + c]
                tile = s6_out[age].crop((t6[0], t6[1], t6[2] + 1, t6[3] + 1)).resize((TW, TH))
                td = ImageDraw.Draw(tile)
                for p in placed:
                    ar = next(a for a in arrows if a["id"] == p["id"])
                    if ar["age"] == age and ar["s6_tile"] == r * 3 + c:
                        tx, ty = p["s6_tail"][0] - t6[0], p["s6_tail"][1] - t6[1]
                        ux, uy = tx - (p["s6_tip"][0] - t6[0]), ty - (p["s6_tip"][1] - t6[1])
                        nn = math.hypot(ux, uy)
                        lx, ly = tx + ux / nn * 26 - 14, ty + uy / nn * 26 - 22
                        td.rectangle([lx - 4, ly + 2, lx + 34 + 12 * (p["id"] >= 10), ly + 44], fill=(255, 255, 255))
                        td.text((lx, ly), str(p["id"]), fill=(220, 0, 0), font=font)
                sheet.paste(tile, (lab + (c + 1) * (TW + gap), y))
        name = "箭头确认_1龄_Fig3E_S6A.png" if age == "1-year" else "箭头确认_2龄_Fig3F_S6B.png"
        sheet.save(os.path.join(out_dir, name))
        print("review sheet", name, sheet.size)
