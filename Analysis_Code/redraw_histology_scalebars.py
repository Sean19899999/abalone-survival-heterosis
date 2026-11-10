"""Redraw legible, labelled scale bars on the histology composites (Figure 3E-F, Supplementary Figure 6A-B).

The original micrographs carry the microscope's scale bars, which are too thin and small to read at print size.
Each bar is drawn thicker and with a larger label on an opaque white box in the bottom-right corner that also covers
the original bar and label. Image content outside the corner box is unchanged.
In "uniform" mode (used for the revision, value_um = 100) every tile carries the same value, and each bar's length is
converted from that tile's own microscope calibration (native_px / native_um), so bar lengths differ where the image
calibrations differ. In "native" mode each bar keeps the value and pixel length of the microscope's own bar.
(With neither mode, a fixed-length bar is drawn and its value is converted from each tile's calibration.)

usage: python3 redraw_histology_scalebars.py SPEC.json
SPEC: {"jobs": [{"src", "dst", "tiles": [[x0,y0,x1,y1], ...], "native_px": [...], "native_um": [...],
                 "bar_px", "font_px", "bar_thick", "font"}]}
"""
import json
import sys
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None


def label_for(um):
    return f"{um:g} µm"


def run(job):
    im = Image.open(job["src"]).convert("RGB")
    d = ImageDraw.Draw(im)
    font = ImageFont.truetype(job["font"], job["font_px"])
    bar_px, thick = job.get("bar_px"), job["bar_thick"]
    pad, gap = round(job["font_px"] * job.get("pad_frac", .12)), round(job["font_px"] * job.get("gap_frac", .08))
    report = []
    for i, (x0, y0, x1, y1) in enumerate(job["tiles"]):
        if job.get("mode") == "native":
            # same value and length as the microscope's own bar; only legibility changes
            bar_px, um = job["native_px"][i], job["native_um"][i]
        elif job.get("mode") == "uniform":
            # one value on every tile; each bar length comes from that tile's own microscope calibration
            um = job["value_um"]
            bar_px = round(job["native_px"][i] * um / job["native_um"][i])
        else:
            px_per_um = job["native_px"][i] / job["native_um"][i]
            um = bar_px / px_per_um
            assert abs(um - round(um)) < 1e-6, (i, um)
        text = label_for(round(um))
        l, t, r, b = d.textbbox((0, 0), text, font=font)
        tw, th = r - l, b - t
        bw = max(tw, bar_px) + 2 * pad
        bh = th + gap + thick + 2 * pad
        bx0, by0 = x1 - bw + 1, y1 - bh + 1
        if job.get("native_x0") and job["native_x0"][i] is not None:
            bx0 = min(bx0, job["native_x0"][i] - pad)  # the box must also hide the whole original bar
        d.rectangle([bx0, by0, x1, y1], fill=(255, 255, 255))
        cx = (bx0 + x1) / 2
        d.text((cx - tw / 2 - l, by0 + pad - t), text, fill=(0, 0, 0), font=font)
        by = by0 + pad + th + gap
        d.rectangle([round(cx - bar_px / 2), by, round(cx - bar_px / 2) + bar_px - 1, by + thick - 1], fill=(0, 0, 0))
        report.append({"tile": i, "value_um": round(um), "box": [bx0, by0, x1, y1]})
    im.save(job["dst"], dpi=im.info.get("dpi", (600, 600)))
    return report


if __name__ == "__main__":
    spec = json.load(open(sys.argv[1], encoding="utf-8"))
    out = {}
    for job in spec["jobs"]:
        out[job["dst"]] = run(job)
        print(job["dst"], [r["value_um"] for r in out[job["dst"]]])
    if len(sys.argv) > 2:
        json.dump(out, open(sys.argv[2], "w"), indent=1)
