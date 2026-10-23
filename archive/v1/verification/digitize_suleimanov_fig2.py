#!/usr/bin/env python3
"""Digitize Suleimanov, Potekhin & Werner 2009 Fig. 2 from the local PDF.

This is a smoke test for the figure-digitization pipeline.  It extracts
magenta and blue curve pixels from the PDF-rendered figure, maps them onto the
published log-log axes, stores the calibrated point cloud as CSV, and creates a
fresh plot from those points for visual comparison.
"""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "refs" / "suleimanov_potekhin_werner_2009_mag_atm.pdf"
FIG_DIR = ROOT / "verification" / "figures" / "suleimanov_2009_fig2"
DATA_DIR = ROOT / "verification" / "data" / "suleimanov_2009_fig2"

SOURCE_CROP = FIG_DIR / "source_crop.png"
REPLOT = FIG_DIR / "digitized_replot.png"
POINTS_CSV = DATA_DIR / "digitized_points.csv"
METADATA_JSON = DATA_DIR / "digitization_metadata.json"

# 300 dpi render of PDF page 3.  Coordinates are in pixels in that render.
PAGE = 3
DPI = 300
CROP_X = 1450
CROP_Y = 160
CROP_W = 820
CROP_H = 1250

# Calibration in source_crop.png pixel coordinates.  These are the plot borders
# of the three stacked panels, calibrated from the major axis ticks in Fig. 2.
X_LEFT = 48
X_RIGHT = 792
LOG10_E_MIN = -2.0
LOG10_E_MAX = math.log10(45.0)

PANELS = [
    {"theta_deg": 5, "y_top": 31, "y_bottom": 413},
    {"theta_deg": 45, "y_top": 413, "y_bottom": 793},
    {"theta_deg": 90, "y_top": 793, "y_bottom": 1176},
]
LOG10_KAPPA_MIN = -11.0
LOG10_KAPPA_MAX = 9.0


def render_and_crop() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    tmp_prefix = Path("/tmp") / "suleimanov_fig2_page"
    subprocess.run(
        [
            "pdftocairo",
            "-png",
            "-r",
            str(DPI),
            "-f",
            str(PAGE),
            "-l",
            str(PAGE),
            str(PDF),
            str(tmp_prefix),
        ],
        check=True,
    )
    page_png = Path(f"{tmp_prefix}-{PAGE}.png")
    with Image.open(page_png).convert("RGB") as im:
        crop = im.crop((CROP_X, CROP_Y, CROP_X + CROP_W, CROP_Y + CROP_H))
        crop.save(SOURCE_CROP)


def pixel_to_data(x: int, y: float, panel: dict) -> tuple[float, float]:
    log_e = LOG10_E_MIN + (x - X_LEFT) / (X_RIGHT - X_LEFT) * (
        LOG10_E_MAX - LOG10_E_MIN
    )
    log_kappa = LOG10_KAPPA_MIN + (panel["y_bottom"] - y) / (
        panel["y_bottom"] - panel["y_top"]
    ) * (LOG10_KAPPA_MAX - LOG10_KAPPA_MIN)
    return 10.0**log_e, 10.0**log_kappa


def classify_pixel(r: int, g: int, b: int) -> str | None:
    # The PDF curves are saturated magenta and blue.  Thresholds include
    # anti-aliased edge pixels but reject black axes/text and white background.
    if r > 170 and b > 150 and g < 120:
        return "magenta"
    if b > 150 and r < 90 and g < 140:
        return "blue"
    return None


def contiguous_segments(values: list[int]) -> list[list[int]]:
    if not values:
        return []
    segments = [[values[0]]]
    for value in values[1:]:
        if value <= segments[-1][-1] + 1:
            segments[-1].append(value)
        else:
            segments.append([value])
    return segments


def extract_points() -> list[dict]:
    with Image.open(SOURCE_CROP).convert("RGB") as im:
        pix = im.load()
        points: list[dict] = []

        for panel_index, panel in enumerate(PANELS, start=1):
            y_top = panel["y_top"]
            y_bottom = panel["y_bottom"]

            for color in ("magenta", "blue"):
                for x in range(X_LEFT, X_RIGHT + 1):
                    ys = []
                    for y in range(y_top, y_bottom + 1):
                        if classify_pixel(*pix[x, y]) == color:
                            ys.append(y)

                    for segment in contiguous_segments(ys):
                        # Single-pixel noise is common around antialiasing.
                        if len(segment) < 2:
                            continue
                        y_mid = sum(segment) / len(segment)
                        energy, opacity = pixel_to_data(x, y_mid, panel)
                        points.append(
                            {
                                "figure": "Suleimanov2009_Fig2",
                                "panel": panel_index,
                                "theta_deg": panel["theta_deg"],
                                "color": color,
                                "x_pixel": x,
                                "y_min_pixel": min(segment),
                                "y_max_pixel": max(segment),
                                "y_pixel": y_mid,
                                "segment_len_px": len(segment),
                                "energy_keV": energy,
                                "opacity_cm2_g": opacity,
                            }
                        )

        return points


def write_points(points: list[dict]) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with POINTS_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(points[0].keys()), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(points)


def plot_points(points: list[dict]) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(5.4, 8.0), sharex=True)
    color_map = {"magenta": "#ff00ff", "blue": "#0040ff"}

    for ax, panel in zip(axes, PANELS):
        for color in ("magenta", "blue"):
            subset = [
                p
                for p in points
                if p["theta_deg"] == panel["theta_deg"] and p["color"] == color
            ]
            ax.scatter(
                [p["energy_keV"] for p in subset],
                [p["opacity_cm2_g"] for p in subset],
                s=1.0,
                c=color_map[color],
                linewidths=0,
                alpha=0.9,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-2, 45.0)
        ax.set_ylim(1e-11, 1e9)
        ax.text(0.015, 2e6, f"theta = {panel['theta_deg']} deg", fontsize=12)
        ax.grid(True, which="both", alpha=0.18, linewidth=0.4)

    axes[1].set_ylabel("Opacity, cm^2 g^-1")
    axes[-1].set_xlabel("Photon energy, keV")
    fig.suptitle("Digitized point cloud from Suleimanov et al. 2009 Fig. 2")
    fig.tight_layout()
    fig.savefig(REPLOT, dpi=180)
    plt.close(fig)


def write_metadata(points: list[dict]) -> None:
    metadata = {
        "source_pdf": str(PDF.relative_to(ROOT)),
        "figure": "Suleimanov et al. 2009 Fig. 2",
        "page": PAGE,
        "dpi": DPI,
        "source_crop": str(SOURCE_CROP.relative_to(ROOT)),
        "digitized_csv": str(POINTS_CSV.relative_to(ROOT)),
        "replot": str(REPLOT.relative_to(ROOT)),
        "crop_pixels": {
            "x": CROP_X,
            "y": CROP_Y,
            "width": CROP_W,
            "height": CROP_H,
        },
        "axis_calibration": {
            "x_left_px": X_LEFT,
            "x_right_px": X_RIGHT,
            "energy_keV_min": 10.0**LOG10_E_MIN,
            "energy_keV_max": 10.0**LOG10_E_MAX,
            "opacity_cm2_g_min": 10.0**LOG10_KAPPA_MIN,
            "opacity_cm2_g_max": 10.0**LOG10_KAPPA_MAX,
            "panels": PANELS,
        },
        "point_count": len(points),
        "method": "color thresholding of saturated magenta/blue curve pixels, then log-axis calibration",
        "limitations": [
            "Curves are not yet semantically separated into FF/ES or mode labels.",
            "Vertical resonance markers are retained but can be filtered using segment_len_px.",
            "Axis calibration is manual from figure border and major tick marks.",
        ],
    }
    METADATA_JSON.write_text(json.dumps(metadata, indent=2) + "\n")


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    render_and_crop()
    points = extract_points()
    if not points:
        raise RuntimeError("No colored curve points extracted")
    write_points(points)
    plot_points(points)
    write_metadata(points)
    print(f"wrote {SOURCE_CROP.relative_to(ROOT)}")
    print(f"wrote {POINTS_CSV.relative_to(ROOT)} ({len(points)} points)")
    print(f"wrote {REPLOT.relative_to(ROOT)}")
    print(f"wrote {METADATA_JSON.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
