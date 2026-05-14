#!/usr/bin/env python3
"""Compare digitized Suleimanov et al. 2009 Fig. 2 points to local opacities.

The digitized data are still a color point cloud, not a fully semantic curve
trace.  This script therefore produces a visual overlay and a preliminary
nearest-curve error metric after excluding the published resonance marker
bands.
"""

from __future__ import annotations

import csv
import statistics
import json
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "verification" / "data" / "suleimanov_2009_fig2"
FIG_DIR = ROOT / "verification" / "figures" / "suleimanov_2009_fig2"

DIGITIZED = DATA_DIR / "digitized_points.csv"
COMPUTED = DATA_DIR / "computed_opacities.csv"
METRICS = DATA_DIR / "opacity_comparison_metrics.json"
OVERLAY = FIG_DIR / "opacity_comparison.png"

PANELS = (5.0, 45.0, 90.0)
COLORS = {"magenta": "#d000d0", "blue": "#1647d8"}
COMPONENT_BY_COLOR = {"magenta": "ff_cm2_g", "blue": "scattering_cm2_g"}
COMPONENT_LABEL = {"magenta": "free_free", "blue": "electron_scattering"}
MODE_LINESTYLE = {1: "-", 2: "--"}

# Fig. 2 marks the proton cyclotron and vacuum resonance as vertical colored
# features.  They are useful annotations, but not opacity curves to score.
EXCLUDED_ENERGY_BANDS_KEV = (
    (0.50, 0.80),
    (5.3, 6.4),
)
MARKER_SEGMENT_LEN_PX = 8


def read_digitized() -> list[dict]:
    with DIGITIZED.open() as f:
        rows = []
        for row in csv.DictReader(f):
            rows.append(
                {
                    "theta_deg": float(row["theta_deg"]),
                    "color": row["color"],
                    "x_pixel": int(float(row.get("x_pixel", 0))),
                    "energy_keV": float(row["energy_keV"]),
                    "opacity_cm2_g": float(row["opacity_cm2_g"]),
                    "segment_len_px": int(float(row.get("segment_len_px", 1))),
                }
            )
    return rows


def read_computed() -> list[dict]:
    with COMPUTED.open() as f:
        rows = []
        for row in csv.DictReader(f):
            rows.append(
                {
                    "theta_deg": float(row["theta_deg"]),
                    "mode": int(row["mode"]),
                    "energy_keV": float(row["energy_keV"]),
                    "ff_cm2_g": float(row["ff_cm2_g"]),
                    "pp_cm2_g": float(row["pp_cm2_g"]),
                    "absorption_cm2_g": float(row["absorption_cm2_g"]),
                    "scattering_cm2_g": float(row["scattering_cm2_g"]),
                    "total_cm2_g": float(row["total_cm2_g"]),
                }
            )
    return rows


def in_excluded_band(energy: float) -> bool:
    return any(lo <= energy <= hi for lo, hi in EXCLUDED_ENERGY_BANDS_KEV)


def interp_log_curve(curve: list[dict], key: str, energy: float) -> float | None:
    if energy < curve[0]["energy_keV"] or energy > curve[-1]["energy_keV"]:
        return None
    lo = 0
    hi = len(curve) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if curve[mid]["energy_keV"] <= energy:
            lo = mid
        else:
            hi = mid

    x0 = math.log10(curve[lo]["energy_keV"])
    x1 = math.log10(curve[hi]["energy_keV"])
    y0 = math.log10(max(curve[lo][key], 1e-300))
    y1 = math.log10(max(curve[hi][key], 1e-300))
    x = math.log10(energy)
    if x1 == x0:
        return y0
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)


def curves_by_theta_mode(computed: list[dict]) -> dict[tuple[float, int], list[dict]]:
    grouped: dict[tuple[float, int], list[dict]] = {}
    for row in computed:
        grouped.setdefault((row["theta_deg"], row["mode"]), []).append(row)
    for curve in grouped.values():
        curve.sort(key=lambda row: row["energy_keV"])
    return grouped


def plot_overlay(digitized: list[dict], computed: list[dict]) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    grouped = curves_by_theta_mode(computed)

    fig, axes = plt.subplots(3, 1, figsize=(6.2, 8.8), sharex=True)
    for ax, theta in zip(axes, PANELS):
        for color, hex_color in COLORS.items():
            pts = [
                p for p in digitized
                if p["theta_deg"] == theta and p["color"] == color
            ]
            ax.scatter(
                [p["energy_keV"] for p in pts],
                [p["opacity_cm2_g"] for p in pts],
                s=1.0,
                c=hex_color,
                alpha=0.22,
                linewidths=0,
            )

        for mode in (1, 2):
            curve = grouped[(theta, mode)]
            x = [row["energy_keV"] for row in curve]
            ax.plot(
                x,
                [row["ff_cm2_g"] for row in curve],
                color=COLORS["magenta"],
                lw=1.5,
                ls=MODE_LINESTYLE[mode],
            )
            ax.plot(
                x,
                [row["scattering_cm2_g"] for row in curve],
                color=COLORS["blue"],
                lw=1.5,
                ls=MODE_LINESTYLE[mode],
            )
            ax.plot(
                x,
                [row["pp_cm2_g"] for row in curve],
                color="0.35",
                lw=0.9,
                ls=MODE_LINESTYLE[mode],
                alpha=0.75,
            )

        for lo, hi in EXCLUDED_ENERGY_BANDS_KEV:
            ax.axvspan(lo, hi, color="0.85", alpha=0.22, lw=0)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-2, 45.0)
        ax.set_ylim(1e-11, 1e9)
        ax.grid(True, which="both", alpha=0.18, linewidth=0.4)
        ax.text(0.014, 1.7e7, f"theta = {int(theta)} deg", fontsize=11)

    axes[0].plot([], [], c=COLORS["magenta"], lw=1.5, label="computed FF")
    axes[0].plot([], [], c=COLORS["blue"], lw=1.5, label="computed ES")
    axes[0].plot([], [], c="0.35", lw=0.9, label="computed pp")
    axes[0].plot([], [], c="black", lw=1.5, ls="-", label="mode 1")
    axes[0].plot([], [], c="black", lw=1.5, ls="--", label="mode 2")
    axes[0].scatter([], [], c="black", s=5, alpha=0.22, label="digitized")
    axes[0].legend(loc="lower left", fontsize=8, frameon=False)

    axes[1].set_ylabel("Opacity, cm^2 g^-1")
    axes[-1].set_xlabel("Photon energy, keV")
    fig.suptitle("Suleimanov et al. 2009 Fig. 2: digitized data vs local opacity code")
    fig.tight_layout()
    fig.savefig(OVERLAY, dpi=190)
    plt.close(fig)


def preliminary_metrics(digitized: list[dict], computed: list[dict]) -> dict:
    grouped = curves_by_theta_mode(computed)
    metric_rows = []

    for theta in PANELS:
        for color in ("magenta", "blue"):
            component = COMPONENT_BY_COLOR[color]
            residual_bins_by_mode: dict[int, dict[int, list[float]]] = {1: {}, 2: {}}
            usable_x = set()

            for point in digitized:
                if point["theta_deg"] != theta or point["color"] != color:
                    continue
                if in_excluded_band(point["energy_keV"]):
                    continue
                if point["segment_len_px"] >= MARKER_SEGMENT_LEN_PX:
                    continue
                usable_x.add(point["x_pixel"])
                log_point = math.log10(point["opacity_cm2_g"])
                candidates = []
                for mode in (1, 2):
                    curve = grouped[(theta, mode)]
                    log_curve = interp_log_curve(curve, component, point["energy_keV"])
                    if log_curve is not None:
                        candidates.append((abs(log_point - log_curve), mode, log_curve))
                if not candidates:
                    continue
                _, mode, log_curve = min(candidates, key=lambda item: item[0])
                residual_bins_by_mode[mode].setdefault(point["x_pixel"], []).append(
                    log_point - log_curve
                )

            for mode, residual_bins in residual_bins_by_mode.items():
                # One representative residual per x-column avoids overweighting
                # antialiasing thickness and multi-pixel colored strokes.
                errors = [statistics.median(vals) for vals in residual_bins.values()]
                if errors:
                    median = statistics.median(errors)
                    abs_errors = sorted(abs(err) for err in errors)
                    median_abs = statistics.median(abs_errors)
                    mad = statistics.median(abs(err - median) for err in errors)
                    robust_scatter = 1.4826 * mad
                    p90_abs = abs_errors[min(len(abs_errors) - 1, int(0.9 * (len(abs_errors) - 1)))]
                    coverage = len(errors) / max(1, len(usable_x))
                else:
                    median = None
                    median_abs = None
                    robust_scatter = None
                    p90_abs = None
                    coverage = 0.0
                metric_rows.append(
                    {
                        "theta_deg": theta,
                        "digitized_color": color,
                        "mapped_component": COMPONENT_LABEL[color],
                        "assigned_mode": mode,
                        "assigned_energy_columns": len(errors),
                        "usable_energy_columns_for_color": len(usable_x),
                        "coverage_fraction": coverage,
                        "median_log10_opacity_residual_dex": median,
                        "median_abs_log10_opacity_error_dex": median_abs,
                        "robust_scatter_1p4826_mad_dex": robust_scatter,
                        "p90_abs_log10_opacity_error_dex": p90_abs,
                    }
                )

    return {
        "source": "Suleimanov et al. 2009 Fig. 2 digitized point cloud",
        "computed_parameters": {
            "temperature_K": 7.0e6,
            "density_g_cm3": 30.0,
            "magnetic_field_G": 1.0e14,
            "angles_deg": list(PANELS),
        },
        "figure_mapping": {
            "magenta": "free-free/true absorption curves",
            "blue": "electron-scattering curves",
            "solid": "one normal-mode branch",
            "dashed": "the other normal-mode branch",
            "local_line_styles": "mode 1 solid, mode 2 dashed for comparison only",
        },
        "excluded_energy_bands_keV": EXCLUDED_ENERGY_BANDS_KEV,
        "method": (
            "For each digitized point outside resonance-marker bands, use its "
            "published color to choose free-free or electron scattering, then "
            "assign it to whichever local mode curve is nearer in log opacity. "
            "Tall pixel segments are also excluded as vertical marker candidates. "
            "Metrics use one median residual per x-column to reduce linewidth bias."
        ),
        "marker_segment_len_px_cut": MARKER_SEGMENT_LEN_PX,
        "limitations": [
            "Digitized curves have not yet been manually separated by line style.",
            "The local opacity code includes vacuum-polarized mode vectors, but this metric excludes the resonance-marker bands.",
            "The local comparison uses Potekhin-Chabrier component cross-sections, not van Adelsberg-Lai tables.",
            "The proton-proton opacity is overplotted but not scored against Fig. 2 labels.",
        ],
        "metrics": metric_rows,
    }


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    digitized = read_digitized()
    computed = read_computed()
    plot_overlay(digitized, computed)
    metrics = preliminary_metrics(digitized, computed)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"wrote {OVERLAY.relative_to(ROOT)}")
    print(f"wrote {METRICS.relative_to(ROOT)}")
    print("figure mapping:", metrics["figure_mapping"])


if __name__ == "__main__":
    main()
