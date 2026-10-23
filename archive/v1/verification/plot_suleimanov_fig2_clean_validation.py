#!/usr/bin/env python3
"""Cleaner branch-level validation view for Suleimanov et al. 2009 Fig. 2.

This script avoids the visually noisy point-cloud overlay by:
- filtering vertical marker segments and resonance-marker energy bands,
- splitting each color/component into upper/lower digitized branches per
  x-column without using the model curves,
- plotting branch overlays and residuals against both computed mode curves.

The output is diagnostic, not a final publication-quality digitization.  It is
intended to expose which apparent differences survive after removing obvious
annotation clutter.
"""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "verification" / "data" / "suleimanov_2009_fig2"
FIG_DIR = ROOT / "verification" / "figures" / "suleimanov_2009_fig2"

DIGITIZED = DATA_DIR / "digitized_points.csv"
COMPUTED = DATA_DIR / "computed_opacities.csv"
OUT_FIG = FIG_DIR / "clean_branch_validation.png"
OUT_METRICS = DATA_DIR / "clean_branch_metrics.json"

PANELS = (5.0, 45.0, 90.0)
COLORS = {"magenta": "#d000d0", "blue": "#1647d8"}
COMPONENT_BY_COLOR = {"magenta": "ff_cm2_g", "blue": "scattering_cm2_g"}
COMPONENT_LABEL = {"magenta": "free_free", "blue": "electron_scattering"}
EXCLUDED_ENERGY_BANDS_KEV = ((0.50, 0.80), (5.3, 6.4))
MARKER_SEGMENT_LEN_PX = 8


def read_digitized() -> list[dict]:
    with DIGITIZED.open() as f:
        out = []
        for row in csv.DictReader(f):
            out.append(
                {
                    "theta_deg": float(row["theta_deg"]),
                    "color": row["color"],
                    "x_pixel": int(float(row["x_pixel"])),
                    "energy_keV": float(row["energy_keV"]),
                    "opacity_cm2_g": float(row["opacity_cm2_g"]),
                    "segment_len_px": int(float(row.get("segment_len_px", 1))),
                }
            )
        return out


def read_computed() -> list[dict]:
    with COMPUTED.open() as f:
        out = []
        for row in csv.DictReader(f):
            out.append(
                {
                    "theta_deg": float(row["theta_deg"]),
                    "mode": int(row["mode"]),
                    "energy_keV": float(row["energy_keV"]),
                    "ff_cm2_g": float(row["ff_cm2_g"]),
                    "scattering_cm2_g": float(row["scattering_cm2_g"]),
                }
            )
        return out


def in_excluded_band(energy: float) -> bool:
    return any(lo <= energy <= hi for lo, hi in EXCLUDED_ENERGY_BANDS_KEV)


def curves_by_theta_mode(computed: list[dict]) -> dict[tuple[float, int], list[dict]]:
    grouped: dict[tuple[float, int], list[dict]] = {}
    for row in computed:
        grouped.setdefault((row["theta_deg"], row["mode"]), []).append(row)
    for curve in grouped.values():
        curve.sort(key=lambda r: r["energy_keV"])
    return grouped


def interp_log_curve(curve: list[dict], key: str, energy: float) -> float | None:
    if energy < curve[0]["energy_keV"] or energy > curve[-1]["energy_keV"]:
        return None
    lo, hi = 0, len(curve) - 1
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


def digitized_branches(digitized: list[dict], theta: float, color: str) -> dict[str, list[dict]]:
    by_x: dict[int, list[dict]] = {}
    for row in digitized:
        if row["theta_deg"] != theta or row["color"] != color:
            continue
        if row["segment_len_px"] >= MARKER_SEGMENT_LEN_PX:
            continue
        if in_excluded_band(row["energy_keV"]):
            continue
        by_x.setdefault(row["x_pixel"], []).append(row)

    upper: list[dict] = []
    lower: list[dict] = []
    for x in sorted(by_x):
        # The PDF line strokes produce at most two real branches per color/x
        # after marker filtering.  Use median values if antialiasing leaves
        # duplicate small segments.
        rows = sorted(by_x[x], key=lambda r: r["opacity_cm2_g"], reverse=True)
        if rows:
            upper.append(_median_point(rows[:1], "upper"))
        if len(rows) >= 2:
            lower.append(_median_point(rows[1:], "lower"))
    return {"upper": upper, "lower": lower}


def _median_point(rows: list[dict], branch: str) -> dict:
    return {
        "branch": branch,
        "x_pixel": rows[0]["x_pixel"],
        "energy_keV": statistics.median(r["energy_keV"] for r in rows),
        "opacity_cm2_g": 10.0 ** statistics.median(
            math.log10(r["opacity_cm2_g"]) for r in rows
        ),
    }


def branch_residuals(branch: list[dict], curve: list[dict], component: str) -> list[float]:
    residuals = []
    for point in branch:
        log_model = interp_log_curve(curve, component, point["energy_keV"])
        if log_model is None:
            continue
        residuals.append(math.log10(point["opacity_cm2_g"]) - log_model)
    return residuals


def summarize(residuals: list[float]) -> dict:
    if not residuals:
        return {
            "n": 0,
            "median_dex": None,
            "median_abs_dex": None,
            "robust_scatter_dex": None,
            "p90_abs_dex": None,
        }
    med = statistics.median(residuals)
    abs_res = sorted(abs(r) for r in residuals)
    mad = statistics.median(abs(r - med) for r in residuals)
    return {
        "n": len(residuals),
        "median_dex": med,
        "median_abs_dex": statistics.median(abs_res),
        "robust_scatter_dex": 1.4826 * mad,
        "p90_abs_dex": abs_res[min(len(abs_res) - 1, int(0.9 * (len(abs_res) - 1)))],
    }


def make_metrics(digitized: list[dict], computed: list[dict]) -> dict:
    curves = curves_by_theta_mode(computed)
    rows = []
    for theta in PANELS:
        for color in ("magenta", "blue"):
            component = COMPONENT_BY_COLOR[color]
            branches = digitized_branches(digitized, theta, color)
            for branch_name, points in branches.items():
                for mode in (1, 2):
                    residuals = branch_residuals(points, curves[(theta, mode)], component)
                    row = {
                        "theta_deg": theta,
                        "color": color,
                        "component": COMPONENT_LABEL[color],
                        "digitized_branch": branch_name,
                        "compared_mode": mode,
                        "digitized_points": len(points),
                    }
                    row.update(summarize(residuals))
                    rows.append(row)
    return {
        "method": (
            "Branches are split by opacity ordering per x-column after masking tall "
            "segments and resonance-marker bands.  Residuals compare each branch to "
            "each computed mode separately; no nearest-mode reassignment is used."
        ),
        "excluded_energy_bands_keV": EXCLUDED_ENERGY_BANDS_KEV,
        "marker_segment_len_px_cut": MARKER_SEGMENT_LEN_PX,
        "metrics": rows,
    }


def best_mode(metrics: dict, theta: float, color: str, branch: str) -> int:
    rows = [
        row for row in metrics["metrics"]
        if row["theta_deg"] == theta and row["color"] == color and row["digitized_branch"] == branch
    ]
    rows = [r for r in rows if r["median_abs_dex"] is not None]
    if not rows:
        return 1
    return min(rows, key=lambda r: r["median_abs_dex"])["compared_mode"]


def plot_clean(digitized: list[dict], computed: list[dict], metrics: dict) -> None:
    curves = curves_by_theta_mode(computed)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 2, figsize=(10.5, 9.0), sharex="col")
    branch_markers = {"upper": "o", "lower": "s"}

    for row_idx, theta in enumerate(PANELS):
        ax = axes[row_idx, 0]
        rax = axes[row_idx, 1]

        for color in ("magenta", "blue"):
            component = COMPONENT_BY_COLOR[color]
            branches = digitized_branches(digitized, theta, color)

            for branch_name, points in branches.items():
                if not points:
                    continue
                mode = best_mode(metrics, theta, color, branch_name)
                curve = curves[(theta, mode)]
                xs = [p["energy_keV"] for p in points]
                ys = [p["opacity_cm2_g"] for p in points]
                label = f"{COMPONENT_LABEL[color]} {branch_name} -> mode {mode}"

                ax.scatter(
                    xs,
                    ys,
                    s=6,
                    c=COLORS[color],
                    marker=branch_markers[branch_name],
                    alpha=0.68,
                    linewidths=0,
                    label=label,
                )
                ax.plot(
                    [p["energy_keV"] for p in curve],
                    [p[component] for p in curve],
                    c=COLORS[color],
                    lw=1.3,
                    ls="-" if mode == 1 else "--",
                    alpha=0.9,
                )

                res_x = []
                res_y = []
                for point in points:
                    log_model = interp_log_curve(curve, component, point["energy_keV"])
                    if log_model is None:
                        continue
                    res_x.append(point["energy_keV"])
                    res_y.append(math.log10(point["opacity_cm2_g"]) - log_model)
                rax.scatter(
                    res_x,
                    res_y,
                    s=6,
                    c=COLORS[color],
                    marker=branch_markers[branch_name],
                    alpha=0.68,
                    linewidths=0,
                )

        for panel_ax in (ax, rax):
            panel_ax.set_xscale("log")
            panel_ax.grid(True, which="both", alpha=0.18, linewidth=0.4)
            for lo, hi in EXCLUDED_ENERGY_BANDS_KEV:
                panel_ax.axvspan(lo, hi, color="0.85", alpha=0.25, lw=0)
        ax.set_yscale("log")
        ax.set_ylim(1e-11, 1e9)
        ax.set_ylabel(f"theta={int(theta)} deg\nopacity")
        rax.axhline(0.0, c="0.25", lw=0.8)
        rax.axhline(0.5, c="0.55", lw=0.6, ls=":")
        rax.axhline(-0.5, c="0.55", lw=0.6, ls=":")
        rax.set_ylim(-3.2, 3.2)
        rax.set_ylabel("digitized - model [dex]")

    axes[0, 0].legend(fontsize=7, loc="lower left", frameon=False, ncols=1)
    axes[-1, 0].set_xlabel("Photon energy, keV")
    axes[-1, 1].set_xlabel("Photon energy, keV")
    axes[0, 0].set_title("Branch-matched overlay, resonance bands masked")
    axes[0, 1].set_title("Residuals against best computed mode")
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=190)
    plt.close(fig)


def main() -> None:
    digitized = read_digitized()
    computed = read_computed()
    metrics = make_metrics(digitized, computed)
    OUT_METRICS.write_text(json.dumps(metrics, indent=2) + "\n")
    plot_clean(digitized, computed, metrics)
    print(f"wrote {OUT_FIG.relative_to(ROOT)}")
    print(f"wrote {OUT_METRICS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
