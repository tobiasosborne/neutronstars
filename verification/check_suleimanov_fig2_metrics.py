#!/usr/bin/env python3
"""Fail/pass gate for the Suleimanov et al. 2009 Fig. 2 opacity comparison."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "verification" / "data" / "suleimanov_2009_fig2" / "opacity_comparison_metrics.json"

TARGETS = [
    {
        "theta_deg": 90.0,
        "digitized_color": "magenta",
        "mapped_component": "free_free",
        "assigned_mode": 1,
    },
    {
        "theta_deg": 5.0,
        "digitized_color": "magenta",
        "mapped_component": "free_free",
        "assigned_mode": 1,
    },
]


def find_row(rows: list[dict], target: dict) -> dict:
    for row in rows:
        if all(row.get(key) == value for key, value in target.items()):
            return row
    raise AssertionError(f"missing metric row: {target}")


def main() -> None:
    payload = json.loads(METRICS.read_text())
    failures = []

    for target in TARGETS:
        row = find_row(payload["metrics"], target)
        checks = {
            "coverage_fraction >= 0.45": row["coverage_fraction"] >= 0.45,
            "abs(median residual) <= 0.70 dex": abs(row["median_log10_opacity_residual_dex"]) <= 0.70,
            "median abs error <= 0.75 dex": row["median_abs_log10_opacity_error_dex"] <= 0.75,
            "robust scatter <= 0.60 dex": row["robust_scatter_1p4826_mad_dex"] <= 0.60,
        }
        if target["theta_deg"] == 90.0 and target["assigned_mode"] == 1:
            checks["p90 abs error <= 1.25 dex"] = row["p90_abs_log10_opacity_error_dex"] <= 1.25

        for label, ok in checks.items():
            if not ok:
                failures.append((target, label, row))

    if failures:
        for target, label, row in failures:
            print(f"FAIL {target}: {label}")
            print(
                "  median={:.3f} dex median_abs={:.3f} dex scatter={:.3f} dex p90={:.3f} dex coverage={:.3f}".format(
                    row["median_log10_opacity_residual_dex"],
                    row["median_abs_log10_opacity_error_dex"],
                    row["robust_scatter_1p4826_mad_dex"],
                    row["p90_abs_log10_opacity_error_dex"],
                    row["coverage_fraction"],
                )
            )
        raise SystemExit(1)

    print("Suleimanov Fig. 2 free-free metric gate passed")


if __name__ == "__main__":
    main()
