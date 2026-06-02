#!/usr/bin/env python3
"""
aggregate_pl_sensitivity.py
===========================

Walk through the pl_mult x migration sweep output directory, read every
*_summary.tsv file, and build a single tidy dataframe ready for a 3x3
buffering heatmap.

Headline diagnostic
-------------------
Prints a 3x3 grid of "rescue effect" per lineage:
  rescue = extinction_rate_without - extinction_rate_with
A POSITIVE value means with-introgression reduced extinction. A
negative or zero value means introgression did not help (or worsened
things, which would be strange and worth investigating).

Use --partial mid-sweep to read what's complete so far.

Usage:
    python 04_futureClimate/aggregate_pl_sensitivity.py
    python 04_futureClimate/aggregate_pl_sensitivity.py --partial
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_LINEAGES = ["pdav-N", "mixed", "rot"]

# Cell directory: pl-0.5__mig-1e-06
CELL_RE = re.compile(
    r"^pl-(?P<pl_mult>[^_]+)__mig-(?P<mig>[^_]+)$"
)
REP_RE = re.compile(r"^(?P<script>with|without)__rep-(?P<rep>\d+)$")


def parse_cell_name(name: str) -> dict | None:
    m = CELL_RE.match(name)
    if not m:
        return None
    return {"pl_mult": float(m.group("pl_mult")),
            "mig":     float(m.group("mig"))}


def parse_rep_name(name: str) -> dict | None:
    m = REP_RE.match(name)
    if not m:
        return None
    return {"script": m.group("script"), "rep": int(m.group("rep"))}


def read_one_rep(rep_dir: Path, cell_params: dict, rep_params: dict) -> list[dict]:
    if not (rep_dir / "DONE").exists():
        return []

    candidates = list(rep_dir.glob("*_summary.tsv"))
    if not candidates:
        return [{"warning": f"DONE but no _summary.tsv in {rep_dir}"}]
    if len(candidates) > 1:
        return [{"warning": f"multiple _summary.tsv in {rep_dir}"}]

    summary = pd.read_csv(candidates[0], sep="\t")
    found = set(summary["lineage"].tolist())
    rows = []
    for lin in EXPECTED_LINEAGES:
        base = {**cell_params, **rep_params, "lineage": lin}
        if lin in found:
            r = summary[summary["lineage"] == lin].iloc[0]
            rows.append({
                **base,
                "fitness_end":   float(r["fitness_end"]),
                "phenotype_end": float(r["phenotype_end"]),
                "optimum_end":   float(r["optimum_end"]),
                "fitness_lag":   float(r["fitness_lag"]),
                "extinct":       False,
            })
        else:
            rows.append({
                **base,
                "fitness_end":   np.nan,
                "phenotype_end": np.nan,
                "optimum_end":   np.nan,
                "fitness_lag":   np.nan,
                "extinct":       True,
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-dir", type=Path,
                    default=Path("output/pl_sensitivity"))
    ap.add_argument("--out", type=Path,
                    default=Path("output/pl_sensitivity/pl_sweep_results.csv"))
    ap.add_argument("--partial", action="store_true")
    args = ap.parse_args()

    if not args.sweep_dir.exists():
        sys.exit(f"ERROR: sweep directory not found: {args.sweep_dir}")

    all_rows: list[dict] = []
    cells_complete = 0
    cells_incomplete = 0
    reps_complete = 0
    reps_incomplete = 0
    warnings: list[str] = []

    for cell_dir in sorted(args.sweep_dir.iterdir()):
        if not cell_dir.is_dir():
            continue
        cell_params = parse_cell_name(cell_dir.name)
        if cell_params is None:
            continue

        rep_complete = 0
        for rep_dir in sorted(cell_dir.iterdir()):
            if not rep_dir.is_dir():
                continue
            rep_params = parse_rep_name(rep_dir.name)
            if rep_params is None:
                continue
            rows = read_one_rep(rep_dir, cell_params, rep_params)
            if not rows:
                reps_incomplete += 1
                continue
            if any("warning" in r for r in rows):
                warnings.extend(r["warning"] for r in rows if "warning" in r)
            else:
                all_rows.extend(rows)
                reps_complete += 1
                rep_complete += 1

        # 10 reps x 2 scripts = 20 reps per cell at full completion
        if rep_complete >= 20:
            cells_complete += 1
        else:
            cells_incomplete += 1

    total_cells = 9
    print(f"[aggregate] cells complete: {cells_complete}/{total_cells} "
          f"({100*cells_complete/total_cells:.1f}%)")
    print(f"[aggregate] cells incomplete: {cells_incomplete}")
    print(f"[aggregate] reps complete: {reps_complete}")
    print(f"[aggregate] reps incomplete: {reps_incomplete}")
    if warnings:
        print(f"[aggregate] {len(warnings)} warning(s):")
        for w in warnings[:5]:
            print(f"    {w}")

    if not all_rows:
        sys.exit("[aggregate] no complete reps found.")

    df = pd.DataFrame(all_rows)
    df = df[["pl_mult", "mig", "rep", "script", "lineage",
             "fitness_end", "phenotype_end", "optimum_end", "fitness_lag",
             "extinct"]]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n[aggregate] wrote {args.out}  ({len(df)} rows)")

    # ---- Headline 3x3 rescue grid per lineage ----
    print("\n[aggregate] RESCUE EFFECT GRID  (extinct_rate_without - extinct_rate_with)")
    print("            positive = with-introgression rescues")
    print()

    for lin in EXPECTED_LINEAGES:
        sub = df[df["lineage"] == lin]
        if sub.empty:
            continue

        # Pivot to a wide pl x mig table for each script
        ext_with = (sub[sub["script"] == "with"]
                    .groupby(["pl_mult", "mig"])["extinct"].mean()
                    .unstack("mig"))
        ext_without = (sub[sub["script"] == "without"]
                       .groupby(["pl_mult", "mig"])["extinct"].mean()
                       .unstack("mig"))
        if ext_with.empty or ext_without.empty:
            continue
        rescue = ext_without - ext_with

        print(f"  --- {lin} ---")
        print(f"  rows = pl_mult (0.5 stronger sel, 2.0 weaker sel)")
        print(f"  cols = migration rate")
        print(rescue.to_string(float_format=lambda x: f"{x:+.2f}"))
        print()

    # ---- Overall summary line per lineage ----
    print("[aggregate] overall extinction rates per (lineage, script):")
    for lin in EXPECTED_LINEAGES:
        sub = df[df["lineage"] == lin]
        if sub.empty:
            continue
        with_r = sub[sub["script"] == "with"]["extinct"].mean()
        wo_r = sub[sub["script"] == "without"]["extinct"].mean()
        diff = wo_r - with_r
        marker = "  <-- BUFFERING" if diff > 0.05 else ("" if abs(diff) < 0.05 else "  <-- WORSE WITH GENE FLOW (investigate)")
        print(f"  {lin:10s}  with={with_r:.3f}  without={wo_r:.3f}  "
              f"rescue={diff:+.3f}{marker}")


if __name__ == "__main__":
    main()
