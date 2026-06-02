#!/usr/bin/env python3
# =============================================================================
# aggregate_pl_sensitivity.py
#
# Collapse the 180-run pl_sensitivity sweep into ONE tidy frame, one row per
# (pl, mig, script, rep, lineage), ready for plotting and stats.
#
# WHY TWO FILES PER RUN. The summary TSV carries fitness_end and fitness_lag but
# NOT n_individuals -- it is blind to demography. Final and founding population
# size live only in the trajectory TSV (last and first generation, column
# n_individuals). To carry both a maladaptation measure (lag) and a persistence
# measure (final n / founding n), we must read both files and join them.
#
# WHAT WE COMPUTE per (run x lineage):
#   fitness_lag    : from summary, as-is (optimum - phenotype, in the SLiM units)
#   fitness_end    : from summary, as-is (mean fitness at the final generation)
#   n_founding     : trajectory, n_individuals at the FIRST generation
#   n_final        : trajectory, n_individuals at the LAST generation
#   persistence    : n_final / n_founding, clipped to [0, 1]  <- main response
#
# DESIGN:
#   * Self-describing: founding size is READ from gen 1 of each run, never
#     hard-coded, so a lineage that founds at a different N is handled correctly.
#   * Pure stdlib + pandas only for the final frame. Trajectory reads pull just
#     the first and last data rows -- we do not hold 3018 rows x 180 files in RAM.
#   * Writes a single tidy CSV. Stats and plots consume that, nothing re-walks.
# =============================================================================

from __future__ import annotations

import argparse
import csv
import re
import sys
from itertools import product
from pathlib import Path

import pandas as pd

PL_LEVELS  = ["0.5", "1", "2"]
MIG_LEVELS = ["1e-06", "1e-05", "0.0001"]
SCRIPTS    = ["with", "without"]
REPS       = [f"{i:02d}" for i in range(1, 11)]

TRAJ_TMPL    = "populus_{w}_introgression_pl_sweep.tsv"
SUMMARY_TMPL = "populus_{w}_introgression_pl_sweep_summary.tsv"

RE_PL     = re.compile(r"pl-([0-9.]+)__")
RE_MIG    = re.compile(r"mig-([0-9.eE+\-]+)")
RE_SCRIPT = re.compile(r"(?:^|/)(without|with)__rep")
RE_REP    = re.compile(r"rep-0*(\d+)")

# Column indices in the trajectory TSV (0-based), per the verified header:
#   0 replicate | 1 generation | 2 lineage | 3 n_individuals | ...
GEN_COL, LIN_COL, N_COL = 1, 2, 3


def parse_coords(path_str: str):
    pl  = (m.group(1) if (m := RE_PL.search(path_str))  else None)
    mig = (m.group(1) if (m := RE_MIG.search(path_str)) else None)
    scr = (m.group(1) if (m := RE_SCRIPT.search(path_str)) else None)
    rm  = RE_REP.search(path_str)
    rep = (f"{int(rm.group(1)):02d}" if rm else None)
    return pl, mig, scr, rep


def traj_founding_and_final(path: Path):
    """Return {lineage: (n_founding, n_final)} by reading only the first and last
    generation blocks. We stream the file once, tracking the min and max
    generation seen per lineage and the n at each -- cheap and RAM-flat."""
    first_gen, last_gen = {}, {}
    first_n, last_n = {}, {}
    with path.open(newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r, None)  # header
        for row in r:
            if len(row) <= N_COL:
                continue
            try:
                g = int(row[GEN_COL]); n = int(row[N_COL])
            except ValueError:
                continue
            ln = row[LIN_COL]
            if ln not in first_gen or g < first_gen[ln]:
                first_gen[ln] = g; first_n[ln] = n
            if ln not in last_gen or g > last_gen[ln]:
                last_gen[ln] = g; last_n[ln] = n
    return {ln: (first_n[ln], last_n[ln]) for ln in first_gen}


def read_summary(path: Path):
    """Return {lineage: {fitness_end, fitness_lag}}."""
    out = {}
    with path.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            out[row["lineage"]] = {
                "fitness_end": float(row["fitness_end"]),
                "fitness_lag": float(row["fitness_lag"]),
            }
    return out


def aggregate(root: Path):
    records, problems = [], []
    expected = set(product(PL_LEVELS, MIG_LEVELS, SCRIPTS, REPS))
    seen = set()

    for d in sorted(root.rglob("*__rep-*")):
        if not d.is_dir():
            continue
        pl, mig, scr, rep = parse_coords(str(d))
        if None in (pl, mig, scr, rep):
            continue
        seen.add((pl, mig, scr, rep))
        w = scr
        traj = d / TRAJ_TMPL.format(w=w)
        summ = d / SUMMARY_TMPL.format(w=w)
        if not traj.exists() or not summ.exists():
            problems.append((d, "missing traj or summary")); continue
        try:
            demo = traj_founding_and_final(traj)
            fit  = read_summary(summ)
        except Exception as e:                       # noqa: BLE001
            problems.append((d, repr(e))); continue

        for ln in fit:
            if ln not in demo:
                problems.append((d, f"lineage {ln} in summary not in trajectory")); continue
            n0, nfin = demo[ln]
            persistence = (nfin / n0) if n0 > 0 else float("nan")
            persistence = min(max(persistence, 0.0), 1.0)
            records.append({
                "pl": pl, "mig": mig, "script": scr, "rep": rep, "lineage": ln,
                "n_founding": n0, "n_final": nfin, "persistence": persistence,
                "fitness_end": fit[ln]["fitness_end"],
                "fitness_lag": fit[ln]["fitness_lag"],
            })

    df = pd.DataFrame.from_records(records)
    # Order the categorical axes meaningfully (not alphabetically).
    if not df.empty:
        df["pl"]     = pd.Categorical(df["pl"], categories=PL_LEVELS, ordered=True)
        df["mig"]    = pd.Categorical(df["mig"], categories=MIG_LEVELS, ordered=True)
        df["script"] = pd.Categorical(df["script"], categories=SCRIPTS, ordered=True)
        df["lineage"] = pd.Categorical(
            df["lineage"], categories=["pdav-N", "mixed", "rot"], ordered=True)
        df = df.sort_values(["pl", "mig", "script", "lineage", "rep"]).reset_index(drop=True)

    missing = expected - seen
    return df, problems, missing


def main():
    ap = argparse.ArgumentParser(description="Aggregate pl_sensitivity sweep into a tidy frame.")
    ap.add_argument("--root", default="output/pl_sensitivity")
    ap.add_argument("--out", default="pl_sensitivity_tidy.csv")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        print(f"ERROR: root not found: {root}", file=sys.stderr); sys.exit(2)

    df, problems, missing = aggregate(root)

    print("=" * 70)
    print("AGGREGATION SUMMARY")
    print("=" * 70)
    print(f"  root            : {root}")
    print(f"  rows (run x lineage) : {len(df)}")
    if not df.empty:
        n_runs = df.groupby(['pl','mig','script','rep'], observed=True).ngroups
        print(f"  distinct runs   : {n_runs}  (expect 180)")
        print(f"  lineages        : {list(df['lineage'].cat.categories)}")
        print(f"  rows/run        : {len(df)//max(n_runs,1)} (expect 3)")
    print(f"  problems        : {len(problems)}")
    print(f"  missing runs    : {len(missing)}")
    for d, why in problems[:10]:
        print(f"     ! {d}: {why}")
    for c in sorted(missing)[:10]:
        print(f"     missing: {c}")

    df.to_csv(args.out, index=False)
    print(f"\n  wrote tidy frame -> {args.out}")

    # A small at-a-glance pivot: mean persistence by (lineage x script), pooled.
    if not df.empty:
        print("\n  mean persistence by lineage x script (pooled over pl, mig, rep):")
        piv = df.pivot_table(index="lineage", columns="script",
                             values="persistence", aggfunc="mean", observed=True)
        print(piv.to_string(float_format=lambda x: f"{x:6.3f}"))


if __name__ == "__main__":
    main()
