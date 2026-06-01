#!/usr/bin/env python3
# =============================================================================
# verify_pl_sensitivity.py  (v2 -- TSV/DONE layout)
#
# Post-run integrity check for the pl_sensitivity sweep.
#
# The sweep ran a 3 x 3 x 2 x 10 design (180 runs):
#     pl_mult   in {0.5, 1, 2}
#     migration in {1e-06, 1e-05, 1e-04}
#     script    in {with, without}      (with/without introgression)
#     rep       in {01 .. 10}
#
# Each run writes a directory:
#   pl_sensitivity/pl-{PL}__mig-{MIG}/{with|without}__rep-{NN}/
#       DONE                                            <- completion sentinel
#       populus_{w}_introgression_pl_sweep.tsv          <- per-generation trajectory
#       populus_{w}_introgression_pl_sweep_summary.tsv  <- one row per lineage
#       stdout.log
#       stderr.log
#
# "FINISHED failed=0" in the driver log says each PROCESS exited cleanly. This
# script checks the stronger claim: that every run left a complete, non-truncated
# set of files on disk. Three independent signals are used, weakest to strongest:
#
#   1. DONE sentinel present?         -> the run says it finished.
#   2. All five expected files there? -> nothing got dropped.
#   3. Trajectory TSV well-formed?    -> not truncated, right shape, parses.
#
# Design choices, stated plainly:
#   * PURE STDLIB. No tskit, no pandas -- runs in any python3, including the base
#     interpreter. (The previous version needed tskit and that is why Pass 2 was
#     skipped.)
#   * SELF-CALIBRATING. We do not hard-code how many generations or lineages a
#     trajectory "should" have. We read every file, take the MODAL (most common)
#     row count, max-generation, and lineage set as the expected values, then
#     flag any file that deviates. Truncation shows up as a short file against
#     the mode. This is robust to a design we don't have memorised.
#   * READ-ONLY. Never deletes or rewrites anything.
# =============================================================================

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

# -----------------------------------------------------------------------------
# EXPECTED GRID
# -----------------------------------------------------------------------------
PL_LEVELS  = ["0.5", "1", "2"]
MIG_LEVELS = ["1e-06", "1e-05", "0.0001"]
SCRIPTS    = ["with", "without"]
REPS       = [f"{i:02d}" for i in range(1, 11)]
EXPECTED_TOTAL = len(PL_LEVELS) * len(MIG_LEVELS) * len(SCRIPTS) * len(REPS)  # 180

SENTINEL     = "DONE"
TRAJ_TMPL    = "populus_{w}_introgression_pl_sweep.tsv"
SUMMARY_TMPL = "populus_{w}_introgression_pl_sweep_summary.tsv"
LOGS = ["stdout.log", "stderr.log"]

TRAJ_HEADER = ["replicate", "generation", "lineage", "n_individuals",
               "mean_fitness", "mean_phenotype", "sd_phenotype", "mean_optimum"]

# -----------------------------------------------------------------------------
# PATH PARSING  (run-dir segment is e.g. "with__rep-08" / "without__rep-04")
# -----------------------------------------------------------------------------
RE_PL     = re.compile(r"pl-([0-9.]+)__")
RE_MIG    = re.compile(r"mig-([0-9.eE+\-]+)")
RE_SCRIPT = re.compile(r"(?:^|/)(without|with)__rep")
RE_REP    = re.compile(r"rep-0*(\d+)")


def parse_coords(path_str: str):
    pl     = (m.group(1) if (m := RE_PL.search(path_str))     else None)
    mig    = (m.group(1) if (m := RE_MIG.search(path_str))    else None)
    script = (m.group(1) if (m := RE_SCRIPT.search(path_str)) else None)
    rep_m  = RE_REP.search(path_str)
    rep    = (f"{int(rep_m.group(1)):02d}" if rep_m else None)
    return pl, mig, script, rep


# -----------------------------------------------------------------------------
# DISCOVERY
# -----------------------------------------------------------------------------
def discover_runs(root: Path):
    run_dirs = {}
    for d in sorted(root.rglob("*__rep-*")):
        if not d.is_dir():
            continue
        coords = parse_coords(str(d))
        if None in coords:
            continue
        run_dirs[coords] = d
    return run_dirs


# -----------------------------------------------------------------------------
# PASS 1 -- COMPLETENESS
# -----------------------------------------------------------------------------
def pass1_completeness(root: Path, run_dirs: dict):
    expected_cells = set(product(PL_LEVELS, MIG_LEVELS, SCRIPTS, REPS))
    present_cells  = set(run_dirs.keys())
    missing_dirs   = expected_cells - present_cells
    extra_dirs     = present_cells - expected_cells

    done_ok, no_sentinel, incomplete = [], [], []
    for cell, d in run_dirs.items():
        w = cell[2]
        expected_files = [SENTINEL, TRAJ_TMPL.format(w=w),
                          SUMMARY_TMPL.format(w=w)] + LOGS
        missing_files = [f for f in expected_files if not (d / f).exists()]
        if missing_files:
            incomplete.append((cell, d, missing_files))
        if not (d / SENTINEL).exists():
            no_sentinel.append((cell, d))
        if not missing_files:
            done_ok.append((cell, d))

    print("=" * 74)
    print("PASS 1 - COMPLETENESS")
    print("=" * 74)
    print(f"  output root         : {root}")
    print(f"  expected runs       : {EXPECTED_TOTAL}")
    print(f"  run dirs found      : {len(run_dirs)}")
    print(f"  fully complete      : {len(done_ok)}   (DONE + both TSVs + both logs)")
    print(f"  missing DONE        : {len(no_sentinel)}")
    print(f"  missing some file   : {len(incomplete)}")
    print(f"  missing run dirs    : {len(missing_dirs)}")
    print(f"  unexpected run dirs : {len(extra_dirs)}")

    complete_cells = {c for c, _ in done_ok}
    print("\n  Complete runs per cell (want 10 in every column):")
    header = f"    {'pl':>4} {'mig':>8} | " + " | ".join(f"{s:>9}" for s in SCRIPTS)
    print(header)
    print("    " + "-" * (len(header) - 4))
    for pl in PL_LEVELS:
        for mig in MIG_LEVELS:
            cells = []
            for script in SCRIPTS:
                n = sum(1 for rep in REPS if (pl, mig, script, rep) in complete_cells)
                flag = "" if n == len(REPS) else "  <-- !"
                cells.append(f"{n:>9}{flag}")
            print(f"    {pl:>4} {mig:>8} | " + " | ".join(cells))

    def _dump(title, items, fmt):
        if items:
            print(f"\n  {title} ({len(items)}):")
            for it in items[:25]:
                print("    " + fmt(it))
            if len(items) > 25:
                print(f"    ... and {len(items) - 25} more")

    _dump("MISSING run dirs", sorted(missing_dirs),
          lambda c: f"pl-{c[0]}  mig-{c[1]}  {c[2]}  rep-{c[3]}")
    _dump("Missing DONE sentinel", no_sentinel, lambda t: f"{t[0]}   {t[1]}")
    _dump("Incomplete file set", incomplete, lambda t: f"{t[0]}  missing: {t[2]}")
    _dump("UNEXPECTED run dirs", sorted(extra_dirs), lambda c: f"{c}")

    return done_ok, missing_dirs, no_sentinel, incomplete, extra_dirs


# -----------------------------------------------------------------------------
# PASS 2 -- TSV INTEGRITY (self-calibrating)
# -----------------------------------------------------------------------------
def read_traj(path: Path):
    with path.open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def pass2_integrity(done_ok, sample):
    print("\n" + "=" * 74)
    print("PASS 2 - TSV INTEGRITY")
    print("=" * 74)

    if sample == 0:
        to_test = list(done_ok)
        mode = "ALL complete runs"
    else:
        by_stratum = defaultdict(list)
        for cell, d in done_ok:
            by_stratum[(cell[0], cell[1], cell[2])].append((cell, d))
        to_test = []
        for items in by_stratum.values():
            to_test.extend(sorted(items)[:sample])
        mode = f"first {sample} rep(s) per (pl,mig,script) stratum"

    print(f"  testing: {mode}  ->  {len(to_test)} trajectory file(s)\n")

    shapes, parse_fail = [], []
    for cell, d in to_test:
        w = cell[2]
        traj = d / TRAJ_TMPL.format(w=w)
        try:
            header, rows = read_traj(traj)
        except Exception as e:                       # noqa: BLE001
            parse_fail.append((cell, traj, repr(e)))
            continue
        header_ok = (header == TRAJ_HEADER)
        ncols = len(header)
        try:
            gens = [int(r[1]) for r in rows if len(r) > 1 and r[1].lstrip("-").isdigit()]
            max_gen = max(gens) if gens else -1
            lineages = tuple(sorted({r[2] for r in rows if len(r) > 2}))
            last_row_ok = bool(rows) and len(rows[-1]) == ncols
        except Exception:                            # noqa: BLE001
            max_gen, lineages, last_row_ok = -1, tuple(), False
        shapes.append([cell, traj, len(rows), max_gen, lineages, header_ok,
                       last_row_ok, ncols])

    if not shapes:
        print("  No trajectory files could be read.")
        return {"parse_fail": parse_fail, "deviations": [], "stderr_flags": []}

    mode_rows = Counter(s[2] for s in shapes).most_common(1)[0][0]
    mode_gen  = Counter(s[3] for s in shapes).most_common(1)[0][0]
    mode_lin  = Counter(s[4] for s in shapes).most_common(1)[0][0]
    print("  calibrated expectation (modal across files):")
    print(f"      rows/trajectory : {mode_rows}")
    print(f"      max generation  : {mode_gen}")
    print(f"      lineages        : {list(mode_lin)}  (n={len(mode_lin)})")
    if mode_gen > 0 and len(mode_lin) > 0:
        implied = mode_gen * len(mode_lin)
        tag = "matches" if implied == mode_rows else "DOES NOT match"
        print(f"      check: max_gen x n_lineages = {implied}  ({tag} modal rows)")
    print()

    deviations = []
    for cell, traj, n_rows, max_gen, lineages, header_ok, last_row_ok, ncols in shapes:
        reasons = []
        if n_rows != mode_rows:  reasons.append(f"rows={n_rows} (exp {mode_rows})")
        if max_gen != mode_gen:  reasons.append(f"max_gen={max_gen} (exp {mode_gen})")
        if lineages != mode_lin: reasons.append(f"lineages={list(lineages)}")
        if not header_ok:        reasons.append("header mismatch")
        if not last_row_ok:      reasons.append("ragged last row (truncated?)")
        if reasons:
            deviations.append((cell, traj, reasons))

    stderr_flags = []
    for cell, d in to_test:
        se = d / "stderr.log"
        if se.exists() and se.stat().st_size > 0:
            stderr_flags.append((cell, se, se.stat().st_size))

    n_clean = len(shapes) - len(deviations)
    print(f"  trajectories parsed   : {len(shapes)} / {len(to_test)}")
    print(f"  parse failures        : {len(parse_fail)}")
    print(f"  matching calibration  : {n_clean}")
    print(f"  deviating             : {len(deviations)}")
    print(f"  non-empty stderr.log  : {len(stderr_flags)}")

    if parse_fail:
        print("\n  PARSE FAILURES (corrupt/unreadable):")
        for cell, p, err in parse_fail:
            print(f"    {cell}\n        {p}\n        {err}")
    if deviations:
        print("\n  DEVIATIONS (shape differs from the modal trajectory):")
        for cell, p, reasons in deviations[:25]:
            print(f"    {cell}: {'; '.join(reasons)}")
            print(f"        {p}")
        if len(deviations) > 25:
            print(f"    ... and {len(deviations) - 25} more")
    if stderr_flags:
        print("\n  NON-EMPTY stderr.log (may be benign warnings -- worth a glance):")
        for cell, p, size in stderr_flags[:25]:
            print(f"    {size:>10,} B   {cell}   {p}")
        if len(stderr_flags) > 25:
            print(f"    ... and {len(stderr_flags) - 25} more")

    return {"parse_fail": parse_fail, "deviations": deviations,
            "stderr_flags": stderr_flags}


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Verify completeness/integrity of the pl_sensitivity sweep (TSV/DONE layout).")
    ap.add_argument("--root", type=str,
                    default="04_futureClimate/output/pl_sensitivity")
    ap.add_argument("--sample", type=int, default=2,
                    help="Trajectories per (pl,mig,script) stratum to open in Pass 2. 0 = all. Default 2.")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        print(f"ERROR: output root not found: {root}", file=sys.stderr)
        sys.exit(2)

    run_dirs = discover_runs(root)
    done_ok, missing_dirs, no_sentinel, incomplete, extra_dirs = \
        pass1_completeness(root, run_dirs)
    integ = pass2_integrity(done_ok, args.sample)

    print("\n" + "=" * 74)
    print("VERDICT")
    print("=" * 74)
    clean = (not missing_dirs and not no_sentinel and not incomplete
             and not extra_dirs and not integ["parse_fail"]
             and not integ["deviations"])
    if clean:
        print("  All 180 runs complete; sampled trajectories match calibration.")
        print("  Safe to proceed to analysis.")
        if integ["stderr_flags"]:
            print(f"  Note: {len(integ['stderr_flags'])} run(s) wrote to stderr -- check if benign.")
        sys.exit(0)
    else:
        print("  Issues found above. Inspect / re-run the affected cells before analysis.")
        sys.exit(1)


if __name__ == "__main__":
    main()
