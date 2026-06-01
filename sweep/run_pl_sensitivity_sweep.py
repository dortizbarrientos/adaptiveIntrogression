#!/usr/bin/env python3
"""
run_pl_sensitivity_sweep.py
===========================

Focused selection-strength sensitivity sweep, designed to run on the
M4 Pro MacBook while the main 5D sweep runs on a different machine.

The question
------------
How does the buffering effect (with vs without introgression) depend on
the strength of stabilising selection?

The grid: 3 x 3 = 9 cells
-------------------------
- pl_mult in {0.5, 1.0, 2.0}    -- multiplier on the empirical pl_vec
- migration m in {1e-6, 1e-5, 1e-4}   -- migration rate

Other parameters fixed at "headline" values:
  mutEffect = 0.05
  C = 50
  gen_time = 15  (Jing's stated typical value for Populus)
  Ne_scale = 1.0

At each cell, 10 replicates per script (with + without) = 180 total runs.
Estimated wall-time on M4 Pro at 10-way parallel: ~90 minutes.

Output: a 3x3 grid of buffering measures (extinction rates, fitness lag)
that can be plotted directly as a heatmap. The empirical Populus point
sits at (pl_mult=1.0, m=1e-5), in the centre of the grid.

Resumability and parallelism follow the same conventions as
run_pilot_sweep.py: DONE markers, ProcessPoolExecutor, per-run logs.

Usage
-----
python 04_futureClimate/run_pl_sensitivity_sweep.py              # full sweep
python 04_futureClimate/run_pl_sensitivity_sweep.py --dry-run    # plan only
python 04_futureClimate/run_pl_sensitivity_sweep.py --parallel 8 # adjust
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

# ============================================================================
# Grid configuration -- focused 3x3
# ============================================================================

PL_MULT_VALUES = [0.5, 1.0, 2.0]
MIG_VALUES = [1e-6, 1e-5, 1e-4]

# Fixed at headline values (the central cell of the M2 sweep)
MUT_EFFECT_FIXED = 0.05
C_FIXED = 50
GEN_TIME_FIXED = 15
NE_SCALE_FIXED = 1.0

N_REPS = 10

# Base seed chosen distinct from the M2 sweep so the two experiments
# are statistically independent (their seeds never coincide)
BASE_SEED = 22222222222222

TIMEOUT_SECONDS = 1800

# Base Ne values per lineage (matches the standard mig.txt and the M2 sweep)
NE_BASE = (50000, 30000, 20000)

# ============================================================================
# Paths
# ============================================================================

ENV_DATA_DIR = "04_futureClimate/data_real"
SCRIPT_WITH = "04_futureClimate/populus_with_introgression.slim"
SCRIPT_WITHOUT = "04_futureClimate/populus_without_introgression.slim"
TRAJECTORY_BUILDER = "04_futureClimate/build_climate_trajectory.py"


@dataclass
class RunSpec:
    cell_dir: Path
    rep_dir: Path
    script_path: str
    script_name: str
    mig_file: Path
    env_file: Path
    mig: float
    pl_mult: float
    rep: int
    seed: int

    def is_done(self) -> bool:
        return (self.rep_dir / "DONE").exists()

    def cell_label(self) -> str:
        return self.cell_dir.name


# ============================================================================
# File generators
# ============================================================================

def make_mig_file(mig_main: float, ne_scale: float, out_path: Path) -> None:
    """
    Generate cell-specific mig.txt; identical convention to the main sweep.

    Fields are tab-separated, matching the convention the SLiM script
    expects (strsplit with sep="\t"). Space-separation looks identical
    in `cat` output but causes the parser to read only the first column,
    leading to subscript-out-of-range crashes on row 2 of the matrix.
    """
    ne = [int(round(n * ne_scale)) for n in NE_BASE]
    mig_direct = mig_main / 10
    lines = [
        f"{ne[0]}\t{mig_main:g}\t{mig_direct:g}",
        f"{mig_main:g}\t{ne[1]}\t{mig_main:g}",
        f"{mig_direct:g}\t{mig_main:g}\t{ne[2]}",
    ]
    out_path.write_text("\n".join(lines) + "\n")


def ensure_trajectory(gen_time: int, sweep_dir: Path) -> Path:
    """Build or reuse the deviation-units trajectory for gen_time=15."""
    out_path = sweep_dir / f"populus_env_gentime{gen_time}.txt"
    if out_path.exists():
        return out_path

    print(f"[sweep] building trajectory for gen_time={gen_time}...")
    cmd = [
        sys.executable, TRAJECTORY_BUILDER,
        "--data-dir", ENV_DATA_DIR,
        "--out", str(out_path),
        "--bioclim", "bio10",
        "--ssp", "585",
        "--gen-time", str(gen_time),
        "--burn-in", "1000",
        "--units", "deviation",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"ERROR building trajectory:\n"
                 f"stdout: {result.stdout}\nstderr: {result.stderr}")
    return out_path


# ============================================================================
# Single-run execution
# ============================================================================

def run_one(spec_dict: dict) -> dict:
    """Execute one SLiM run, with DONE-marker resumability."""
    rep_dir = Path(spec_dict["rep_dir"])
    cell = spec_dict["cell"]
    script = spec_dict["script"]
    rep = spec_dict["rep"]

    if (rep_dir / "DONE").exists():
        return {"ok": True, "elapsed": 0.0, "rep_dir": str(rep_dir),
                "cell": cell, "script": script, "rep": rep, "err": None,
                "skipped": True}

    rep_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "slim",
        "-s", str(spec_dict["seed"]),
        "-d", "run='pl_sweep'",
        "-d", f"wkdir='{rep_dir}'",
        "-d", f"migMatrix='{spec_dict['mig_file']}'",
        "-d", f"Environment='{spec_dict['env_file']}'",
        "-d", "popscale=1",
        "-d", f"mutEffect={spec_dict['mutEffect']}",
        "-d", f"C={spec_dict['C']}",
        "-d", f"pl_mult={spec_dict['pl_mult']}",
        spec_dict["script_path"],
    ]

    start = time.time()
    try:
        with open(rep_dir / "stdout.log", "w") as out, \
             open(rep_dir / "stderr.log", "w") as err:
            result = subprocess.run(cmd, stdout=out, stderr=err,
                                     timeout=TIMEOUT_SECONDS)
        elapsed = time.time() - start
        if result.returncode == 0:
            (rep_dir / "DONE").touch()
            return {"ok": True, "elapsed": elapsed, "rep_dir": str(rep_dir),
                    "cell": cell, "script": script, "rep": rep, "err": None,
                    "skipped": False}
        return {"ok": False, "elapsed": elapsed, "rep_dir": str(rep_dir),
                "cell": cell, "script": script, "rep": rep,
                "err": f"exit code {result.returncode}", "skipped": False}
    except subprocess.TimeoutExpired:
        return {"ok": False, "elapsed": time.time() - start, "rep_dir": str(rep_dir),
                "cell": cell, "script": script, "rep": rep,
                "err": "TIMEOUT", "skipped": False}
    except Exception as e:
        return {"ok": False, "elapsed": time.time() - start, "rep_dir": str(rep_dir),
                "cell": cell, "script": script, "rep": rep,
                "err": str(e), "skipped": False}


# ============================================================================
# Job-list builder
# ============================================================================

def cell_label(pl_mult: float, mig: float) -> str:
    return f"pl-{pl_mult:g}__mig-{mig:g}"


def build_specs(sweep_dir: Path) -> tuple[list[RunSpec], Path]:
    env_file = ensure_trajectory(GEN_TIME_FIXED, sweep_dir)
    specs: list[RunSpec] = []
    cell_idx = 0
    for pl_mult in PL_MULT_VALUES:
        for mig in MIG_VALUES:
            label = cell_label(pl_mult, mig)
            cell_dir = sweep_dir / label
            cell_dir.mkdir(parents=True, exist_ok=True)

            mig_file = cell_dir / "mig.txt"
            make_mig_file(mig, NE_SCALE_FIXED, mig_file)

            for rep in range(1, N_REPS + 1):
                seed = BASE_SEED + cell_idx * 1000 + rep
                for script_name, script_path in [
                    ("with", SCRIPT_WITH),
                    ("without", SCRIPT_WITHOUT),
                ]:
                    rep_dir = cell_dir / f"{script_name}__rep-{rep:02d}"
                    specs.append(RunSpec(
                        cell_dir=cell_dir, rep_dir=rep_dir,
                        script_path=script_path, script_name=script_name,
                        mig_file=mig_file, env_file=env_file,
                        mig=mig, pl_mult=pl_mult, rep=rep, seed=seed,
                    ))
            cell_idx += 1
    return specs, env_file


def spec_to_dict(s: RunSpec) -> dict:
    return {
        "rep_dir": str(s.rep_dir),
        "cell": s.cell_label(),
        "script": s.script_name,
        "rep": s.rep,
        "seed": s.seed,
        "mig_file": str(s.mig_file),
        "env_file": str(s.env_file),
        "mutEffect": MUT_EFFECT_FIXED,
        "C": C_FIXED,
        "pl_mult": s.pl_mult,
        "script_path": s.script_path,
    }


# ============================================================================
# Main
# ============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep-dir", type=Path,
                    default=Path("04_futureClimate/output/pl_sensitivity"),
                    help="Output directory")
    ap.add_argument("--parallel", type=int, default=10,
                    help="Concurrent SLiM processes (default: 10 for M4 Pro)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for p in [SCRIPT_WITH, SCRIPT_WITHOUT, TRAJECTORY_BUILDER, ENV_DATA_DIR]:
        if not Path(p).exists():
            sys.exit(f"ERROR: required input not found: {p}\n"
                     f"Run from the repo root (~/adaptiveIntrogression).")

    args.sweep_dir.mkdir(parents=True, exist_ok=True)

    print(f"[sweep] enumerating jobs and building trajectory...")
    specs, env_file = build_specs(args.sweep_dir)
    total = len(specs)
    already_done = sum(1 for s in specs if s.is_done())
    remaining = total - already_done

    config = {
        "axes": {"pl_mult_values": PL_MULT_VALUES, "mig_values": MIG_VALUES},
        "fixed": {
            "mutEffect": MUT_EFFECT_FIXED,
            "C": C_FIXED,
            "gen_time": GEN_TIME_FIXED,
            "ne_scale": NE_SCALE_FIXED,
        },
        "n_reps": N_REPS,
        "base_seed": BASE_SEED,
        "cells": total // (N_REPS * 2),
        "total_runs": total,
        "parallel": args.parallel,
        "trajectory_file": str(env_file),
        "start_time": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(args.sweep_dir / "sweep_config.json", "w") as f:
        json.dump(config, f, indent=2)

    # Estimate ~4 minutes per run at popscale=1
    est_h = (remaining * 254 / args.parallel) / 3600
    print(f"[sweep] cells: 3 x 3 = 9 (pl_mult x migration)")
    print(f"[sweep] total runs: {total}")
    print(f"[sweep] already done: {already_done}")
    print(f"[sweep] remaining: {remaining}")
    print(f"[sweep] parallelism: {args.parallel}")
    print(f"[sweep] estimated wall-time: {est_h:.1f} hours")
    print(f"[sweep] output: {args.sweep_dir}")

    if args.dry_run:
        print(f"[sweep] DRY RUN -- not executing.")
        print(f"\n[sweep] first 5 jobs:")
        for s in specs[:5]:
            print(f"   {s.cell_label()}  script={s.script_name} rep={s.rep:02d} seed={s.seed}")
        return

    log_path = args.sweep_dir / "sweep.log"
    completed = already_done
    failed = 0
    timings = []
    sweep_start = time.time()

    with ProcessPoolExecutor(max_workers=args.parallel) as ex, \
         open(log_path, "a", buffering=1) as logf:

        def log(msg):
            ts = time.strftime("%Y-%m-%d %H:%M:%S")
            line = f"[{ts}] {msg}"
            print(line)
            logf.write(line + "\n")

        log(f"START total={total} done={already_done} remaining={remaining} "
            f"parallel={args.parallel}")

        pending = [s for s in specs if not s.is_done()]
        futures = {ex.submit(run_one, spec_to_dict(s)): s for s in pending}

        for fut in as_completed(futures):
            try:
                result = fut.result()
            except Exception as e:
                spec = futures[fut]
                log(f"FAIL_INTERNAL cell={spec.cell_label()} "
                    f"script={spec.script_name} rep={spec.rep} err={e}")
                failed += 1
                completed += 1
                continue

            completed += 1
            if result["ok"]:
                if not result.get("skipped"):
                    timings.append(result["elapsed"])
                avg = (sum(timings) / len(timings)) if timings else 0.0
                wall = time.time() - sweep_start
                done_in_session = completed - already_done
                eta_s = ((total - completed) * (wall / max(done_in_session, 1))
                         if done_in_session > 0 else 0)
                log(f"OK   {result['cell']} script={result['script']} "
                    f"rep={result['rep']:02d} t={result['elapsed']:.0f}s  "
                    f"progress={completed}/{total} ({100*completed/total:.1f}%) "
                    f"avg={avg:.0f}s eta={eta_s/60:.1f}min")
            else:
                failed += 1
                log(f"FAIL {result['cell']} script={result['script']} "
                    f"rep={result['rep']:02d} err={result['err']}")

        wall_min = (time.time() - sweep_start) / 60
        log(f"FINISHED completed={completed}/{total} failed={failed} "
            f"wall={wall_min:.1f}min")
        if failed > 0:
            log(f"NOTE: {failed} run(s) failed. Re-run to retry; "
                f"successful runs are skipped via DONE markers.")


if __name__ == "__main__":
    main()
