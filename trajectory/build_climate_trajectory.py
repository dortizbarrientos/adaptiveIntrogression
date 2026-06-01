#!/usr/bin/env python3
"""
build_climate_trajectory.py
===========================

Build the per-generation climate trajectory file consumed by the SLiM
polygenic-fitness simulation (Populus rebuild, NCC-style architecture).

The output format mirrors the Fulvetta `env_example_data.txt`:
- Tab-separated, no header.
- One row per population (in fixed lineage order).
- One column per generation.
- Values are the climate variable (default BIO10, mean temperature of
  warmest quarter) experienced by that population in that generation.

Design choices, stated rather than buried:

1. **Which climate variable.**  BIO10 (mean T of warmest quarter) is the
   default because (a) it was Jing's specific suggestion in her reply, and
   (b) the bird paper used it for the same reasons (warm-season heat stress
   is the main physiological challenge under warming). Other BIOs can be
   specified on the command line.

2. **How the per-generation trajectory is built.**  We linearly interpolate
   between three anchor points: present-day (year 2020, the baseline of the
   bioclim data); the 2061-2080 midpoint (year 2070); and the 2081-2100
   midpoint (year 2090).  These are the only future bracket centres for which
   CMIP6 averages are available in the Populus dataset.  We then convert
   calendar years to generations using the user-supplied generation time
   (default 5 years for Populus).  Linear interpolation between centres is
   the same choice the Fulvetta `env_example_data.txt` makes ("linear
   climate change input used to simulate future warming" per their README).

3. **How we collapse 82 populations into 3 lineages.**  The Fulvetta model
   has three populations matching three species.  We have ~82 sampling
   points clustered into three lineages (north P. davidiana, admixed
   bridge, source P. rotundifolia) per the INFO410 assignments.  For each
   lineage we compute the *mean* of the chosen bioclim variable over its
   member populations, weighted by sample size if requested.  Each lineage
   then has a single trajectory.

4. **GCM and SSP handling.**  The default is to average across the three
   GCMs (BCC-CSM2-MR, CMCC-ESM2, MPI-ESM1-2-HR) for SSP585 (highest-emission
   scenario, matching the main-text headline).  Alternative SSPs and
   single-GCM trajectories can be requested for sensitivity analysis.

5. **Units: absolute vs deviation.**  The trajectory file can be written in
   either units. Default is `--units deviation`: each lineage's trajectory
   is expressed as deviation from its own present-day baseline (the mean
   value during the burn-in period). Burn-in values are 0; warming values
   rise. This matches the SLiM model's expectation that founders (which
   carry no mutations and therefore phenotype z ~ 0) sit near the optimum
   during burn-in. Use `--units absolute` for diagnostics if you want to
   see the raw temperatures.

Input files (expected in --data-dir, defaults match the Populus repo layout):
    now_env.csv             pop_id, bio1, bio2, ..., bio19   (present-day)
    future_env/<scenario>.csv   pop_id, bio1, ...           (one per SSPxperiodxGCM)
    INFO410.csv             individual_id, pop_id, lineage   (assignments)

Output:
    populus_env.txt         tab-separated, lineage x generation matrix
                            (3 rows, N_generations columns)
    populus_env_trajectory.png  diagnostic figure showing the three trajectories
    populus_env_metadata.json   what was averaged, generations covered, units, etc.

Usage:
    python build_climate_trajectory.py \\
        --data-dir 04_futureClimate/data_real/ \\
        --out 04_futureClimate/data_real/populus_env.txt \\
        --bioclim bio10 \\
        --ssp 585 \\
        --gen-time 5 \\
        --burn-in 1000 \\
        --units deviation  # default; use 'absolute' for diagnostics
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas required. Install with: pip install pandas")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("ERROR: matplotlib required. Install with: pip install matplotlib")


YEAR_PRESENT = 2020
YEAR_2070 = 2070
YEAR_2090 = 2090

DEFAULT_LINEAGE_ORDER = ["pdav-N", "mixed", "rot"]


@dataclass
class LineageTrajectory:
    lineage: str
    bioclim_variable: str
    present: float
    year2070: float
    year2090: float
    n_populations: int

    def value_at_year(self, year: float) -> float:
        if year <= YEAR_PRESENT:
            return self.present
        if year <= YEAR_2070:
            frac = (year - YEAR_PRESENT) / (YEAR_2070 - YEAR_PRESENT)
            return self.present + frac * (self.year2070 - self.present)
        if year <= YEAR_2090:
            frac = (year - YEAR_2070) / (YEAR_2090 - YEAR_2070)
            return self.year2070 + frac * (self.year2090 - self.year2070)
        return self.year2090


LINEAGE_LABEL_MAP = {
    "P. davidiana-N": "pdav-N",
    "admixed":        "mixed",
    "P. rotundifolia": "rot",
}


def load_population_assignments(info_csv: Path) -> pd.DataFrame:
    info = _read_csv_autosep(info_csv)
    if {"pop_id", "lineage"}.issubset(info.columns):
        pop_col, lineage_col = "pop_id", "lineage"
        translate = False
    elif {"pop", "group2"}.issubset(info.columns):
        pop_col, lineage_col = "pop", "group2"
        translate = True
    else:
        sys.exit(
            f"ERROR: {info_csv} must contain either "
            f"(`pop_id`, `lineage`) or (`pop`, `group2`). "
            f"Found: {list(info.columns)}"
        )
    pop_lineage = (
        info[[pop_col, lineage_col]]
        .rename(columns={pop_col: "pop_id", lineage_col: "lineage"})
        .drop_duplicates()
    )
    if translate:
        unknown = set(pop_lineage["lineage"]) - set(LINEAGE_LABEL_MAP)
        if unknown:
            sys.exit(
                f"ERROR: unrecognised lineage labels in {info_csv}: {sorted(unknown)}. "
                f"Expected one of {sorted(LINEAGE_LABEL_MAP)}. "
                f"Extend LINEAGE_LABEL_MAP at the top of this script if needed."
            )
        pop_lineage["lineage"] = pop_lineage["lineage"].map(LINEAGE_LABEL_MAP)
    multi = pop_lineage.groupby("pop_id")["lineage"].nunique()
    bad = multi[multi > 1]
    if len(bad) > 0:
        sys.exit(
            f"ERROR: {len(bad)} populations map to multiple lineages: "
            f"{bad.index.tolist()}"
        )
    print(f"[trajectory] loaded {len(pop_lineage)} populations across "
          f"{pop_lineage['lineage'].nunique()} lineages:")
    for lin, n in pop_lineage['lineage'].value_counts().sort_index().items():
        print(f"    {lin}: {n} populations")
    return pop_lineage.reset_index(drop=True)


def _read_csv_autosep(csv_path: Path) -> pd.DataFrame:
    return pd.read_csv(csv_path, sep=None, engine="python")


def load_bioclim_table(csv_path: Path, bioclim_col: str) -> pd.Series:
    df = _read_csv_autosep(csv_path)
    if "pop_id" in df.columns:
        pop_col = "pop_id"
    elif "pop" in df.columns:
        pop_col = "pop"
    else:
        sys.exit(
            f"ERROR: {csv_path} missing population identifier column. "
            f"Expected one of `pop_id` or `pop`. Found: {list(df.columns)}"
        )
    cols_lower = {c.lower(): c for c in df.columns}
    if bioclim_col.lower() not in cols_lower:
        sys.exit(
            f"ERROR: {csv_path} missing bioclim variable `{bioclim_col}`. "
            f"Available: {[c for c in df.columns if c.lower().startswith('bio')]}"
        )
    actual_col = cols_lower[bioclim_col.lower()]
    return df.set_index(pop_col)[actual_col]


def average_future_bioclim(
    future_dir: Path,
    period: str,
    ssp: str,
    bioclim_col: str,
    gcms: Iterable[str] | None = None,
) -> pd.Series:
    YEAR_FROM_PERIOD = {"2061-2080": "2080", "2081-2100": "2100"}
    year_token = YEAR_FROM_PERIOD.get(period, period)
    patterns = [
        f"*ssp{ssp}*{period}*.csv",
        f"*ssp{ssp}*{year_token}*.csv",
        f"*/*ssp{ssp}*{period}*.csv",
        f"*/*ssp{ssp}*{year_token}*.csv",
    ]
    matches: list[Path] = []
    for pat in patterns:
        matches.extend(future_dir.glob(pat))
    matches = sorted(set(matches))
    if not matches:
        all_csvs = sorted(future_dir.rglob("*.csv"))
        sample = [str(p.relative_to(future_dir)) for p in all_csvs[:8]]
        sys.exit(
            f"ERROR: no files match SSP={ssp} period={period} (year_token={year_token}) "
            f"under {future_dir}. Sample of files present: {sample}"
        )
    if gcms is not None:
        matches = [m for m in matches if any(g in m.name for g in gcms)]
        if not matches:
            sys.exit(f"ERROR: no files match SSP={ssp} period={period} gcms={gcms}")
    print(f"[trajectory] SSP{ssp} {period}: averaging across {len(matches)} GCMs:")
    for m in matches:
        print(f"    {m.relative_to(future_dir)}")
    series_list = [load_bioclim_table(m, bioclim_col) for m in matches]
    return pd.concat(series_list, axis=1).mean(axis=1)


def build_lineage_trajectories(
    pop_lineage: pd.DataFrame,
    present_values: pd.Series,
    values_2070: pd.Series,
    values_2090: pd.Series,
    bioclim_col: str,
    lineage_order: list[str],
    weight_by_sample_size: bool = False,
) -> list[LineageTrajectory]:
    trajectories: list[LineageTrajectory] = []
    common_pops = (
        set(present_values.index) & set(values_2070.index) & set(values_2090.index)
    )
    print(f"[trajectory] populations present in all three time points: {len(common_pops)}")
    for lineage in lineage_order:
        lin_pops = set(pop_lineage[pop_lineage["lineage"] == lineage]["pop_id"])
        usable = lin_pops & common_pops
        if not usable:
            sys.exit(
                f"ERROR: lineage `{lineage}` has no populations with full bioclim data. "
                f"Assigned populations: {sorted(lin_pops)[:5]}..."
            )
        idx = list(usable)
        present_mean = present_values.loc[idx].mean()
        y2070_mean = values_2070.loc[idx].mean()
        y2090_mean = values_2090.loc[idx].mean()
        traj = LineageTrajectory(
            lineage=lineage,
            bioclim_variable=bioclim_col,
            present=float(present_mean),
            year2070=float(y2070_mean),
            year2090=float(y2090_mean),
            n_populations=len(usable),
        )
        trajectories.append(traj)
        print(
            f"[trajectory] {lineage}: n={traj.n_populations}, "
            f"present={traj.present:.2f}, 2070={traj.year2070:.2f}, "
            f"2090={traj.year2090:.2f}  (delta={traj.year2090-traj.present:+.2f})"
        )
    return trajectories


def expand_to_generation_matrix(
    trajectories: list[LineageTrajectory],
    gen_time: float,
    burn_in_generations: int,
    year_start: int = YEAR_PRESENT,
    year_end: int = YEAR_2090,
) -> tuple[np.ndarray, list[int]]:
    n_warming = int(np.ceil((year_end - year_start) / gen_time))
    total_gens = burn_in_generations + n_warming + 1
    matrix = np.zeros((len(trajectories), total_gens), dtype=float)
    for t_idx, traj in enumerate(trajectories):
        matrix[t_idx, :burn_in_generations] = traj.present
        for g in range(n_warming + 1):
            year = year_start + g * gen_time
            matrix[t_idx, burn_in_generations + g] = traj.value_at_year(year)
    gen_indices = list(range(total_gens))
    print(
        f"[trajectory] matrix shape: {matrix.shape}  "
        f"(burn-in: {burn_in_generations} gens, warming+plateau: {n_warming + 1} gens)"
    )
    return matrix, gen_indices


def convert_to_deviation(
    matrix: np.ndarray,
    burn_in_generations: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert an absolute-units climate trajectory matrix into deviations
    from each lineage's burn-in baseline.

    Why we do this. The SLiM polygenic-fitness model treats founder
    individuals as carrying no mutations (mutations arise during the run),
    so their initial phenotype z is ~ 0 (just environmental noise). If we
    feed the model a trajectory in absolute degrees C -- where the baseline
    sits at ~19 C for pdav-N -- the founders are 19 trait-units away from
    the optimum and fitness collapses to the cushion value. By centring
    each lineage's trajectory on its own present-day baseline, founders at
    z=0 match theta ~ 0 during burn-in, fitness can stabilise near 1, and
    warming becomes the only signal that perturbs the system.

    The baseline for each lineage is the mean trajectory value across its
    burn-in window. Because the burn-in is held constant in
    `expand_to_generation_matrix`, this is exactly the lineage's
    present-day value, but we compute it from the matrix to remain robust
    to any future changes in how burn-in is constructed.

    Parameters
    ----------
    matrix : np.ndarray
        Absolute-units matrix, shape (n_lineages, n_generations).
    burn_in_generations : int
        Number of leading columns considered the equilibrium baseline.

    Returns
    -------
    deviation_matrix : np.ndarray
        Same shape as `matrix`, but each row centred on its burn-in mean.
        Burn-in columns are 0; warming columns rise.
    baselines : np.ndarray
        Per-lineage baselines (the values subtracted), shape (n_lineages,).
        Recorded in the metadata for reversibility.
    """
    if burn_in_generations <= 0:
        sys.exit(
            "ERROR: deviation units require a positive burn-in window "
            "to compute a baseline. Use --burn-in N (N > 0)."
        )
    if burn_in_generations > matrix.shape[1]:
        sys.exit(
            f"ERROR: burn-in length ({burn_in_generations}) exceeds matrix "
            f"length ({matrix.shape[1]}). Cannot compute baseline."
        )
    baselines = matrix[:, :burn_in_generations].mean(axis=1)
    deviation_matrix = matrix - baselines[:, np.newaxis]
    bi_means = deviation_matrix[:, :burn_in_generations].mean(axis=1)
    if not np.allclose(bi_means, 0.0, atol=1e-10):
        sys.exit(
            f"ERROR: deviation conversion failed sanity check. "
            f"Burn-in means after centring: {bi_means}"
        )
    print(f"[trajectory] converted to deviation units; baselines subtracted: "
          f"{[f'{b:.3f}' for b in baselines]}")
    return deviation_matrix, baselines


def write_slim_env_file(matrix: np.ndarray, out_path: Path) -> None:
    np.savetxt(out_path, matrix, fmt="%.6f", delimiter="\t")
    print(f"[trajectory] wrote {out_path}  ({matrix.shape[0]} rows, {matrix.shape[1]} cols)")


def plot_trajectories(
    matrix: np.ndarray,
    lineage_order: list[str],
    bioclim_col: str,
    burn_in_generations: int,
    units: str,
    out_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    gens = np.arange(matrix.shape[1])
    colors = {"pdav-N": "#C2566E", "mixed": "#6E94C4", "rot": "#574071"}
    for i, lineage in enumerate(lineage_order):
        ax.plot(
            gens,
            matrix[i, :],
            label=lineage,
            color=colors.get(lineage, f"C{i}"),
            lw=2,
        )
    if burn_in_generations > 0:
        ax.axvspan(0, burn_in_generations, color="grey", alpha=0.08, label="burn-in")
    if units == "deviation":
        ax.axhline(0, color="black", lw=0.5, ls=":")
    ax.set_xlabel("generation")
    if units == "deviation":
        ax.set_ylabel(f"{bioclim_col} deviation from present-day (degC)")
    else:
        ax.set_ylabel(f"{bioclim_col} (lineage mean, absolute)")
    ax.set_title(f"Per-lineage climate trajectory ({units} units, input to SLiM)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[trajectory] wrote diagnostic plot {out_path}")


def write_metadata(metadata: dict, out_path: Path) -> None:
    with open(out_path, "w") as f:
        json.dump(metadata, f, indent=2, default=str)
    print(f"[trajectory] wrote metadata {out_path}")


def create_mock_data(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "future_env").mkdir(exist_ok=True)
    rng = np.random.default_rng(42)
    pops = [f"P{i:02d}" for i in range(1, 10)]
    lineages = (["pdav-N"] * 3) + (["mixed"] * 3) + (["rot"] * 3)
    info = pd.DataFrame({"pop_id": pops, "lineage": lineages, "individual_id": pops})
    info.to_csv(out_dir / "INFO410.csv", index=False)
    present_bio10 = 15.0 + rng.normal(0, 1.5, len(pops))
    now = pd.DataFrame({"pop_id": pops, "bio10": present_bio10})
    for b in ["bio1", "bio5", "bio12"]:
        now[b] = rng.normal(10, 3, len(pops))
    now.to_csv(out_dir / "now_env.csv", index=False)
    for gcm in ["BCC-CSM2-MR", "CMCC-ESM2", "MPI-ESM1-2-HR"]:
        for ssp, warming_2090 in [("245", 2.5), ("585", 3.5)]:
            for period, frac in [("2061-2080", 0.6), ("2081-2100", 1.0)]:
                bio10 = present_bio10 + warming_2090 * frac + rng.normal(0, 0.3, len(pops))
                fut = pd.DataFrame({"pop_id": pops, "bio10": bio10})
                for b in ["bio1", "bio5", "bio12"]:
                    fut[b] = rng.normal(10, 3, len(pops))
                fut.to_csv(
                    out_dir / "future_env" / f"{gcm}_ssp{ssp}_{period}.csv", index=False
                )
    print(f"[trajectory] wrote mock data to {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build per-generation climate trajectory for SLiM fitness model"
    )
    ap.add_argument("--data-dir", type=Path, required=False,
                    help="Directory containing now_env.csv, future_env/, INFO410.csv")
    ap.add_argument("--out", type=Path, default=Path("populus_env.txt"),
                    help="Output env file (Fulvetta format)")
    ap.add_argument("--bioclim", default="bio10",
                    help="Bioclim variable (default: bio10, mean T warmest quarter)")
    ap.add_argument("--ssp", choices=["245", "585"], default="585",
                    help="Shared Socioeconomic Pathway (default 585)")
    ap.add_argument("--gen-time", type=float, default=5.0,
                    help="Generation time in years (default 5, for Populus)")
    ap.add_argument("--burn-in", type=int, default=1000,
                    help="Number of equilibrium generations before warming")
    ap.add_argument("--gcms", nargs="*", default=None,
                    help="Restrict to specific GCMs (default: all available)")
    ap.add_argument("--lineage-order", nargs=3, default=DEFAULT_LINEAGE_ORDER,
                    help=f"Lineage order for SLiM subpop rows "
                         f"(default {DEFAULT_LINEAGE_ORDER})")
    ap.add_argument("--units", choices=["deviation", "absolute"], default="deviation",
                    help="Trajectory units. 'deviation' (default) centres each "
                         "lineage on its present-day baseline so burn-in values "
                         "are 0; 'absolute' keeps raw temperatures (diagnostic only).")
    ap.add_argument("--mock", action="store_true",
                    help="Generate mock data in --data-dir and use it (for testing)")
    args = ap.parse_args()

    if args.mock:
        if args.data_dir is None:
            args.data_dir = Path("/tmp/populus_mock_data")
        create_mock_data(args.data_dir)

    if args.data_dir is None or not args.data_dir.exists():
        sys.exit(f"ERROR: --data-dir required and must exist (got {args.data_dir})")

    info_path = args.data_dir / "INFO410.csv"
    pop_lineage = load_population_assignments(info_path)

    present = load_bioclim_table(args.data_dir / "now_env.csv", args.bioclim)

    future_dir = args.data_dir / "future_env"
    if not future_dir.exists():
        sys.exit(f"ERROR: {future_dir} not found")
    values_2070 = average_future_bioclim(future_dir, "2061-2080", args.ssp,
                                          args.bioclim, args.gcms)
    values_2090 = average_future_bioclim(future_dir, "2081-2100", args.ssp,
                                          args.bioclim, args.gcms)

    trajectories = build_lineage_trajectories(
        pop_lineage, present, values_2070, values_2090,
        args.bioclim, args.lineage_order,
    )

    matrix_abs, _ = expand_to_generation_matrix(
        trajectories, args.gen_time, args.burn_in,
    )

    if args.units == "deviation":
        matrix, baselines = convert_to_deviation(matrix_abs, args.burn_in)
        baselines_record = baselines.tolist()
    else:
        matrix = matrix_abs
        baselines_record = None

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_slim_env_file(matrix, args.out)
    plot_trajectories(
        matrix, args.lineage_order, args.bioclim, args.burn_in, args.units,
        args.out.with_suffix(".png"),
    )
    write_metadata(
        {
            "bioclim": args.bioclim,
            "ssp": args.ssp,
            "gen_time_years": args.gen_time,
            "burn_in_generations": args.burn_in,
            "units": args.units,
            "baselines_subtracted": baselines_record,
            "year_anchors": {"present": YEAR_PRESENT,
                             "2070": YEAR_2070, "2090": YEAR_2090},
            "lineage_order": args.lineage_order,
            "lineage_anchors": [
                {"lineage": t.lineage,
                 "n_populations": t.n_populations,
                 "present": t.present, "y2070": t.year2070, "y2090": t.year2090,
                 "delta_2090": t.year2090 - t.present}
                for t in trajectories
            ],
            "matrix_shape": list(matrix.shape),
        },
        args.out.with_suffix(".json"),
    )

    print("\n[trajectory] DONE")


if __name__ == "__main__":
    main()
