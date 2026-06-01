# SLiM polygenic-fitness model — first draft

Two SLiM 4 scripts implementing the with/without-introgression contrast for
the *Populus* fitness rebuild. Mirrors the architecture of the
[Fulvetta repository](https://github.com/willright28/Project-for-three-Fulvetta-species)
(Zhang et al. 2025, *Nature Climate Change*), adapted for the *Populus*
three-lineage geometry.

## The contrast

| Script | What it does | Differs from the other in... |
|---|---|---|
| `populus_with_introgression.slim` | Honours all migration rates from `mig.txt` throughout | (the reference) |
| `populus_without_introgression.slim` | Zeroes migration INTO pdav-N at start of warming | exactly one inserted block + output filenames |

The two scripts are byte-identical apart from these explicit differences.
Any difference in their outputs is therefore attributable to the
introgression channel and nothing else.

## What's pending before you can run

Two inputs needed before SLiM can be invoked:

1. **`mig.txt`** — migration matrix. Diagonal = Ne per lineage, off-diagonal
   = per-generation migration rates. Format identical to the Fulvetta
   `mig.txt`. For *Populus*, plausible literature placeholders are below.

2. **`populus_env.txt`** — the per-generation climate trajectory. Already built
   on real data via `build_climate_trajectory.py`; lives in
   `04_futureClimate/data_real/populus_env.txt`.

A literature-justified `mig.txt` for *Populus* (placeholder pending Yupeng):
```
50000   1e-5    1e-6
1e-5    30000   1e-5
1e-6    1e-5    20000
```
Row/column order: **pdav-N, mixed, rot**. Ne values are mid-range from
*Populus* literature (Ingvarsson 2008 et seq.). Migration rates are
deliberately modest — strong enough to matter, weak enough that the
introgression effect must accumulate over many generations. **These swap
out when Yupeng's fastsimcoal2 fit arrives.**

## Default parameter values

Embedded in both scripts as `defineConstant` defaults, all overridable via
`-d` on the command line:

| Parameter | Default | Notes |
|---|---|---|
| `C` (number of climate-adaptive loci) | 50 | Sweep range: 20–200 |
| `mutEffect` (SD of effect size) | 0.05 | Sweep range: 0.01–0.20 |
| `mu` (mutation rate) | 5e-9 | Mid-range *Populus* literature (2.5e-9 to 7e-9) |
| `pl` (selection-function width) | 2.5 | Roughly the SD of BIO10 across *Populus* populations |
| `esd` (environmental noise SD) | 0.1 | Modest relative to mutEffect |
| `warming_start` (gen warming begins) | 1000 | Matches the burn-in length in `populus_env.txt` |

## Running a single replicate

From the repo root, with `mig.txt` and `populus_env.txt` in place:

```bash
mkdir -p 04_futureClimate/output

slim -d "run='r001'" \
     -d "wkdir='04_futureClimate/output'" \
     -d "migMatrix='04_futureClimate/data_real/mig.txt'" \
     -d "Environment='04_futureClimate/data_real/populus_env.txt'" \
     04_futureClimate/populus_with_introgression.slim
```

The same invocation for `populus_without_introgression.slim` produces the
companion run. Output files live in `04_futureClimate/output/`:

- `populus_with_introgression_r001.tsv` — per-generation per-lineage data
- `populus_with_introgression_r001_summary.tsv` — end-of-warming summary

## Per-generation output

One row per (replicate, generation, lineage):

```
replicate  generation  lineage  n_individuals  mean_fitness  mean_phenotype  sd_phenotype  mean_optimum
```

`mean_optimum` is included so the climate trajectory can be verified after
the fact without reading the SLiM script.

## End-of-warming summary

One row per (replicate, lineage) at the final generation:

```
replicate  lineage  fitness_end  phenotype_end  optimum_end  fitness_lag
```

`fitness_lag` is `optimum_end - phenotype_end` — how far behind the moving
target the lineage finished. This is the diagnostic that reveals *why* one
lineage buffers more than another, not just *that* it does.

## What's NOT in this first draft

- **Real GEA loadings.** Effect sizes are drawn from `N(0, mutEffect)` per
  Fulvetta default. When Yupeng sends per-SNP loadings, a one-function
  change will swap these in as Option B (sensitivity analysis).
- **Pipeline integration.** The driver shell script and the analytic harness
  (Item C) come next.
- **A multi-replicate runner.** Running N replicates with seeded RNG and
  aggregating outputs is one shell loop away; deferred until the single
  replicate is validated.

## Limiting cases worth checking (Item C deliverables)

These are the standard analytic anchors the harness will validate against:

- **No climate change** (set `populus_env.txt` constant): fitness should
  stay near 1.0 throughout in both scripts.
- **No migration in `mig.txt`** (all off-diagonal entries zero): with- and
  without-introgression scripts should give identical results.
- **Large Ne limit**: stochastic SLiM should track the deterministic
  breeder-equation expectation.

If the SLiM scripts fail any of these, we know to debug before trusting
the headline buffering result.
