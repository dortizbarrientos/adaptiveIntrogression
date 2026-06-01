# build_climate_trajectory.py

Builds the per-generation climate file consumed by the SLiM polygenic-fitness
simulation (NCC-style rebuild). Output format mirrors the Fulvetta
`env_example_data.txt` exactly: tab-separated, one row per lineage, one column
per generation.

## Quick start

Test it works on this machine, with mock data:

```bash
python build_climate_trajectory.py --mock --data-dir /tmp/populus_mock_data \
       --out /tmp/populus_env.txt --burn-in 50
```

Run on the real *Populus* data:

```bash
python build_climate_trajectory.py \
    --data-dir 04_futureClimate/data_real \
    --out    04_futureClimate/data_real/populus_env.txt \
    --bioclim bio10 \
    --ssp 585 \
    --gen-time 5 \
    --burn-in 1000
```

## What it does, in plain terms

1. Reads each population's present-day BIO10 (or whatever bioclim variable you
   pass to `--bioclim`).
2. Reads BIO10 under the future scenario you pick (`--ssp 585` for the severe
   case, `--ssp 245` for the moderate one), averaged across the three GCMs.
3. Collapses 82 populations into 3 lineages using `INFO410.csv`'s lineage
   assignments (`pdav-N`, `mixed`, `rot`).
4. For each lineage, linearly interpolates BIO10 from present-day (2020)
   through the 2070 midpoint to the 2090 midpoint, at one step per generation
   (5 years).
5. Prepends a burn-in period during which climate is stationary, so the SLiM
   model can reach equilibrium variance before warming starts.

## Expected input layout

```
04_futureClimate/data_real/
├── INFO410.csv                            # pop_id, lineage  (at minimum)
├── now_env.csv                            # pop_id, bio1, bio5, bio10, ...
└── future_env/
    ├── BCC-CSM2-MR_ssp245_2061-2080.csv
    ├── BCC-CSM2-MR_ssp245_2081-2100.csv
    ├── BCC-CSM2-MR_ssp585_2061-2080.csv
    ├── BCC-CSM2-MR_ssp585_2081-2100.csv
    ├── CMCC-ESM2_ssp585_2061-2080.csv
    └── ...                                # 12 files total (3 GCMs × 2 SSPs × 2 periods)
```

If Yupeng's filenames don't match this convention, the glob `*ssp{585}*{2081-2100}*.csv`
will fail; adjust the pattern in `average_future_bioclim()`.

## Outputs

- `populus_env.txt` — the SLiM input file
- `populus_env.png` — diagnostic trajectory plot (always inspect this)
- `populus_env.json` — metadata recording what was averaged

## Lineage order

The default order is `pdav-N, mixed, rot`, matching subpop 0/1/2 in the SLiM
model. If you flip the order, the SLiM model's subpop labels flip too.

## Design choices, made explicit

- **BIO10** (mean temperature of warmest quarter) is the default because Jing
  suggested it in her reply, and the bird paper used it for the same reason
  (warm-season heat stress drives mortality in the species studied).
- **Linear interpolation** between present, 2070, and 2090. The Fulvetta repo
  uses the same choice. A natural extension is to use the full transient
  CMIP6 trajectory if Yupeng can share it; this script's interpolation logic
  would change accordingly.
- **GCMs averaged**, not picked. Single-GCM trajectories are available with
  `--gcms BCC-CSM2-MR` etc. for sensitivity.
- **Lineage mean**, equal weighting per population. Sample-size weighting is
  a one-line change if it matters.

## Sanity checks built in

- If a population maps to multiple lineages, the script exits.
- If any lineage has zero populations with complete bioclim data, it exits.
- If the requested bioclim column is missing from a file, it exits and lists
  what's available.
- The diagnostic PNG should always show: (i) flat lines through the burn-in,
  (ii) monotonic warming after, (iii) plausible end values (the present plus
  ~3-4°C under SSP585 for the warmest scenario).
