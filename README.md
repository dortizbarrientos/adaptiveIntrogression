# Adaptive introgression and climate buffering in *Populus*

We model the fitness consequences of the adaptive-haplotype introgression at *GA2OX8* reported in the accompanying empirical study. We are currently rebuilding the simulation along the lines of recent hybridising-bird work (Zhang et al. 2025, *Nature Climate Change*); this README describes both what we have today and what we are still building.

We ask one question: **does adaptive introgression measurably buffer population fitness as the climate warms?** We answer it with a polygenic fitness simulation that climate projections drive generation by generation, contrasting two demographically matched runs — one in which gene flow connects *P. rotundifolia* to northern *P. davidiana*, one in which we sever it. We hold the mechanistic *GA2OX8* relay simulation that raised this question in `99_reserve/`, per the lead authors' editorial decision, but we keep it available should reviewers ask.

---

## Status: rebuild in progress

| Component | State | Location |
|---|---|---|
| Per-generation climate trajectory builder | **We run this on real data** | `04_futureClimate/build_climate_trajectory.py` |
| SLiM polygenic model (with/without gene flow) | We are designing it, mirroring the [Fulvetta repo](https://github.com/willright28/Project-for-three-Fulvetta-species) | `04_futureClimate/` (pending) |
| Analytic harness for the new model | We are building it | `04_futureClimate/` (pending) |
| Buffering result + main figure | We will produce these from the harness | `figures/` (pending) |
| GA2OX8 relay simulation | We completed and reserved it | `99_reserve/` |
| Previous calibration-based fitness model | We superseded it with the current rebuild | `superseded/` |

We use **SSP585, 2081–2100** as the headline scenario throughout, with BIO10 (mean temperature of the warmest quarter) as the focal bioclimatic variable. We run other scenarios for sensitivity.

---

## Quick start

Our trajectory builder already runs end-to-end on the real data. From a clean clone:

```bash
git clone https://github.com/dortizbarrientos/adaptiveIntrogression.git
cd adaptiveIntrogression
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python 04_futureClimate/build_climate_trajectory.py \
    --data-dir 04_futureClimate/data_real \
    --out     04_futureClimate/data_real/populus_env.txt \
    --bioclim bio10 --ssp 585 --gen-time 5 --burn-in 1000
```

The script reads the real per-population bioclim values (current plus future across three GCMs and two SSPs), collapses 82 populations into three lineages (pdav-N, mixed, rot), and interpolates a per-generation climate trajectory at five-year generation steps. It writes a 3 × 1015-generation matrix in the format the SLiM model consumes, a diagnostic plot, and a metadata JSON.

We are adding the remaining components of the pipeline (SLiM model, analytic harness, figure assembly).

---

## How we organise the repository

```
adaptiveIntrogression/
├── 04_futureClimate/    the active rebuild
│   ├── build_climate_trajectory.py   builds per-generation climate from CMIP6 bioclim
│   ├── README_climate_trajectory.md  documents the trajectory builder
│   └── data_real/                    real Populus data (the lead authors kindly shared)
│       ├── INFO410.csv               410 individuals × population × lineage
│       ├── now_env.csv               current bioclim, 82 populations
│       ├── future_env/               future bioclim (3 GCMs × 2 SSPs × 2 periods)
│       ├── hap-freq.csv              GA2OX8 haplotype frequencies
│       ├── vulnerability/            genetic offsets per population
│       ├── recombination.csv         local recombination around GA2OX8
│       ├── target_variants.xlsx      15 climate-associated variants
│       ├── pda410-ga2ox8-20kb.recode.vcf   variant calls near the locus
│       ├── populus_env.txt           per-generation climate trajectory (we built this)
│       ├── populus_env.png           diagnostic plot
│       └── populus_env.json          build metadata
│
├── 99_reserve/          GA2OX8 mechanistic simulation (we hold this in reserve)
│   ├── 01_analytic/     validates against known theory
│   ├── 02_simulation/   simulates the SLiM cold-spot relay
│   ├── 03_analysis/     aggregates and builds the composite figure
│   ├── data/            grid end-state metrics (we commit these)
│   └── figures/         the relay + weight composite figure
│
├── superseded/          previous calibration-based fitness model
│   ├── 04_futureClimate/   the calibrated approach (we replaced this with exogenous CMIP6)
│   ├── figures/            its figures
│   └── ...                 retired plumbing scripts
│
├── docs/                README_pipeline.md, tutorial.html
├── figures/             rebuild figures (we will populate this)
├── tools/               repo-toolkit (we use this as an interactive command helper)
│
├── pipeline_lib.sh      shared stage definitions (carry-over; we will update for the rebuild)
├── run_all.sh           entry point (carry-over; we will update for the rebuild)
├── requirements.txt
├── LICENSE              MIT
└── CITATION.cff
```

---

## What we changed, and why

Three problems with our previous approach motivated this rebuild.

**We calibrated the optimum shift to the isolated lineage's genomic offset.** Because we defined the two lineages in part by their offset difference, the buffering result partly restated that difference rather than emerging as an independent prediction. In the rebuild, we drop the calibration entirely: we take the per-generation climate trajectory from CMIP6 projections, apply it identically to both lineages, and let buffering emerge as a prediction the model produces rather than as a target we tune it to reach.

**We anchored variation to a single locus.** Our previous model contrasted variance based on one locus's haplotype frequency. In the rebuild, we construct a polygenic trait from real climate-associated SNPs across the genome (each individual's breeding value is the sum across these loci), and we let each lineage's standing variation follow from its actual genotype matrix and a fitted demographic model.

**We compared abstractly.** Our previous comparison set with and without against an assumed variance contrast. In the rebuild, we run two demographically matched SLiM simulations that share the same climate trajectory and differ only in whether we activate the introgression channel between source and recipient. This matches the clean two-model contrast Zhang et al. use in the bird paper.

We mirror the architecture of the [Qu lab's Fulvetta repository](https://github.com/willright28/Project-for-three-Fulvetta-species), adapting it for the *Populus* three-lineage geometry (source, admixed bridge, northern recipient) rather than the bird paper's parapatric pair.

---

## Requirements

We use Python 3.10+ (we test the trajectory builder on 3.10–3.14). The trajectory builder depends on NumPy, pandas, Matplotlib, and SciPy (see `requirements.txt`). The polygenic model will require SLiM 4 (macOS: `brew install slim`).

## Data

The study's lead authors kindly shared the real data we use in `04_futureClimate/data_real/`. It includes the per-population genomic-vulnerability offsets, current and future bioclimatic variables (three GCMs × two SSPs × two periods), the *GA2OX8* haplotype frequencies, recombination estimates, and the variant calls around the locus. See `04_futureClimate/data_real/Data-list.docx` for the column-level description. The paper's data-availability statement determines where these data live relative to the study's primary archives (NGDC/NCBI).

## Reproducibility

The trajectory builder runs internal sanity checks and aborts with a clear message if any fails: we require uniquely mappable population assignments, every required bioclim column, and at least one population with complete bioclim data per lineage. We will apply the same standard to the SLiM model and analytic harness as they land.

We git-ignore regenerable outputs (`.venv/`, `__pycache__/`). We commit the trajectory builder, its outputs (`populus_env.{txt,png,json}`), and reference materials in `99_reserve/` and `superseded/`.

## License and citation

We release this under the MIT License (`LICENSE`). If you use this code, please cite the accompanying paper and this repository (`CITATION.cff`).
