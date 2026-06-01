# Cleaning the tree and publishing to GitHub

Two phases: prune the working tree to what matters, then push it as a repository.
Run everything from the repository root (the folder that contains `run_all.sh`).

## 0. Drop in the scaffolding

Copy these files into the repository root:

```
.gitignore   LICENSE   CITATION.cff   README.md   cleanup.sh
```

Replace `docs/supplementary.tex` with the repo-relative copy provided here (its
figure and code-listing paths now point at `../figures/` and the sibling source
folders, so it compiles from inside `docs/`).

## 1. Prune (safe: dry run first)

```bash
bash cleanup.sh            # DRY RUN — shows what it would do, changes nothing
bash cleanup.sh --apply    # do it: lifts the aggregated CSV to the root, moves the
                           # plane/validation figures into figures/, and removes the
                           # 216 raw SLiM logs + phase plots + summaries (regenerable)
```

Use `--apply --archive` instead if you would rather tar the raw grid to
`results_grid_raw.tar.gz` (gitignored) than delete it. The raw grid is fully
reproducible with `bash run_all.sh sim aggregate`, so deleting it is safe; only the
small `grid_endstate_metrics.csv` is needed to rebuild every figure.

## 2. Build the composite figure (so figures/ is complete)

```bash
bash run_all.sh figure     # writes figures/figure_results.{png,pdf} from the committed CSV
```

## 3. Initialise git and make the first commit

```bash
git init
git add .
git status                 # confirm results_grid/ and caches are NOT staged
git commit -m "Initial commit: GA2OX8 adaptive-introgression modelling pipeline"
```

## 4. Create the GitHub repository and push

**With the GitHub CLI (`gh`)** — one command:

```bash
gh repo create populus_pipeline --public --source=. --remote=origin --push
```

**Without `gh`** — create an empty repo named `populus_pipeline` on github.com
(no README/license, since you already have them), then:

```bash
git branch -M main
git remote add origin https://github.com/USERNAME/populus_pipeline.git
git push -u origin main
```

## 5. Finish the metadata

- In `README.md` and `CITATION.cff`, replace `USERNAME` with your GitHub handle.
- Update `CITATION.cff` and the `LICENSE` copyright line with the final author list.
- On acceptance, add the paper's DOI to `CITATION.cff`.

## Notes

- `.gitignore` keeps `results_grid/`, `results_coldspot/`, caches, and LaTeX aux out
  of version control. The committed `grid_endstate_metrics.csv` is what lets the
  figure rebuild without re-running SLiM.
- If you later regenerate the grid, nothing new gets committed unless you change the
  aggregated CSV — exactly what you want.
- GitHub rejects single files over 100 MB. Nothing here approaches that once the raw
  grid is excluded; the largest tracked file is a figure (well under 1 MB).
