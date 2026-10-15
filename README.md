# lagged_mlm

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23383372.svg)](https://doi.org/10.5281/zenodo.23383372)

Simulations behind "False positives in multilevel models for longitudinal data"
(`false_positives_mlm.R`, published at <https://johnflournoy.science/lagged_mlm/> from `docs/`).

## Environment

Everything runs in one Docker image, built from `docker/Dockerfile`: R 4.3.2, CRAN packages from the
Posit Package Manager snapshot of 2024-04-01, and riclpmr at a fixed commit. `modal_sim.py` builds its
Modal image from the same file.

```bash
docker build -t lagged-mlm docker/
```

## Commands

Tests:

```bash
docker run --rm -v $PWD:/w -w /w lagged-mlm Rscript tests/test_factorial_generator.R
docker run --rm -v $PWD:/w -w /w lagged-mlm Rscript tests/test_summarize.R
python tests/test_modal_manifest.py
```

The factorial simulations, locally (resumable; rerun the same command after an interruption) or on Modal:

```bash
python run_local.py --reps 2500 --n 50,150 --waves 6,12 --cores 16
modal run --detach modal_sim.py --reps 2500
```

Summarize a finished run (writes `summary_cells.csv` / `.RDS` into the run directory):

```bash
docker run --rm -v $PWD:/w -w /w lagged-mlm Rscript summarize_factorial.R factorial_sims_r2500
```

Render the page and publish it:

```bash
docker run --rm --user $(id -u):$(id -g) -e HOME=/tmp -v $PWD:/r -w /r lagged-mlm \
  Rscript -e 'rmarkdown::render("false_positives_mlm.R", knit_root_dir = "/r")'
cp false_positives_mlm.html docs/index.html
```

The page header (`false_positives_mlm_head.html`) carries the site's Simple Analytics script, marked `data-external="1"` so pandoc links it instead of inlining it.

The Python launchers need `modal` (for `modal_sim.py`, and imported by `run_local.py` for the shared
task selection and manifest checks).
