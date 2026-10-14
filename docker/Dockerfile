# The one R environment for lagged_mlm: simulations (sim_factorial.R,
# run_local.py, and modal_sim.py, which builds its image from this file),
# tests, summaries, and rendering false_positives_mlm.R.
#
#   docker build -t lagged-mlm docker/
#
# R 4.3.2 with CRAN packages from the Posit Package Manager snapshot of
# 2024-04-01 and riclpmr at a fixed commit. Changing the snapshot or the
# riclpmr commit changes the run manifest (modal_sim.manifest), so existing
# result directories will refuse to take new tasks; that is intended.
FROM rocker/r-ver:4.3.2

ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

RUN /rocker_scripts/install_pandoc.sh

RUN apt-get update && apt-get install -y --no-install-recommends \
      libpng16-16 libjpeg-turbo8 libxml2 libglpk40 libcairo2 \
    && rm -rf /var/lib/apt/lists/*

ARG PPM=https://packagemanager.posit.co/cran/__linux__/jammy/2024-04-01
ARG RICLPMR_SHA=76822ab8f1d4ddf27a810c359489df69256abd86
ARG PKGS="Matrix lme4 lmerTest lavaan nlme dplyr tidyr simsem semPlot ggplot2 rmarkdown knitr jsonlite tidySEM"

RUN Rscript -e "options(repos = c(CRAN = '${PPM}')); install.packages(strsplit('${PKGS}', ' ')[[1]])"
RUN Rscript -e "install.packages('https://github.com/jflournoy/riclpmr/archive/${RICLPMR_SHA}.tar.gz', repos = NULL, type = 'source')"
# install.packages() only warns when a package fails; loading each one here
# makes the build fail instead.
RUN Rscript -e "for (p in c(strsplit('${PKGS}', ' ')[[1]], 'riclpmr')) { suppressMessages(library(p, character.only = TRUE)); cat(p, as.character(packageVersion(p)), '\n') }" \
    && Rscript -e "stopifnot(rmarkdown::pandoc_available())"
