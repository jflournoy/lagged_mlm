# Null simulations with a stable trait component: the same lagged-MLM
# specifications as ar_stuff.R, but with data generated from the RI-CLPM
# (make_riclpm_lavaan) instead of the CLPM, plus the RI-CLPM itself fit with
# riclpmr as the correctly specified benchmark. Also runs the naive model
# (y on lagged x with a random intercept and nothing else), including on the
# CLPM data (trait variance 0) used for the original simulations.
#
# All nlme fits use nlme's default optimizer (nlminb). The 2018 run_model()
# used lmeControl(opt = "optim"), which stops short of the maximum for the
# centered AR(1) model at N around 300; a fit that fails to converge is
# recorded as a failure, never refit with another optimizer.
#
# make_riclpm_lavaan() defines x and y as indicators of the random intercepts
# and within-person deviations but does not fix their residual variances, and
# lavaan::simulateData() then gives each observed variable a residual variance
# of 1 (measurement error the model never mentions). The generator below adds
# those residual variances fixed to 0, so the data are exactly an RI-CLPM.
#
# Usage: Rscript sim_riclpm.R <reps_per_cell> <cores> [out_dir]
# (out_dir defaults to riclpm_sims_r<reps_per_cell>, next to this script)
#
# Each cell (trait variance x specification x N) is split into chunks of
# chunk_size replications. Every chunk sets its own seed and saves its own RDS,
# so the run can be stopped and resumed, and is reproducible regardless of how
# chunks are scheduled across cores. A replication whose model fails to fit is
# recorded with its error message (est and p are NA); failures are counted and
# reported in the summary, never dropped.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript sim_riclpm.R <reps_per_cell> <cores> [out_dir]")
reps_per_cell <- as.integer(args[[1]])
cores <- as.integer(args[[2]])
if (is.na(reps_per_cell) || reps_per_cell < 1) stop("reps_per_cell must be a positive integer")
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1) stop("run this with Rscript so the script's own directory can be found")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
out_dir <- if (length(args) >= 3) args[[3]] else file.path(script_dir, sprintf("riclpm_sims_r%d", reps_per_cell))
if (is.na(cores) || cores < 1) stop("cores must be a positive integer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(simsem)
  library(nlme)
  library(lmerTest)
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(lavaan)
  library(riclpmr)
  library(parallel)
})

# Reuse the exact function definitions from the write-up rather than copies.
ex <- parse(file.path(script_dir, "ar_stuff.R"))
needed <- c("make_clpm_lavaan", "make_riclpm_lavaan", "generate_data", "run_model")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% needed) eval(e, globalenv())
}
missing_fns <- needed[!vapply(needed, exists, logical(1))]
if (length(missing_fns) > 0) stop("not found in ar_stuff.R: ", paste(missing_fns, collapse = ", "))

N_seq <- c(seq(30, 150, 40), seq(300, 1050, 150))
trait_vars <- c(0, 0.5, 1)
chunk_size <- min(500L, reps_per_cell)
n_chunks <- ceiling(reps_per_cell / chunk_size)

specs <- list(
  naive                    = NULL,
  sea                      = NULL,
  wcen_x_raw_ylag          = NULL,
  no_wcen_ar               = list(wcen = FALSE, use_lag_y = FALSE, use_lme4 = FALSE, wcen_re = FALSE),
  no_wcen_ylag             = list(wcen = FALSE, use_lag_y = TRUE,  use_lme4 = FALSE, wcen_re = FALSE),
  wcen_ar                  = list(wcen = TRUE,  use_lag_y = FALSE, use_lme4 = FALSE, wcen_re = FALSE),
  wcen_ylag                = list(wcen = TRUE,  use_lag_y = TRUE,  use_lme4 = FALSE, wcen_re = FALSE),
  wcen_ylag_l4             = list(wcen = TRUE,  use_lag_y = TRUE,  use_lme4 = TRUE,  wcen_re = FALSE),
  wcen_ylag_l4_wcenre      = list(wcen = TRUE,  use_lag_y = TRUE,  use_lme4 = TRUE,  wcen_re = TRUE),
  riclpm                   = NULL
)

riclpm_fit_text <- riclpm_text(list(x = paste0("x", 1:12), y = paste0("y", 1:12)),
                               constrain_over_waves = TRUE)

# The RI-CLPM benchmark: fit to the wide data; y_x is the constrained
# cross-lagged path from x at t-1 to y at t (the analogue of x_lag above).
run_riclpm <- function(n, gen_mod) {
  d <- simulateData(model = gen_mod, sample.nobs = n, empirical = FALSE)
  fit <- lavriclpm(riclpmModel = riclpm_fit_text, data = d)
  if (!lavInspect(fit, "converged")) stop("RI-CLPM did not converge")
  pe <- parameterEstimates(fit)
  row <- pe[!is.na(pe$label) & pe$label == "y_x", ]
  if (nrow(row) == 0) stop("RI-CLPM fit has no y_x parameter")
  row <- row[1, ]
  if (!is.finite(row$est) || !is.finite(row$pvalue)) stop("RI-CLPM fit returned no estimate or p-value for y_x")
  data.frame(n = n, x_lag_p = row$pvalue, x_lag_bias = row$est)
}

# The naive model: no lagged outcome, no residual autocorrelation.
run_naive <- function(n, gen_mod) {
  d <- generate_data(n = n, gen_mod = gen_mod, wcen = FALSE)
  fit <- lme(y ~ 1 + x_lag, random = ~ 1 | id, data = d, na.action = na.omit)
  tt <- summary(fit)$tTable
  data.frame(n = n, x_lag_p = tt["x_lag", "p-value"], x_lag_bias = tt["x_lag", "Value"])
}

# The within/between model used for the month-level lagged analyses in the SEA
# papers: lagged person-mean-centered x, the person mean (grand-mean centered),
# and time, with a random intercept; no prior outcome and no AR structure. (The
# papers fit it in brms with weak normal(0, 10) priors; lme gives essentially the
# same estimates.)
run_sea <- function(n, gen_mod) {
  d <- generate_data(n = n, gen_mod = gen_mod, wcen = TRUE)
  d$time_c <- d$wave - mean(d$wave)
  fit <- lme(y ~ 1 + wcen_x_lag + gcen_x + time_c, random = ~ 1 | id, data = d, na.action = na.omit)
  tt <- summary(fit)$tTable
  data.frame(n = n, x_lag_p = tt["wcen_x_lag", "p-value"], x_lag_bias = tt["wcen_x_lag", "Value"])
}

# Centered x (within and between) with the raw, uncentered lagged outcome and a
# random intercept: what you get if you center x but not the lagged y.
run_wcen_x_raw_ylag <- function(n, gen_mod) {
  d <- generate_data(n = n, gen_mod = gen_mod, wcen = TRUE)
  fit <- lme(y ~ 1 + y_lag + wcen_x_lag + gcen_x, random = ~ 1 | id, data = d, na.action = na.omit)
  tt <- summary(fit)$tTable
  data.frame(n = n, x_lag_p = tt["wcen_x_lag", "p-value"], x_lag_bias = tt["wcen_x_lag", "Value"])
}

fix_ov_resid <- paste(c(paste0("x", 1:12, " ~~ 0*x", 1:12),
                        paste0("y", 1:12, " ~~ 0*y", 1:12)), collapse = "\n")

make_generator <- function(trait_var) {
  if (trait_var == 0) return(make_clpm_lavaan())
  paste(make_riclpm_lavaan(generating_params = list(s = .7, cl = 0, r = 1,
                                                    v = trait_var, cv = 0, co = .5)),
        fix_ov_resid, sep = "\n")
}

# Errors become rows with est and p NA and the message in `error`; warnings
# (e.g., singular fits, lavaan variance warnings) are counted in `n_warnings`
# and the estimate is kept, so the rates are not conditioned on clean fits.
one_rep <- function(n, gen_mod, spec_name) {
  n_warn <- 0L
  tryCatch({
    r <- withCallingHandlers({
      if (spec_name == "riclpm") run_riclpm(n, gen_mod) else
        if (spec_name == "naive") run_naive(n, gen_mod) else
        if (spec_name == "sea") run_sea(n, gen_mod) else
        if (spec_name == "wcen_x_raw_ylag") run_wcen_x_raw_ylag(n, gen_mod) else
          do.call(run_model, c(list(n = n, gen_mod = gen_mod), specs[[spec_name]]))
    }, warning = function(w) {
      n_warn <<- n_warn + 1L
      invokeRestart("muffleWarning")
    })
    if (!is.finite(r$x_lag_bias) || !is.finite(r$x_lag_p)) stop("fit returned a non-finite estimate or p-value")
    r$error <- NA_character_
    r$n_warnings <- n_warn
    r
  }, error = function(e) {
    data.frame(n = n, x_lag_p = NA_real_, x_lag_bias = NA_real_,
               error = conditionMessage(e), n_warnings = n_warn)
  })
}

# One seed per chunk, fixed by the chunk's position in a deterministic
# expand.grid order, before any reordering for scheduling. Specs added after the
# first full run get their own seed range each (`later_specs`, spec = seed base),
# so adding one never changes the seeds of chunks that already exist.
later_specs <- c(sea = 30261012L, wcen_x_raw_ylag = 40261012L)
make_tasks <- function(spec_names, seed_base) {
  tk <- expand.grid(chunk = seq_len(n_chunks), n = N_seq, spec = spec_names,
                    trait_var = trait_vars, stringsAsFactors = FALSE)
  # Trait variance 0 is the CLPM generator used in 2018; every spec is rerun on
  # it here (with nlme's default optimizer) except the RI-CLPM fit, whose
  # random-intercept variance would sit at its boundary of zero.
  tk <- tk[tk$trait_var > 0 | tk$spec != "riclpm", ]
  tk$seed <- seed_base + seq_len(nrow(tk))
  tk
}
tasks <- do.call(rbind, c(list(make_tasks(setdiff(names(specs), names(later_specs)), 20261012L)),
                          lapply(names(later_specs), function(sp) make_tasks(sp, later_specs[[sp]]))))
if (anyDuplicated(tasks$seed)) stop("seed collision between chunks")
tasks$file <- file.path(out_dir, sprintf("v%s_%s_n%04d_c%03d.RDS",
                                         tasks$trait_var, tasks$spec, tasks$n, tasks$chunk))
# Biggest jobs first so the slow cells don't straggle at the end.
tasks <- tasks[order(-tasks$n, tasks$spec == "riclpm"), ]

tasks$reps <- ifelse(tasks$chunk < n_chunks, chunk_size, reps_per_cell - chunk_size * (n_chunks - 1))

# A chunk on disk counts as done only if it reads cleanly and has the expected
# number of rows and the expected seed; anything else stops the run here,
# before any compute, rather than being reused.
check_chunk <- function(f, reps, seed) {
  r <- tryCatch(readRDS(f), error = function(e) stop("unreadable chunk ", f, ": ", conditionMessage(e)))
  if (nrow(r) != reps || !all(r$seed == seed)) {
    stop("chunk ", f, " has ", nrow(r), " rows (expected ", reps, ") or the wrong seed; delete it or use another out_dir")
  }
  invisible(TRUE)
}
done <- file.exists(tasks$file)
invisible(mapply(check_chunk, tasks$file[done], tasks$reps[done], tasks$seed[done]))
todo <- tasks[!done, ]
message(sprintf("%d chunks total, %d to run, %d reps per cell, %d cores",
                nrow(tasks), nrow(todo), reps_per_cell, cores))

run_task <- function(i) {
  t <- todo[i, ]
  set.seed(t$seed)
  gen_mod <- make_generator(t$trait_var)
  t0 <- proc.time()[["elapsed"]]
  res <- bind_rows(replicate(t$reps, one_rep(t$n, gen_mod, t$spec), simplify = FALSE))
  res$spec <- t$spec
  res$trait_var <- t$trait_var
  res$chunk <- t$chunk
  res$seed <- t$seed
  res$chunk_secs <- proc.time()[["elapsed"]] - t0
  tmp <- paste0(t$file, ".tmp")
  saveRDS(res, tmp)
  if (!file.rename(tmp, t$file)) stop("could not move ", tmp, " into place")
  sprintf("%s done", basename(t$file))
}

started <- Sys.time()
out <- mclapply(seq_len(nrow(todo)), run_task, mc.cores = cores, mc.preschedule = FALSE)
crashed <- vapply(out, inherits, logical(1), what = "try-error")
if (any(crashed)) {
  stop(sum(crashed), " chunks crashed:\n",
       paste(vapply(out[crashed], as.character, character(1)), collapse = "\n"))
}
message("elapsed: ", format(Sys.time() - started))
missing_chunks <- tasks$file[!file.exists(tasks$file)]
if (length(missing_chunks) > 0) {
  stop(length(missing_chunks), " chunks were not written (a worker may have been killed):\n",
       paste(head(missing_chunks, 20), collapse = "\n"))
}
invisible(mapply(check_chunk, tasks$file, tasks$reps, tasks$seed))

all_res <- bind_rows(lapply(tasks$file, readRDS))
summary_tab <- all_res %>%
  group_by(trait_var, spec, n) %>%
  summarize(reps = n(),
            failed = sum(!is.na(error)),
            n_valid = sum(!is.na(x_lag_p)),
            reps_with_warnings = sum(n_warnings > 0),
            fp_rate = mean(x_lag_p < .05, na.rm = TRUE),
            mc_se = sqrt(fp_rate * (1 - fp_rate) / n_valid),
            mean_est = mean(x_lag_bias, na.rm = TRUE),
            sd_est = sd(x_lag_bias, na.rm = TRUE),
            secs_per_rep = sum(chunk_secs[!duplicated(chunk)]) / n(),
            .groups = "drop")
if (any(summary_tab$reps != reps_per_cell)) stop("some cells do not have ", reps_per_cell, " replications")
if (any(summary_tab$n_valid + summary_tab$failed != summary_tab$reps)) stop("replications that are neither valid nor recorded as failed")
saveRDS(all_res, file.path(out_dir, "riclpm_sims_all.RDS"))
saveRDS(summary_tab, file.path(out_dir, "riclpm_sims_summary.RDS"))
print(as.data.frame(summary_tab), digits = 3)
if (any(summary_tab$failed > 0)) {
  message("WARNING: some fits failed; see the 'failed' column and the 'error' column in riclpm_sims_all.RDS")
}
if (any(summary_tab$reps_with_warnings > 0)) {
  message("NOTE: some fits raised warnings (estimates kept); see 'reps_with_warnings' and the 'n_warnings' column")
}
