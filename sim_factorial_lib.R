# Functions for the factorial round of the lagged-MLM simulations: a data
# generator with the number of waves, the x and y stabilities, the x -> y
# cross-lag, the within-wave x-y correlation, the trait variance and the trait
# x-y correlation as parameters, and the model specifications, each returning
# estimate, standard error, df and p for the term(s) it is about.
#
# At the 2018 design (12 waves, both stabilities .7, cross-lags 0, within-wave
# covariance .5, trait covariance 0) the generator text is identical to
# make_clpm_lavaan() / make_riclpm_lavaan() + the fixed residual variances used
# by sim_riclpm.R; tests/test_factorial_generator.R checks that, and that
# make_long() gives the same long data as ar_stuff.R's generate_data().
#
# Data are drawn as Z %*% chol(Sigma) + mu from the moments simulateData()
# itself would use (taken from it once per task), not with simulateData() per
# replication: MASS::mvrnorm's eigendecomposition is not reproducible across
# CPUs and BLAS builds, so the same seed gave different data on Modal and
# locally; a Cholesky factor of the same Sigma is.
#
# Every specification in a condition is fit to the same simulated datasets
# (one seed per condition x N x chunk), so spec-to-spec differences are not
# Monte Carlo noise.
#
# Everything with one row per person-wave is fit with lme4 (via lmerTest, for
# Satterthwaite df), except the two AR(1) specifications: lme4 has no residual
# correlation structures, so those stay on nlme with its default optimizer.

suppressMessages({
  library(lavaan)
  library(nlme)
  library(lme4)
  library(lmerTest)
  library(riclpmr)
})

num <- function(v) paste0(format(v, digits = 15, scientific = FALSE, drop0trailing = TRUE), "*")

# p: the generating parameters, a list with waves, sx, sy, clxy (x[t-1] -> y[t]),
# clyx (y[t-1] -> x[t]), co (within-wave residual covariance), v (trait
# variance; 0 gives the CLPM), tcor (trait correlation; covariance tcor * v).
make_generator_text <- function(p) {
  W <- p$waves
  if (p$v == 0) {
    if (p$tcor != 0) stop("trait correlation is undefined when the trait variance is 0")
    xs <- paste0("x", 1:W); ys <- paste0("y", 1:W)
    return(paste(
      paste0(xs[-1], " ~ ", num(p$sx), xs[-W], collapse = "\n"),
      paste0(ys[-1], " ~ ", num(p$sy), ys[-W], collapse = "\n"),
      paste0(c(xs, ys), " ~~ ", num(1), c(xs, ys), collapse = "\n"),
      paste0(xs, " ~~ ", num(p$co), ys, collapse = "\n"),
      paste0(xs[-1], " ~ ", num(p$clyx), ys[-W], collapse = "\n"),
      paste0(ys[-1], " ~ ", num(p$clxy), xs[-W], collapse = "\n"),
      sep = "\n"))
  }
  latent_means <- paste(paste0(c("kappa =~ ", "omega =~ "),
                               c(paste(paste0("1*x", 1:W), collapse = " + "),
                                 paste(paste0("1*y", 1:W), collapse = " + "))),
                        collapse = "\n")
  intercepts <- paste(c(paste0("x", 1:W, " ~ mu", 1:W, "*1"),
                        paste0("y", 1:W, " ~ pi", 1:W, "*1")), collapse = "\n")
  latent_covar <- paste0("\nkappa ~~ ", num(p$v), "kappa #variance\nomega ~~ ", num(p$v),
                         "omega #variance\nkappa ~~ ", num(p$tcor * p$v), "omega #covariance")
  latent_resids <- paste(c(paste0("p", 1:W, " =~ 1*x", 1:W),
                           paste0("q", 1:W, " =~ 1*y", 1:W)), collapse = "\n")
  regressions <- paste(c(paste0("p", W:2, " ~ ", num(p$sx), "p", (W - 1):1, " + ", num(p$clyx), "q", (W - 1):1),
                         paste0("q", W:2, " ~ ", num(p$sy), "q", (W - 1):1, " + ", num(p$clxy), "p", (W - 1):1)),
                       collapse = "\n")
  first_wave <- paste0("\np1 ~~ p1\nq1 ~~ q1\np1 ~~ ", num(p$co), "q1")
  residvar <- paste(c(paste0("p", 2:W, " ~~ ", num(1), "p", 2:W),
                      paste0("q", 2:W, " ~~ ", num(1), "q", 2:W)), collapse = "\n")
  contemporaneous <- paste0("p", 2:W, " ~~ ", num(p$co), "q", 2:W, collapse = "\n")
  # make_riclpm_lavaan() leaves the observed residual variances free, and
  # simulateData() then gives each a variance of 1; fix them at 0 so the data
  # are exactly an RI-CLPM (as in sim_riclpm.R).
  fix_ov_resid <- paste(c(paste0("x", 1:W, " ~~ 0*x", 1:W), paste0("y", 1:W, " ~~ 0*y", 1:W)),
                        collapse = "\n")
  paste(latent_means, intercepts, latent_covar, latent_resids, regressions,
        first_wave, residvar, contemporaneous, fix_ov_resid, sep = "\n")
}

# The model-implied mean and covariance of x1..xW, y1..yW, exactly as
# simulateData() computes them, and the Cholesky factor used for every draw.
implied_moments <- function(gen_text, W) {
  fit <- attr(lavaan::simulateData(model = gen_text, sample.nobs = 10, return.fit = TRUE), "fit")
  imp <- lavInspect(fit, "implied")
  ov <- c(paste0("x", 1:W), paste0("y", 1:W))
  if (!setequal(rownames(imp$cov), ov)) stop("implied covariance has unexpected variables: ", paste(rownames(imp$cov), collapse = " "))
  # Neither generator gives x or y a nonzero mean (the CLPM has no mean
  # structure; the RI-CLPM's intercepts are free, which simulateData() sets to
  # 0), so the draws are mean zero. Anything else is a generator change that
  # simulate_wide() would have to handle.
  if (!is.null(imp$mean) && any(imp$mean != 0)) stop("generator implies nonzero means: ", paste(imp$mean, collapse = " "))
  Sigma <- imp$cov[ov, ov]
  list(Sigma = Sigma, R = chol(Sigma), names = ov)
}

simulate_wide <- function(n, mom) {
  z <- matrix(rnorm(n * length(mom$names)), n) %*% mom$R
  colnames(z) <- mom$names
  as.data.frame(z)
}

# Long data, one row per person-wave, sorted by id then wave, with the same
# columns (and values) as ar_stuff.R's generate_data(wcen = TRUE), plus time_c.
make_long <- function(wide, W) {
  x <- as.matrix(wide[paste0("x", 1:W)])
  y <- as.matrix(wide[paste0("y", 1:W)])
  n <- nrow(x)
  lagm <- function(m) cbind(NA_real_, m[, -W, drop = FALSE])
  long <- function(m) as.vector(t(m))
  mx <- rowMeans(x); my <- rowMeans(y)
  wx <- x - mx; wy <- y - my
  data.frame(id = rep(seq_len(n), each = W), wave = rep(seq_len(W), times = n),
             x = long(x), y = long(y), x_lag = long(lagm(x)), y_lag = long(lagm(y)),
             wcen_x = long(wx), gcen_x = rep(mx - mean(x), each = W), wcen_x_lag = long(lagm(wx)),
             wcen_y = long(wy), gcen_y = rep(my - mean(y), each = W), wcen_y_lag = long(lagm(wy)),
             time_c = rep(seq_len(W) - (W + 1) / 2, times = n))
}

lmer_terms <- function(fit, terms) {
  cs <- coef(summary(fit, ddf = "Satterthwaite"))
  data.frame(term = terms, est = cs[terms, "Estimate"], se = cs[terms, "Std. Error"],
             df = cs[terms, "df"], p = cs[terms, "Pr(>|t|)"], row.names = NULL)
}
lme_terms <- function(fit, terms) {
  tt <- summary(fit)$tTable
  data.frame(term = terms, est = tt[terms, "Value"], se = tt[terms, "Std.Error"],
             df = tt[terms, "DF"], p = tt[terms, "p-value"], row.names = NULL)
}

# The RI-CLPM fit is riclpmr::lavriclpm(), with the model text parsed once per
# number of waves instead of on every fit (12% of the fit's time). The lavaan
# options are lavriclpm()'s own; tests/test_factorial_generator.R checks that
# this gives results identical to calling lavriclpm().
riclpm_lavaan_opts <- list(int.ov.free = FALSE, int.lv.free = FALSE, auto.fix.first = FALSE,
                           auto.fix.single = FALSE, auto.cov.lv.x = FALSE, auto.cov.y = FALSE,
                           auto.var = FALSE)
riclpm_text_for <- function(W) riclpm_text(list(x = paste0("x", 1:W), y = paste0("y", 1:W)), constrain_over_waves = TRUE)
riclpm_partables <- new.env()
fit_riclpm <- function(wide, W) {
  key <- as.character(W)
  if (is.null(riclpm_partables[[key]])) {
    riclpm_partables[[key]] <- do.call(lavaanify, c(list(model = riclpm_text_for(W)), riclpm_lavaan_opts))
  }
  fit <- do.call(lavaan, c(list(model = riclpm_partables[[key]], data = wide), riclpm_lavaan_opts))
  if (!lavInspect(fit, "converged")) stop("RI-CLPM did not converge")
  pe <- parameterEstimates(fit)
  row <- pe[!is.na(pe$label) & pe$label == "y_x", ]
  if (nrow(row) == 0) stop("RI-CLPM fit has no y_x parameter")
  row <- row[1, ]
  data.frame(term = "y_x", est = row$est, se = row$se, df = NA_real_, p = row$pvalue)
}

lmer_spec <- function(f, terms) function(d, w, W) lmer_terms(lmer(f, data = d, REML = TRUE), terms)
lme_ar_spec <- function(f, terms) function(d, w, W) {
  lme_terms(lme(f, random = ~ 1 | id, data = d, correlation = corAR1(form = ~ wave | id), na.action = na.omit), terms)
}

# Each spec takes the long data, the wide data and the number of waves, and
# returns one row per reported term. `roles` says what each term estimates,
# which decides the cells where its true value is zero (see is_null()):
#   lag     - the lagged x -> y effect
#   within  - the concurrent within-person x-y association
#   between - the between-person (person-mean) x-y association
specs <- list(
  naive            = list(fit = lmer_spec(y ~ 1 + x_lag + (1 | id), "x_lag"), roles = "lag"),
  sea              = list(fit = lmer_spec(y ~ 1 + wcen_x_lag + gcen_x + time_c + (1 | id), "wcen_x_lag"), roles = "lag"),
  wcen_x_raw_ylag  = list(fit = lmer_spec(y ~ 1 + y_lag + wcen_x_lag + gcen_x + (1 | id), "wcen_x_lag"), roles = "lag"),
  no_wcen_ylag     = list(fit = lmer_spec(y ~ 1 + y_lag + x_lag + (1 | id), "x_lag"), roles = "lag"),
  wcen_ylag        = list(fit = lmer_spec(y ~ 1 + wcen_y_lag + wcen_x_lag + gcen_x + (1 | id), "wcen_x_lag"), roles = "lag"),
  wcen_ylag_wcenre = list(fit = lmer_spec(y ~ 1 + wcen_y_lag + wcen_x_lag + gcen_x + (1 + wcen_x_lag || id), "wcen_x_lag"),
                          roles = "lag"),
  no_wcen_ar       = list(fit = lme_ar_spec(y ~ 1 + x_lag, "x_lag"), roles = "lag"),
  wcen_ar          = list(fit = lme_ar_spec(y ~ 1 + wcen_x_lag + gcen_x, "wcen_x_lag"), roles = "lag"),
  riclpm           = list(fit = function(d, w, W) fit_riclpm(w, W), roles = "lag"),
  # The concurrent within/between model, and the same with the prior outcome added
  # (raw, then person-mean centered): does x go with y "over and above" y[t-1]?
  sea_concurrent           = list(fit = lmer_spec(y ~ 1 + wcen_x + gcen_x + time_c + (1 | id), c("wcen_x", "gcen_x")),
                                  roles = c("within", "between")),
  sea_concurrent_ylag      = list(fit = lmer_spec(y ~ 1 + wcen_x + gcen_x + y_lag + time_c + (1 | id), c("wcen_x", "gcen_x")),
                                  roles = c("within", "between")),
  sea_concurrent_wcen_ylag = list(fit = lmer_spec(y ~ 1 + wcen_x + gcen_x + wcen_y_lag + time_c + (1 | id), c("wcen_x", "gcen_x")),
                                  roles = c("within", "between"))
)

# Whether a term's true value is zero in a condition. The lagged effect is
# zero when x[t-1] -> y[t] is. A concurrent within-person association also
# arises from a within-wave correlation or from either cross-lag (x[t-1] -> y[t]
# with a stable x makes x[t] and y[t] covary). A between-person association
# arises from the trait correlation, and, through person means of a finite
# number of waves, from everything that makes the within-person parts covary.
is_null <- function(role, cond) {
  switch(role,
         lag     = cond$clxy == 0,
         within  = cond$co == 0 && cond$clxy == 0 && cond$clyx == 0,
         between = cond$tcor == 0 && cond$co == 0 && cond$clxy == 0 && cond$clyx == 0,
         stop("unknown term role: ", role))
}

# Specs fit in a condition: all, except the RI-CLPM when there is no trait
# variance (its random-intercept variance would sit at the boundary of zero).
specs_for <- function(v) if (v == 0) setdiff(names(specs), "riclpm") else names(specs)

# The grid, one row per (condition, N, chunk); every spec is run on each
# task's datasets. Seeds depend only on the row's position in the
# deterministic condition grid and the chunk number, never on reps per cell,
# so a pilot's chunk 1 is the start of the full run's chunk 1.
seed_base <- 60261013L
chunk_size_max <- 250L
make_grid <- function(reps_per_cell) {
  cond <- expand.grid(n = c(50, 150, 600), sy = c(0, .3, .7), co = c(0, .5), waves = c(6, 12, 24),
                      clxy = c(0, .1), v = c(0, .5, 1), tcor = c(0, .3), stringsAsFactors = FALSE)
  cond$sx <- .7
  cond$clyx <- 0
  cond$cond <- seq_len(nrow(cond))
  # The trait correlation is undefined without trait variance.
  cond <- cond[!(cond$v == 0 & cond$tcor > 0), ]
  chunk_size <- min(chunk_size_max, reps_per_cell)
  n_chunks <- ceiling(reps_per_cell / chunk_size)
  if (n_chunks >= 1000) stop("more than 999 chunks per cell would collide seeds")
  g <- merge(cond, data.frame(chunk = seq_len(n_chunks)), by = NULL)
  g <- g[order(g$cond, g$chunk), ]
  g$reps <- ifelse(g$chunk < n_chunks, chunk_size, reps_per_cell - chunk_size * (n_chunks - 1))
  g$seed <- seed_base + g$cond * 1000L + g$chunk
  if (anyDuplicated(g$seed)) stop("seed collision between chunks")
  g$task <- seq_len(nrow(g))
  rownames(g) <- NULL
  g
}

# Errors become rows with NA estimates and the message in `error`; warnings and
# messages (e.g. lme4's singular-fit message) are counted and the estimate is
# kept, so the rates are not conditioned on clean fits.
fit_spec <- function(spec_name, d, w, W) {
  n_warn <- 0L; n_msg <- 0L
  sp <- specs[[spec_name]]
  na_rows <- function(msg) data.frame(term = NA_character_, role = sp$roles, est = NA_real_, se = NA_real_,
                                      df = NA_real_, p = NA_real_, error = msg)
  r <- tryCatch({
    r <- withCallingHandlers(sp$fit(d, w, W),
                             warning = function(w) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") },
                             message = function(m) { n_msg <<- n_msg + 1L; invokeRestart("muffleMessage") })
    if (!all(is.finite(r$est)) || !all(is.finite(r$se)) || !all(is.finite(r$p))) {
      stop("fit returned a non-finite estimate, SE or p-value")
    }
    # The RI-CLPM is a z test and has no df; every other spec must return a finite one.
    if (!all(is.finite(r$df) | (spec_name == "riclpm" & is.na(r$df)))) stop("fit returned a non-finite df")
    r$role <- sp$roles
    r$error <- NA_character_
    r
  }, error = function(e) na_rows(conditionMessage(e)))
  r$spec <- spec_name
  r$n_warnings <- n_warn; r$n_messages <- n_msg
  r
}

run_task <- function(t) {
  set.seed(t$seed)
  W <- t$waves
  mom <- implied_moments(make_generator_text(as.list(t[c("waves", "sx", "sy", "clxy", "clyx", "co", "v", "tcor")])), W)
  sp_names <- specs_for(t$v)
  spec_secs <- setNames(numeric(length(sp_names)), sp_names)
  t0 <- proc.time()[["elapsed"]]
  res <- vector("list", t$reps)
  for (i in seq_len(t$reps)) {
    w <- simulate_wide(t$n, mom)
    d <- make_long(w, W)
    res[[i]] <- do.call(rbind, lapply(sp_names, function(s) {
      s0 <- proc.time()[["elapsed"]]
      r <- fit_spec(s, d, w, W)
      spec_secs[[s]] <<- spec_secs[[s]] + proc.time()[["elapsed"]] - s0
      r
    }))
    res[[i]]$rep <- i
  }
  res <- do.call(rbind, res)
  res$null <- vapply(res$role, is_null, logical(1), cond = t)
  res$spec_secs <- spec_secs[res$spec]
  # Every grid column is carried into the results, so a factor added to the
  # grid cannot be left out of them.
  if (any(names(t) %in% names(res))) stop("grid and result columns collide: ", paste(intersect(names(t), names(res)), collapse = " "))
  cbind(t[rep(1, nrow(res)), ], res, chunk_secs = proc.time()[["elapsed"]] - t0, row.names = NULL)
}
