# Gate for sim_factorial_lib.R: at the 2018 design, the new generator must
# produce exactly the old generator text, and make_long() the same long data as
# false_positives_mlm.R's generate_data(); the Cholesky draws must have the implied
# moments; the generating effects must come back out of the data in the
# direction their names say; and the null flags and grid must be as designed.
#
#   Rscript tests/test_factorial_generator.R   (from the repo root)

suppressMessages({ library(dplyr); library(tidyr) })
source("sim_factorial_lib.R")
ex <- parse("false_positives_mlm.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% c("make_clpm_lavaan", "make_riclpm_lavaan", "generate_data")) eval(e, globalenv())
}
fails <- 0L
check <- function(ok, what) {
  cat(if (ok) "ok  " else "FAIL", what, "\n")
  if (!ok) fails <<- fails + 1L
}

base <- list(waves = 12, sx = .7, sy = .7, clxy = 0, clyx = 0, co = .5, tcor = 0)
old_fix <- paste(c(paste0("x", 1:12, " ~~ 0*x", 1:12), paste0("y", 1:12, " ~~ 0*y", 1:12)), collapse = "\n")
check(identical(make_generator_text(c(base, v = 0)), make_clpm_lavaan()), "CLPM text identical at the 2018 design")
for (v in c(.5, 1)) {
  old <- paste(make_riclpm_lavaan(generating_params = list(s = .7, cl = 0, r = 1, v = v, cv = 0, co = .5)),
               old_fix, sep = "\n")
  check(identical(make_generator_text(c(base, v = v)), old), sprintf("RI-CLPM text identical at v = %s", v))
}

for (v in c(0, 1)) {
  gen <- make_generator_text(c(base, v = v))
  set.seed(11); old <- as.data.frame(generate_data(n = 40, gen_mod = gen, wcen = TRUE))
  set.seed(11); new <- make_long(lavaan::simulateData(model = gen, sample.nobs = 40, empirical = FALSE), 12)
  cols <- c("id", "wave", "x", "y", "x_lag", "y_lag", "wcen_x", "gcen_x", "wcen_x_lag", "wcen_y", "gcen_y", "wcen_y_lag")
  check(nrow(old) == nrow(new) && isTRUE(all.equal(old[cols], new[cols], check.attributes = FALSE, tolerance = 1e-12)),
        sprintf("make_long() matches generate_data() at v = %s", v))
}

# The parameters reach the text: a nonzero x -> y effect shows up only on y ~ x,
# and the trait covariance is tcor * v.
g <- make_generator_text(list(waves = 6, sx = .7, sy = .3, clxy = .1, clyx = 0, co = 0, v = .5, tcor = .3))
check(grepl("q6 ~ 0.3*q5 + 0.1*p5", g, fixed = TRUE) && grepl("p6 ~ 0.7*p5 + 0*q5", g, fixed = TRUE), "6-wave stabilities and cross-lags")
check(grepl("kappa ~~ 0.15*omega", g, fixed = TRUE) && !grepl("x7", g, fixed = TRUE), "trait covariance and wave count")
check(inherits(try(make_generator_text(list(waves = 6, sx = .7, sy = .7, clxy = 0, clyx = 0, co = 0, v = 0, tcor = .3)), silent = TRUE), "try-error"),
      "trait correlation without trait variance is an error")

# The Cholesky draws have the implied moments, mean zero and covariance Sigma
# (n = 200,000; sampling SD of a covariance near 1 is about .003).
for (pp in list(c(base, v = 0), list(waves = 6, sx = .7, sy = .3, clxy = .1, clyx = 0, co = .5, v = 1, tcor = .3))) {
  mom <- implied_moments(make_generator_text(pp), pp$waves)
  set.seed(3); w <- simulate_wide(2e5, mom)
  check(max(abs(cov(w) - mom$Sigma)) < .03 && max(abs(colMeans(w))) < .02,
        sprintf("draws reproduce the implied moments (waves = %d, v = %s)", pp$waves, pp$v))
}

# Effects come out where their names put them. CLPM: regressing y[t] on y[t-1]
# and x[t-1] recovers sy and clxy, and a y -> x effect does not leak into it.
recover <- function(pp, n = 1e5) {
  set.seed(5); d <- make_long(simulate_wide(n, implied_moments(make_generator_text(pp), pp$waves)), pp$waves)
  coef(lm(y ~ y_lag + x_lag, data = d))[c("y_lag", "x_lag")]
}
b <- recover(list(waves = 6, sx = .7, sy = .3, clxy = .1, clyx = 0, co = 0, v = 0, tcor = 0))
check(abs(b[["y_lag"]] - .3) < .01 && abs(b[["x_lag"]] - .1) < .01, "CLPM: y[t] ~ y[t-1] + x[t-1] recovers sy = .3, clxy = .1")
b <- recover(list(waves = 6, sx = .7, sy = .3, clxy = 0, clyx = .1, co = 0, v = 0, tcor = 0))
check(abs(b[["x_lag"]]) < .01, "CLPM: a y -> x effect does not show up as x -> y")
# RI-CLPM with unequal stabilities: the riclpm spec's y_x is the x -> y effect.
pp <- list(waves = 6, sx = .7, sy = .3, clxy = .1, clyx = 0, co = .5, v = .5, tcor = .3)
set.seed(7); w <- simulate_wide(2e4, implied_moments(make_generator_text(pp), 6))
r <- specs$riclpm$fit(NULL, w, 6)
check(abs(r$est - .1) < .02, sprintf("RI-CLPM fit recovers clxy = .1 with sx != sy (got %.3f)", r$est))

# The parsed-once RI-CLPM fit gives exactly what riclpmr::lavriclpm() gives.
same <- vapply(1:3, function(i) {
  set.seed(20 + i); w <- simulate_wide(300, implied_moments(make_generator_text(pp), 6))
  a <- parameterEstimates(lavriclpm(riclpmModel = riclpm_text_for(6), data = w))
  a <- a[!is.na(a$label) & a$label == "y_x", ][1, ]
  b <- specs$riclpm$fit(NULL, w, 6)
  identical(c(a$est, a$se, a$pvalue), c(b$est, b$se, b$p))
}, logical(1))
check(all(same), "RI-CLPM fit identical to lavriclpm() on 3 datasets")

# Null flags.
cnd <- function(...) modifyList(list(clxy = 0, clyx = 0, co = 0, tcor = 0), list(...))
check(is_null("lag", cnd(co = .5, tcor = .3)) && !is_null("lag", cnd(clxy = .1)), "lag null depends only on clxy")
check(is_null("within", cnd(tcor = .3)) && !is_null("within", cnd(co = .5)) && !is_null("within", cnd(clxy = .1)),
      "within null: no within-wave cor and no cross-lag; trait cor allowed")
check(is_null("between", cnd()) && !is_null("between", cnd(tcor = .3)) && !is_null("between", cnd(co = .5)),
      "between null: nothing correlated anywhere")
check(all(vapply(specs, function(s) all(s$roles %in% c("lag", "within", "between")), logical(1))), "every spec term has a known role")

# A task's results carry every grid column, one row per spec x term x rep.
g3 <- make_grid(2)
tsk <- g3[g3$n == 50 & g3$waves == 6 & g3$v == .5 & g3$tcor == .3 & g3$clxy == .1 & g3$sy == .3 & g3$co == .5, ]
tsk$extra_factor <- 42
res <- run_task(tsk)
check(all(names(tsk) %in% names(res)) && all(res$extra_factor == 42), "results carry every grid column")
n_terms <- sum(vapply(specs[specs_for(tsk$v)], function(s) length(s$roles), integer(1)))
check(nrow(res) == 2 * n_terms && all(is.na(res$error)), sprintf("task gives %d rows (2 reps x %d terms) and no errors", 2 * n_terms, n_terms))
check(identical(res$null[res$role == "lag"], rep(FALSE, sum(res$role == "lag"))) &&
      all(!res$null[res$role %in% c("within", "between")]), "null flags in a non-null condition are all FALSE")

grid <- make_grid(20)
check(nrow(grid) == 540 && !anyDuplicated(grid$seed), "grid has 540 condition x N tasks with distinct seeds")
g10k <- make_grid(10000)
check(identical(g10k$seed[g10k$chunk == 1], grid$seed), "pilot seeds are the full run's chunk-1 seeds")

if (fails > 0) stop(fails, " check(s) failed")
cat("all checks passed\n")
