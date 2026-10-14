# Summaries of the factorial simulations, one row per cell (condition x N) x
# spec x term role. Each task file is reduced to exact sums (count, sum of
# estimates, sum of squared estimates, rejections, ...) and the sums are added
# over a cell's chunks, so the full results never have to sit in memory at
# once and the statistics are the same as computing them on the pooled reps.
#
# Rates are reported two ways: among fits that returned an estimate
# (`reject_rate`) and over all replications, counting a failed fit as a
# non-rejection (`reject_rate_all`). Failures are counted, never dropped.

alpha <- .05

# Per-task sums, keyed by spec x role (a failed fit has no term name, but its
# role is known). The rows a task must have are checked against the grid row.
task_sums <- function(d, g_row, roles_per_spec) {
  if (!identical(unique(d$seed), g_row$seed)) stop("task ", g_row$task, ": seed does not match the grid")
  want_specs <- specs_for(g_row$v)
  if (!setequal(unique(d$spec), want_specs)) {
    stop("task ", g_row$task, ": specs ", paste(sort(unique(d$spec)), collapse = " "), " but the grid row needs ",
         paste(sort(want_specs), collapse = " "))
  }
  if (anyDuplicated(d[c("rep", "spec", "role")])) stop("task ", g_row$task, ": duplicated rep x spec x role rows")
  expected_rows <- g_row$reps * sum(roles_per_spec[want_specs])
  if (nrow(d) != expected_rows || length(unique(d$rep)) != g_row$reps) {
    stop("task ", g_row$task, ": ", nrow(d), " rows / ", length(unique(d$rep)), " reps, expected ",
         expected_rows, " rows / ", g_row$reps, " reps")
  }
  ok <- is.na(d$error)
  key <- paste(d$spec, d$role, sep = "|")
  s <- function(x) tapply(x, key, sum)
  out <- data.frame(
    key = names(s(rep(1, nrow(d)))),
    n_reps = as.vector(s(rep(1, nrow(d)))),
    n_failed = as.vector(s(!ok)),
    n_warn = as.vector(s(d$n_warnings > 0)),
    n_msg = as.vector(s(d$n_messages > 0)),
    n_valid = as.vector(s(ok)),
    n_reject = as.vector(s(ok & !is.na(d$p) & d$p < alpha)),
    sum_est = as.vector(s(ifelse(ok, d$est, 0))),
    sum_est2 = as.vector(s(ifelse(ok, d$est^2, 0))),
    sum_se = as.vector(s(ifelse(ok, d$se, 0))),
    row.names = NULL)
  terms <- tapply(d$term[ok], key[ok], function(x) paste(unique(x), collapse = "/"))
  out$term <- unname(terms[out$key])
  out$null <- as.vector(tapply(d$null, key, function(x) {
    if (length(unique(x)) != 1) stop("task ", g_row$task, ": mixed null flags within a spec x role")
    x[1]
  })[out$key])
  out$spec <- sub("\\|.*", "", out$key)
  out$role <- sub(".*\\|", "", out$key)
  out$key <- NULL
  out$task <- g_row$task
  out$spec_secs <- as.vector(tapply(d$spec_secs, key, `[`, 1)[paste(out$spec, out$role, sep = "|")])
  out
}

# Add the chunks of each cell together and turn sums into statistics.
cell_summary <- function(sums, grid, reps_per_cell) {
  cond_cols <- c("cond", "n", "waves", "sx", "sy", "clxy", "clyx", "co", "v", "tcor")
  s <- merge(sums, grid[c("task", cond_cols)], by = "task")
  if (nrow(s) != nrow(sums)) stop("task sums do not all match a grid row")
  by <- c(cond_cols, "spec", "role")
  num <- c("n_reps", "n_failed", "n_warn", "n_msg", "n_valid", "n_reject", "sum_est", "sum_est2", "sum_se", "spec_secs")
  a <- aggregate(s[num], s[by], sum)
  first <- function(col) aggregate(s[col], s[by], function(x) {
    u <- unique(x[!is.na(x)])
    if (length(u) > 1) paste(u, collapse = "/") else if (length(u) == 0) NA else u
  })
  a <- merge(a, first("term"), by = by)
  a <- merge(a, first("null"), by = by)
  if (any(a$n_reps != reps_per_cell)) {
    bad <- a[a$n_reps != reps_per_cell, c(by, "n_reps")]
    stop(nrow(bad), " cell x spec rows do not have ", reps_per_cell, " reps (a chunk is missing); first: ",
         paste(capture.output(print(head(bad, 3))), collapse = "\n"))
  }
  nv <- a$n_valid
  a$reject_rate <- a$n_reject / nv
  a$reject_rate_all <- a$n_reject / a$n_reps
  a$mc_se <- sqrt(a$reject_rate * (1 - a$reject_rate) / nv)
  a$fail_rate <- a$n_failed / a$n_reps
  a$warn_rate <- a$n_warn / a$n_reps
  a$msg_rate <- a$n_msg / a$n_reps
  a$mean_est <- a$sum_est / nv
  ss <- a$sum_est2 - a$sum_est^2 / nv
  # Cancellation can make a true zero slightly negative; anything more is a bug.
  if (any(ss < -1e-8 * pmax(a$sum_est2, 1), na.rm = TRUE)) stop("negative sum of squares in sd_est")
  a$sd_est <- sqrt(pmax(ss, 0) / (nv - 1))
  # Ratio of the empirical SD of the estimates to the mean of the reported SEs
  # (not their root mean square); 1 means the SEs are right on average.
  a$mean_se <- a$sum_se / nv
  a$sd_over_se <- a$sd_est / a$mean_se
  a$secs_per_rep <- a$spec_secs / a$n_reps
  a[order(a$cond, a$n, a$spec, a$role), setdiff(names(a), c("sum_est", "sum_est2", "sum_se", "spec_secs"))]
}
