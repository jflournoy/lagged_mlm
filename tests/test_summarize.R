# The summary's sums-then-statistics route must equal computing the statistics
# on the pooled replications, failures must be counted and never dropped, and
# a cell with a missing chunk must stop the summary.
#
#   Rscript tests/test_summarize.R   (from the repo root)

suppressMessages(source("sim_factorial_lib.R"))
source("summarize_factorial_lib.R")
fails <- 0L
check <- function(ok, what) { cat(if (ok) "ok  " else "FAIL", what, "\n"); if (!ok) fails <<- fails + 1L }

grid <- make_grid(6)
g <- grid[grid$v == 0, ][1:2, ]
g$reps <- 3L
grid <- grid[grid$task %in% g$task, ]
grid$reps <- 3L
roles_per_spec <- vapply(specs, function(s) length(s$roles), integer(1))
set.seed(1)
fake_task <- function(gr) {
  sp <- specs_for(gr$v)
  do.call(rbind, lapply(seq_len(gr$reps), function(r) do.call(rbind, lapply(sp, function(s) {
    ro <- specs[[s]]$roles
    data.frame(task = gr$task, seed = gr$seed, spec = s, role = ro, term = paste0("t_", ro), rep = r,
               est = rnorm(length(ro)), se = runif(length(ro), .5, 1.5), p = runif(length(ro)),
               df = 10, error = NA_character_, n_warnings = 0L, n_messages = 0L,
               null = TRUE, spec_secs = 1)
  }))))
}
d1 <- fake_task(g[1, ]); d2 <- fake_task(g[2, ])
# Make one fit fail in task 1 and one carry a warning.
d1$error[1] <- "boom"; d1$est[1] <- NA; d1$se[1] <- NA; d1$p[1] <- NA; d1$term[1] <- NA
d1$n_warnings[which(d1$spec == d1$spec[1] & d1$role == d1$role[1] & d1$rep == 2)] <- 1L
grid$cond[grid$task == g$task[2]] <- grid$cond[grid$task == g$task[1]]
grid$n[grid$task == g$task[2]] <- grid$n[grid$task == g$task[1]]
grid$chunk <- 1:2
sums <- rbind(task_sums(d1, g[1, ], roles_per_spec), task_sums(d2, g[2, ], roles_per_spec))
cells <- cell_summary(sums, grid, 6L)
pooled <- rbind(d1, d2)
x <- pooled[pooled$spec == d1$spec[1] & pooled$role == d1$role[1], ]
row <- cells[cells$spec == d1$spec[1] & cells$role == d1$role[1], ]
ok <- is.na(x$error)
check(nrow(row) == 1 && row$n_reps == 6 && row$n_failed == 1 && row$n_valid == 5, "failed fit counted, not dropped")
check(isTRUE(all.equal(row$mean_est, mean(x$est[ok]))) && isTRUE(all.equal(row$sd_est, sd(x$est[ok]))) &&
      isTRUE(all.equal(row$mean_se, mean(x$se[ok]))), "mean, sd and mean SE equal the pooled computation")
check(isTRUE(all.equal(row$reject_rate, mean(x$p[ok] < .05))) && isTRUE(all.equal(row$reject_rate_all, sum(x$p[ok] < .05) / 6)),
      "rejection rates over valid fits and over all reps")
check(row$warn_rate == 1 / 6, "warning counted")
check(inherits(try(cell_summary(sums[sums$task == g$task[1], ], grid, 6L), silent = TRUE), "try-error"),
      "a cell missing a chunk stops the summary")
bad <- d1[-1, ]
check(inherits(try(task_sums(bad, g[1, ], roles_per_spec), silent = TRUE), "try-error"), "a task with a missing row stops")
bad <- d1; bad$seed <- bad$seed + 1
check(inherits(try(task_sums(bad, g[1, ], roles_per_spec), silent = TRUE), "try-error"), "a task with the wrong seed stops")
bad <- d1[d1$spec != "naive", ]
check(inherits(try(task_sums(bad, g[1, ], roles_per_spec), silent = TRUE), "try-error"), "a task missing a whole spec stops")
bad <- d1; i <- which(bad$spec == "sea" & bad$rep == 1); bad$rep[i] <- 2L
bad <- bad[!(bad$spec == "sea" & bad$rep == 3), ]
bad <- rbind(bad, bad[bad$spec == "naive" & bad$rep == 3, ][0, ])
check(inherits(try(task_sums(bad, g[1, ], roles_per_spec), silent = TRUE), "try-error"), "duplicated rep x spec x role rows stop")
# An all-failed spec in a cell stays in the summary, with NaN statistics.
af1 <- d1; af2 <- d2
for (k in c("af1", "af2")) {
  x <- get(k); j <- x$spec == "riclpm" | x$spec == "naive"
  x$error[j] <- "boom"; x$est[j] <- NA; x$se[j] <- NA; x$p[j] <- NA; x$term[j] <- NA
  assign(k, x)
}
cs <- cell_summary(rbind(task_sums(af1, g[1, ], roles_per_spec), task_sums(af2, g[2, ], roles_per_spec)), grid, 6L)
nv <- cs[cs$spec == "naive", ]
check(nrow(nv) == 1 && nv$n_failed == 6 && nv$fail_rate == 1 && is.nan(nv$reject_rate) && nv$reject_rate_all == 0,
      "an all-failed cell is kept, fail_rate 1, reject_rate NaN")

if (fails > 0) stop(fails, " check(s) failed")
cat("all summary checks passed\n")
