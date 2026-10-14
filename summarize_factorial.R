# Summarize a factorial run directory written by run_local.py or modal_sim.py.
#
#   Rscript summarize_factorial.R <run_dir> [--partial]
#
# Writes <run_dir>/summary_cells.csv and .RDS (one row per cell x spec x term
# role; see summarize_factorial_lib.R) and prints the null-cell false-positive
# rates. The tasks to expect are the union of the launchers' selection records
# (selection_*.json); a selected task with no file stops the summary, as does
# R code that differs from the hashes in manifest.json. --partial summarizes
# only the cells whose tasks are all present, for a look at a run still in
# progress, and marks every output as partial.

options(width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2 || (length(args) == 2 && args[[2]] != "--partial")) {
  stop("usage: Rscript summarize_factorial.R <run_dir> [--partial]")
}
run_dir <- args[[1]]
partial <- length(args) == 2
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
suppressMessages(source(file.path(script_dir, "sim_factorial_lib.R")))
source(file.path(script_dir, "summarize_factorial_lib.R"))

mf <- jsonlite::fromJSON(file.path(run_dir, "manifest.json"), simplifyVector = FALSE)
code_now <- vapply(names(mf$code_sha256), function(f) {
  sub(" .*", "", system2("sha256sum", file.path(script_dir, f), stdout = TRUE))
}, character(1))
if (!identical(unname(code_now), unname(unlist(mf$code_sha256)))) {
  stop("R code differs from the code that wrote these results (manifest.json code_sha256): ",
       paste(names(code_now)[code_now != unlist(mf$code_sha256)], collapse = " "))
}
grid <- read.csv(file.path(run_dir, "grid.csv"))
# read.csv types integers as integer and numbers as double, as make_grid does
# not always, so compare values, not storage modes.
if (!isTRUE(all.equal(make_grid(mf$reps), grid, check.attributes = FALSE))) {
  stop("grid.csv does not match make_grid(", mf$reps, ") from the current code")
}
roles_per_spec <- vapply(specs, function(s) length(s$roles), integer(1))

sel_files <- list.files(run_dir, "^selection_.*\\.json$", full.names = TRUE)
if (length(sel_files) == 0) stop("no selection_*.json in ", run_dir, ": which tasks were launched is unknown")
expected <- sort(unique(unlist(lapply(sel_files, function(f) jsonlite::fromJSON(f)$tasks))))
if (!all(expected %in% grid$task)) stop("selection records name tasks not in the grid")

files <- list.files(run_dir, "^task_\\d{5}\\.RDS$", full.names = TRUE)
ids <- as.integer(sub("^task_(\\d{5})\\.RDS$", "\\1", basename(files)))
if (!all(ids %in% expected)) stop("task files that no selection record names: ", paste(head(setdiff(ids, expected), 20), collapse = " "))
cell_key <- function(g) paste(g$cond, g$n)
missing_tasks <- setdiff(expected, ids)
use <- grid[grid$task %in% expected, ]
if (length(missing_tasks) > 0) {
  if (!partial) stop(length(missing_tasks), " selected tasks have no file (run still going, or tasks lost); rerun with --partial to look anyway")
  bad_cells <- unique(cell_key(grid[grid$task %in% missing_tasks, ]))
  use <- use[!cell_key(use) %in% bad_cells, ]
}
banner <- if (partial) sprintf("*** PARTIAL: %d of %d selected tasks done; %d cells left out; NOT a final summary ***",
                               length(ids), length(expected), length(unique(cell_key(grid[grid$task %in% missing_tasks, ])))) else ""
if (partial) message(banner)
if (nrow(use) == 0) stop("no complete cells to summarize")
message(length(unique(cell_key(use))), " complete cells from ", nrow(use), " task files")

sums <- do.call(rbind, lapply(seq_len(nrow(use)), function(i) {
  g <- use[i, ]
  task_sums(readRDS(file.path(run_dir, sprintf("task_%05d.RDS", g$task))), g, roles_per_spec)
}))
cells <- cell_summary(sums, grid, mf$reps)
cells$partial <- partial
out_base <- file.path(run_dir, if (partial) "summary_cells_PARTIAL" else "summary_cells")
saveRDS(cells, paste0(out_base, ".RDS"))
write.csv(cells, paste0(out_base, ".csv"), row.names = FALSE)
message("wrote ", out_base, ".csv/.RDS (", nrow(cells), " rows)")

# Pooled view of the null cells: the plain mean of cell rates over whatever
# cells are in each row (not a design-weighted marginal). Rows pool different
# cells: the RI-CLPM is only fit with trait variance, and each role has its own
# null set, so `cells` says how many went in. A NaN anywhere stops the table
# rather than dropping the cell.
pooled_table <- function(nc, title) {
  vals <- c("reject_rate", "reject_rate_all", "sd_over_se", "fail_rate")
  if (anyNA(nc[vals])) stop(title, ": null cells with NA/NaN statistics (all-failed fits?): ",
                             paste(capture.output(print(head(nc[!complete.cases(nc[vals]), c("spec", "role", "cond", "n")]))), collapse = "\n"))
  nc$col <- paste0("sy", nc$sy, "_w", nc$waves)
  grp <- interaction(nc$role, nc$spec, nc$col, drop = TRUE, sep = "|")
  agg <- do.call(rbind, lapply(split(nc, grp), function(x) data.frame(
    row = paste(x$role[1], x$spec[1]), col = x$col[1], cells = nrow(x),
    reject_rate = mean(x$reject_rate), reject_rate_all = mean(x$reject_rate_all),
    sd_over_se = mean(x$sd_over_se), fail_rate = mean(x$fail_rate))))
  wide <- function(v, digits) {
    w <- reshape(agg[c("row", "col", v)], idvar = "row", timevar = "col", direction = "wide")
    names(w) <- sub(paste0("^", v, "\\."), "", names(w))
    w <- w[order(w$row), c("row", sort(setdiff(names(w), "row")))]
    w[-1] <- round(w[-1], digits)
    w
  }
  cat("\n==", title, "==\n")
  if (nzchar(banner)) cat(banner, "\n")
  cat("-- false-positive rate (among fits that returned an estimate) --\n"); print(wide("reject_rate", 3), row.names = FALSE)
  cat("-- false-positive rate counting failed fits as non-rejections --\n"); print(wide("reject_rate_all", 3), row.names = FALSE)
  cat("-- sd(est) / mean(SE) --\n"); print(wide("sd_over_se", 2), row.names = FALSE)
  cat("-- failure rate --\n"); print(wide("fail_rate", 4), row.names = FALSE)
  cat("-- cells pooled per entry --\n"); print(wide("cells", 0), row.names = FALSE)
}
nc <- cells[cells$null, ]
pooled_table(nc, "null cells, all")
pooled_table(nc[nc$v > 0, ], "null cells with trait variance (the set the RI-CLPM is fit on)")
worst_fail <- head(cells[order(-cells$fail_rate), c("spec", "role", "n", "waves", "sy", "co", "v", "tcor", "clxy", "fail_rate")], 10)
cat("\n== highest failure rates ==\n")
if (nzchar(banner)) cat(banner, "\n")
print(worst_fail, row.names = FALSE, digits = 3)
