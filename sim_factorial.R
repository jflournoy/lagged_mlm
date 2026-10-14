# Command-line front end for the factorial simulations (functions in
# sim_factorial_lib.R). Built to be fanned out by modal_sim.py, but runs the
# same anywhere.
#
#   Rscript sim_factorial.R grid <reps_per_cell> <out.csv>
#   Rscript sim_factorial.R run  <reps_per_cell> <out_dir> <task_id>[,<task_id>...]
#   Rscript sim_factorial.R runpar <reps_per_cell> <out_dir> <cores> <task_id>[,<task_id>...]
#
# `run` loads the packages once, then runs each listed task (one chunk of one
# cell) and saves it as <out_dir>/task_<id>.RDS. The package load time is
# stored in every file (`load_secs`) so it can be separated from compute.
# `runpar` does the same with <cores> separate R worker processes (a PSOCK
# cluster, not forks: forked workers segfaulted at random in 45 of 2,400 tasks
# of the first 2,500-rep run, and those tasks reran cleanly unforked), each task
# handed to the next free worker in the order given (run_local.py passes
# heaviest first). A task that errors or delivers nothing is reported and the
# run exits nonzero after the rest.

t_start <- proc.time()[["elapsed"]]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1) stop("run this with Rscript so the script's own directory can be found")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
if (length(args) < 3) stop("usage: sim_factorial.R grid <reps_per_cell> <out.csv> | run <reps_per_cell> <out_dir> <ids>")
reps_per_cell <- as.integer(args[[2]])
if (is.na(reps_per_cell) || reps_per_cell < 1) stop("reps_per_cell must be a positive integer")

source(file.path(script_dir, "sim_factorial_lib.R"))
load_secs <- proc.time()[["elapsed"]] - t_start
grid <- make_grid(reps_per_cell)

parse_ids <- function(x) {
  ids <- as.integer(strsplit(x, ",")[[1]])
  if (anyNA(ids) || !all(ids %in% grid$task)) stop("task ids not in the grid: ", x)
  ids
}
save_task <- function(id, out_dir) {
  res <- run_task(grid[grid$task == id, ])
  res$load_secs <- load_secs
  res$host <- Sys.info()[["nodename"]]
  f <- file.path(out_dir, sprintf("task_%05d.RDS", id))
  saveRDS(res, paste0(f, ".tmp"))
  if (!file.rename(paste0(f, ".tmp"), f)) stop("could not move ", f, ".tmp into place")
  f
}

if (args[[1]] == "grid") {
  write.csv(grid, args[[3]], row.names = FALSE)
} else if (args[[1]] == "run") {
  if (length(args) != 4) stop("run needs <reps_per_cell> <out_dir> <task_ids>")
  dir.create(args[[3]], showWarnings = FALSE, recursive = TRUE)
  for (id in parse_ids(args[[4]])) save_task(id, args[[3]])
} else if (args[[1]] == "runpar") {
  if (length(args) != 5) stop("runpar needs <reps_per_cell> <out_dir> <cores> <task_ids>")
  out_dir <- args[[3]]
  cores <- as.integer(args[[4]])
  if (is.na(cores) || cores < 1) stop("cores must be a positive integer")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ids <- parse_ids(args[[5]])
  started <- Sys.time()
  cl <- parallel::makePSOCKcluster(cores, outfile = "")
  parallel::clusterExport(cl, c("script_dir", "reps_per_cell", "out_dir"), envir = environment())
  parallel::clusterEvalQ(cl, {
    t0 <- proc.time()[["elapsed"]]
    suppressMessages(source(file.path(script_dir, "sim_factorial_lib.R")))
    load_secs <- proc.time()[["elapsed"]] - t0
    grid <- make_grid(reps_per_cell)
    NULL
  })
  parallel::clusterExport(cl, "save_task", envir = environment())
  out <- parallel::parLapplyLB(cl, ids, function(id) {
    tryCatch({
      f <- save_task(id, out_dir)
      message(format(Sys.time(), "%H:%M:%S"), " ", basename(f))
      f
    }, error = function(e) structure(conditionMessage(e), class = "task_error"))
  }, chunk.size = 1)
  parallel::stopCluster(cl)
  delivered <- vapply(out, function(x) is.character(x) && !inherits(x, "task_error") &&
                        file.exists(x), logical(1))
  message(sum(delivered), " tasks done, ", sum(!delivered), " failed, elapsed ", format(Sys.time() - started))
  if (any(!delivered)) {
    stop(sum(!delivered), " tasks failed:\n", paste(ids[!delivered], vapply(out[!delivered], function(x) {
      if (is.null(x)) "no result" else as.character(x)
    }, character(1)), collapse = "\n"))
  }
} else {
  stop("unknown mode: ", args[[1]])
}
