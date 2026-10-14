"""
Fan the factorial lagged-MLM simulations (sim_factorial.R) out over Modal.

One task is one chunk of reps for one condition x N, with every spec fit to
the same datasets. Layout, built around the Starter plan's cap of 100
containers (not cores):

- run_task: one task per input. Each worker container reserves CPU_CORES
  physical cores and takes PROCS inputs at once (one single-threaded R process
  per hyperthread), so 99 containers run ~1,600 R processes. Modal hands each
  freed slot the next input, so long and short tasks balance across slots.
- dispatch: runs on Modal, maps run_task over the tasks (heaviest first), and
  is the only writer of results to the Volume (Modal advises against many
  concurrent writers on a v1 Volume). Because it runs remotely, the run does
  not depend on this terminal staying connected.
- The run's Volume directory holds manifest.json (reps, grid hash, R code
  hashes, package snapshot); a run may only add to a directory whose manifest
  matches, so results from different code or designs can never mix. Tasks
  already in the directory are skipped, so rerunning the same command resumes.

    modal run --detach modal_sim.py --reps 5000
    modal run --detach modal_sim.py --reps 2000 --n 50,150 --waves 6,12   # a subset of cells
    modal run modal_sim.py --reps 50 --run-name check_r50 --only 1,2,3 --max-workers 2
    modal volume get lagged-mlm-factorial <run-name> <local dir>
"""

import csv
import hashlib
import io
import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import modal

HERE = Path(__file__).parent

PPM = "https://packagemanager.posit.co/cran/__linux__/jammy/2024-04-01"
RICLPMR_SHA = "76822ab8f1d4ddf27a810c359489df69256abd86"
PKGS = ["Matrix", "lme4", "lmerTest", "lavaan", "nlme", "dplyr", "tidyr"]

r_vec = ", ".join(f'"{p}"' for p in PKGS)
image = (
    modal.Image.from_registry("rocker/r-ver:4.3.2", add_python="3.12")
    .env({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    .run_commands(
        f"Rscript -e 'options(repos = c(CRAN = \"{PPM}\")); install.packages(c({r_vec}))'",
        f"Rscript -e 'install.packages(\"https://github.com/jflournoy/riclpmr/archive/{RICLPMR_SHA}.tar.gz\", repos = NULL, type = \"source\")'",
        # install.packages() only warns when a package fails; make the build fail instead.
        f"Rscript -e 'for (p in c({r_vec}, \"riclpmr\")) {{ suppressMessages(library(p, character.only = TRUE)); "
        f"cat(p, as.character(packageVersion(p)), \"\\n\") }}'",
    )
    .add_local_file(HERE / "sim_factorial_lib.R", "/sim/sim_factorial_lib.R")
    .add_local_file(HERE / "sim_factorial.R", "/sim/sim_factorial.R")
)

app = modal.App("lagged-mlm-factorial", image=image)
VOLUME_NAME = "lagged-mlm-factorial"
vol = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

CPU_CORES = 8.0  # physical cores reserved per worker container
# Ceiling on physical cores a worker may use. Modal bills the higher of the
# reservation and actual usage, and an uncapped container can burst above its
# reservation: in the first packed check (2026-10-13) that roughly doubled the
# bill with no gain in speed. Set with SIM_CPU_LIMIT; printed on every run.
CPU_LIMIT = float(os.environ.get("SIM_CPU_LIMIT", "8.0"))
if CPU_LIMIT < CPU_CORES:
    raise SystemExit(f"SIM_CPU_LIMIT ({CPU_LIMIT}) is below the reservation ({CPU_CORES})")
PROCS = 16       # R processes per container: one per hyperthread
MEM_MB = 8192    # reservation; peak R process is ~650 MB (N = 600, 24 waves, RI-CLPM), and
                 # Modal bills the higher of reservation and use, with no limit set
MAX_WORKERS = 99  # the workspace cap is 100 containers; one is the dispatcher
# The code that determines results (generator, specs, seeds, task output).
# sim_factorial.R only launches tasks and writes files, so changing how tasks
# are run does not change what a task computes.
CODE_FILES = ["sim_factorial_lib.R"]


@app.function(cpu=1.0, memory=512)
def grid_csv(reps: int) -> str:
    r = subprocess.run(["Rscript", "/sim/sim_factorial.R", "grid", str(reps), "/tmp/grid.csv"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"grid failed:\n{r.stderr[-3000:]}")
    return Path("/tmp/grid.csv").read_text()


@app.function(cpu=(CPU_CORES, CPU_LIMIT), memory=MEM_MB, max_containers=MAX_WORKERS, scaledown_window=2, timeout=2 * 3600)
@modal.concurrent(max_inputs=PROCS)
def run_task(reps: int, task_id: int) -> dict:
    t0 = time.time()
    # Own directory per input: inputs share the container, and Modal reuses
    # warm containers.
    out = Path(tempfile.mkdtemp(prefix=f"task{task_id}_"))
    r = subprocess.run(["Rscript", "/sim/sim_factorial.R", "run", str(reps), str(out), str(task_id)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"task {task_id} failed:\n{r.stderr[-3000:]}")
    f = out / f"task_{task_id:05d}.RDS"
    data = f.read_bytes()
    shutil.rmtree(out)
    return {"task": task_id, "name": f.name, "data": data, "wall_secs": time.time() - t0}


@app.function(cpu=0.25, memory=1024, timeout=24 * 3600, volumes={"/results": vol})
def dispatch(reps: int, run_name: str, ids: list[int], max_workers: int) -> dict:
    run_task.update_autoscaler(max_containers=max_workers)
    d = Path("/results") / run_name
    done = {int(p.name[5:10]) for p in d.glob("task_*.RDS")}
    todo = [i for i in ids if i not in done]
    print(f"{len(ids)} tasks, {len(done)} already done, {len(todo)} to run", flush=True)
    started = time.time()
    last_commit, log, failed = started, [], []
    for res in run_task.map([reps] * len(todo), todo, order_outputs=False,
                            return_exceptions=True, wrap_returned_exceptions=False):
        if isinstance(res, Exception):
            failed.append(repr(res)[:2000])
            print("TASK FAILED:", repr(res)[:500], flush=True)
            continue
        tmp = d / (res["name"] + ".tmp")
        tmp.write_bytes(res["data"])
        os.replace(tmp, d / res["name"])
        log.append({"task": res["task"], "wall_secs": res["wall_secs"]})
        if time.time() - last_commit > 30:
            vol.commit()
            last_commit = time.time()
        if len(log) % 25 == 0:
            print(f"{len(log)}/{len(todo)} done, {len(failed)} failed, {time.time() - started:.0f}s", flush=True)
    summary = {"reps": reps, "elapsed_secs": time.time() - started, "n_done": len(log),
               "n_failed": len(failed), "tasks": log, "failed": failed}
    (d / f"dispatch_{int(started)}.json").write_text(json.dumps(summary))
    vol.commit()
    return {k: summary[k] for k in ("elapsed_secs", "n_done", "n_failed", "failed")}


def manifest(reps: int, grid_text: str) -> dict:
    sha = lambda b: hashlib.sha256(b).hexdigest()
    return {"reps": reps, "grid_sha256": sha(grid_text.encode()),
            "code_sha256": {f: sha((HERE / f).read_bytes()) for f in CODE_FILES},
            "ppm_snapshot": PPM, "riclpmr_sha": RICLPMR_SHA}


def check_manifest(have, has_task_files: bool, want: dict, where: str) -> bool:
    """Whether the run may use this directory; True means a manifest must be
    written. A mismatched manifest, or task files with none, stops the run."""
    if have is not None:
        if have != want:
            diff = sorted(k for k in set(have) | set(want) if have.get(k) != want.get(k))
            raise SystemExit(f"{where} holds results from a different run (differs in: {', '.join(diff)}); "
                             "use another --run-name")
        return False
    if has_task_files:
        raise SystemExit(f"{where} has task files but no manifest.json, so their origin is unknown; use another --run-name")
    return True


def select_tasks(grid: list[dict], only: str, n: str, waves: str) -> list[dict]:
    """Filter the full grid. Seeds and the manifest belong to the full grid,
    so cells left out now can be run later into the same directory with the
    same seeds. Values not in the grid are an error, not an empty selection."""
    for col, arg in (("n", n), ("waves", waves)):
        if arg:
            want = {float(v) for v in arg.split(",")}
            have = {float(r[col]) for r in grid}
            if not want <= have:
                raise SystemExit(f"--{col} values not in the grid: {sorted(want - have)} (grid has {sorted(have)})")
            grid = [r for r in grid if float(r[col]) in want]
    if only:
        wanted = {int(i) for i in only.split(",")}
        missing = wanted - {int(r["task"]) for r in grid}
        if missing:
            raise SystemExit(f"task ids not in the (filtered) grid: {sorted(missing)}")
        grid = [r for r in grid if int(r["task"]) in wanted]
    return grid


def selection_record(grid: list[dict], filters: dict) -> dict:
    """Which tasks a launch selected; the summary expects every one of them."""
    return {"written": time.strftime("%Y%m%dT%H%M%S"), "filters": filters,
            "tasks": sorted(int(r["task"]) for r in grid)}


def cost_weight(row):
    # Rough relative cost, only for ordering: grows with person-waves and reps.
    return int(row["n"]) * int(row["waves"]) * int(row["reps"])


@app.local_entrypoint()
def main(reps: int = 20, run_name: str = "", only: str = "", n: str = "", waves: str = "",
         max_workers: int = MAX_WORKERS):
    if not 1 <= max_workers <= MAX_WORKERS:
        raise SystemExit(f"--max-workers must be between 1 and {MAX_WORKERS}")
    run_name = run_name or f"factorial_sims_r{reps}"
    where = f"volume {VOLUME_NAME}:/{run_name}"
    grid_text = grid_csv.remote(reps)
    want = manifest(reps, grid_text)
    entries = [Path(e.path).name for e in vol.listdir("/")]
    if run_name in entries:
        names = [Path(e.path).name for e in vol.listdir(run_name)]
        have = json.loads(b"".join(vol.read_file(f"{run_name}/manifest.json"))) if "manifest.json" in names else None
        has_tasks = any(n.startswith("task_") for n in names)
    else:
        have, has_tasks = None, False
    if check_manifest(have, has_tasks, want, where):
        with vol.batch_upload() as up:
            up.put_file(io.BytesIO(json.dumps(want, indent=1).encode()), f"{run_name}/manifest.json")
            up.put_file(io.BytesIO(grid_text.encode()), f"{run_name}/grid.csv")
    grid = select_tasks(list(csv.DictReader(io.StringIO(grid_text))), only, n, waves)
    rec = selection_record(grid, {"n": n, "waves": waves, "only": only})
    with vol.batch_upload() as up:
        up.put_file(io.BytesIO(json.dumps(rec).encode()), f"{run_name}/selection_{rec['written']}.json")
    ids = [int(r["task"]) for r in sorted(grid, key=cost_weight, reverse=True)]
    call = dispatch.spawn(reps, run_name, ids, max_workers)
    print(f"workers reserve {CPU_CORES} cores, limit {CPU_LIMIT}, {MEM_MB} MB, {PROCS} R processes each", flush=True)
    print(f"dispatched {len(ids)} tasks to {where} with up to {max_workers} workers x {PROCS} R processes; "
          f"call {call.object_id}", flush=True)
    print(f"download: modal volume get {VOLUME_NAME} {run_name} <local dir>", flush=True)
    print(call.get())
