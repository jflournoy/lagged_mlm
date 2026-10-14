"""
Run the factorial simulations on this machine, in the jflournoy/verse-cmdstan
container with the R library pinned to the same package snapshot as the Modal
image. Same task selection, manifest guard and resume rule as modal_sim.py.

    ~/code/replicate-ai-playground/.venv/bin/python run_local.py --reps 2500 --n 50,150 --waves 6,12 --cores 16
"""
import argparse
import csv
import io
import json
import subprocess
import sys
from pathlib import Path

from modal_sim import check_manifest, cost_weight, manifest, select_tasks, selection_record

HERE = Path(__file__).resolve().parent
IMAGE = "jflournoy/verse-cmdstan"
RLIB = Path("/tmp/claude-1000/-home-jflournoy-code-johnflournoy-science/8e8cb46d-f036-4bab-b24a-19f138b97b61/scratchpad/rlib_sim")


def docker_r(*args, check=True):
    if not RLIB.is_dir():
        raise SystemExit(f"R library not found at {RLIB}; build it with the 2024-04-01 PPM snapshot first")
    cmd = ["docker", "run", "--rm", "-e", "OPENBLAS_NUM_THREADS=1", "-e", "OMP_NUM_THREADS=1",
           "-v", f"{HERE}:/w", "-v", f"{RLIB}:/rlib", "-w", "/w", "--user", f"{_uid()}",
           "-e", "R_LIBS=/rlib", IMAGE, "Rscript", "sim_factorial.R", *map(str, args)]
    return subprocess.run(cmd, check=check)


def _uid():
    import os
    return f"{os.getuid()}:{os.getgid()}"


def write_selection(out: Path, grid: list[dict], filters: dict):
    rec = selection_record(grid, filters)
    (out / f"selection_{rec['written']}.json").write_text(json.dumps(rec))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--cores", type=int, required=True)
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--n", default="")
    ap.add_argument("--waves", default="")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    out = Path(a.out_dir) if a.out_dir else HERE / f"factorial_sims_r{a.reps}"
    out.mkdir(parents=True, exist_ok=True)
    grid_file = out / ".grid_tmp.csv"
    docker_r("grid", a.reps, f"/w/{grid_file.relative_to(HERE)}")
    grid_text = grid_file.read_text()
    grid_file.unlink()
    want = manifest(a.reps, grid_text)
    mf = out / "manifest.json"
    have = json.loads(mf.read_text()) if mf.exists() else None
    if check_manifest(have, any(out.glob("task_*.RDS")), want, str(out)):
        mf.write_text(json.dumps(want, indent=1))
        (out / "grid.csv").write_text(grid_text)
    grid = select_tasks(list(csv.DictReader(io.StringIO(grid_text))), a.only, a.n, a.waves)
    write_selection(out, grid, {"n": a.n, "waves": a.waves, "only": a.only})
    todo = [r for r in sorted(grid, key=cost_weight, reverse=True)
            if not (out / f"task_{int(r['task']):05d}.RDS").exists()]
    print(f"{len(grid)} tasks selected, {len(todo)} to run on {a.cores} cores -> {out}", flush=True)
    if not todo:
        return
    r = docker_r("runpar", a.reps, f"/w/{out.relative_to(HERE)}", a.cores, ",".join(r["task"] for r in todo), check=False)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
