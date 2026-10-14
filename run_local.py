"""
Run the factorial simulations on this machine, in the lagged-mlm image built
from docker/Dockerfile (the same environment modal_sim.py builds on Modal).
Same task selection, manifest guard and resume rule as modal_sim.py.

    ~/code/replicate-ai-playground/.venv/bin/python run_local.py --reps 2500 --n 50,150 --waves 6,12 --cores 16
"""
import argparse
import csv
import io
import json
import os
import subprocess
import sys
from pathlib import Path

from modal_sim import check_manifest, cost_weight, manifest, select_tasks, selection_record

HERE = Path(__file__).resolve().parent
IMAGE = "lagged-mlm"


def check_image():
    r = subprocess.run(["docker", "image", "inspect", IMAGE], capture_output=True)
    if r.returncode != 0:
        raise SystemExit(f"docker image '{IMAGE}' not found; build it with: docker build -t {IMAGE} {HERE / 'docker'}")


def docker_r(*args, check=True):
    cmd = ["docker", "run", "--rm", "-v", f"{HERE}:/w", "-w", "/w", "--user", f"{os.getuid()}:{os.getgid()}",
           "-e", "HOME=/tmp", IMAGE, "Rscript", "sim_factorial.R", *map(str, args)]
    return subprocess.run(cmd, check=check)


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
    check_image()
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
