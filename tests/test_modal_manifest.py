"""The run-directory guard in modal_sim.py: a run may only add to a directory
whose manifest matches the current run.

    ~/code/replicate-ai-playground/.venv/bin/python tests/test_modal_manifest.py   (from the repo root)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from modal_sim import check_manifest, manifest  # noqa: E402


def refuses(have, has_tasks, want):
    try:
        check_manifest(have, has_tasks, want, "test")
    except SystemExit:
        return True
    return False


want = manifest(20, "task,n\n1,50\n")
assert check_manifest(None, False, want, "test") is True, "empty directory: write the manifest"
assert check_manifest(dict(want), True, want, "test") is False, "same run with tasks: resume, keep manifest"
assert refuses(manifest(20, "task,n\n1,150\n"), True, want), "different grid refused"
assert refuses(manifest(40, "task,n\n1,50\n"), False, want), "different reps refused"
assert refuses(dict(want, code_sha256={k: "0" for k in want["code_sha256"]}), True, want), "different R code refused"
assert refuses(dict(want, riclpmr_sha="x"), True, want), "different package pin refused"
assert refuses(None, True, want), "task files without a manifest refused"
print("all manifest checks passed")

from modal_sim import select_tasks  # noqa: E402

grid = [{"task": str(i + 1), "n": n, "waves": w} for i, (n, w) in
        enumerate((n, w) for n in ("50", "150", "600") for w in ("6", "12", "24"))]
sub = select_tasks(grid, "", "50,150", "6,12")
assert [r["task"] for r in sub] == ["1", "2", "4", "5"], "n and waves filters keep exactly the crossed cells"
assert [r["task"] for r in select_tasks(grid, "", "", "")] == [r["task"] for r in grid], "no filter keeps every task"
for args in (("", "75", ""), ("", "", "8"), ("9", "50", "")):
    try:
        select_tasks(grid, *args)
    except SystemExit:
        continue
    raise AssertionError(f"filter {args} should have been refused")
print("all task-selection checks passed")
