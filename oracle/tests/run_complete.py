"""Run independent native validation cases in isolated Python processes."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CASES = (
    [
        "uniform11",
        "uniform59",
        "phase11_character",
        "phase11_cofactor",
        "prepare11_character",
        "prepare11_cofactor",
        "prepare59_character",
        "prepare59_cofactor",
        "selector",
        "signed_oop",
        "group11_M1",
        "group11_M2",
        "group59_M1",
    ]
    + [f"failure{i}" for i in range(6)]
    + ["failure_mid_chain"]
)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--jobs", type=int, default=1)
parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
args = parser.parse_args()
assert args.jobs > 0
out = ROOT / "results/complete_validation"
out.mkdir(parents=True, exist_ok=True)


def run(case):
    started = time.monotonic()
    with (out / f"{case}.log").open("w") as log:
        r = subprocess.run(
            [sys.executable, str(ROOT / "tests/test_complete.py"), "--case", case],
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    row = dict(case=case, returncode=r.returncode, seconds=time.monotonic() - started)
    print(json.dumps(row), flush=True)
    return row


with ThreadPoolExecutor(max_workers=args.jobs) as pool:
    results = [
        f.result() for f in as_completed([pool.submit(run, c) for c in args.cases])
    ]
summary = {"runs": results}
(out / "last_run.json").write_text(json.dumps(summary, indent=2) + "\n")
sys.exit(0 if all(r["returncode"] == 0 for r in results) else 1)
