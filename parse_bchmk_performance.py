#!/usr/bin/env python3
"""
Parse GROMACS log files in benchmark directories and extract ns/day performance.

This script can live anywhere (e.g. /scratch/alpine/ivta1597/LCA_boltz_models)
and will read logs from the TARGET_ROOT path below.

Expected directory naming pattern (example):
  prod_md_benchmark_28np_1omp_1cpuspertask_50ps_rep2

Expected log file inside each directory:
  prod_md.log

Output:
  np_values (sorted)
  perf_rep1, perf_rep2, perf_rep3 (numpy arrays aligned to np_values)
"""

from __future__ import annotations
import re
from pathlib import Path
import numpy as np

# ---- CONFIG ----
TARGET_ROOT = Path("/scratch/alpine/ivta1597/LCA_boltz_models/benchmark/seq14/1ntomp/100ps")
LOG_NAME = "prod_md.log"
# ----------------

DIR_RE = re.compile(r"^prod_md_benchmark_(\d+)np_.*_rep(\d+)$")
PERF_RE = re.compile(r"^\s*Performance:\s*([0-9]*\.?[0-9]+)\s+([0-9]*\.?[0-9]+)\s*$")

def extract_ns_per_day(log_path: Path) -> float:
    """Extracts the ns/day performance value from a GROMACS log file.

    Args:
        log_path: Path to a prod_md.log file.

    Returns:
        ns/day, from the log's last "Performance:" line.

    Raises:
        ValueError: No "Performance:" line was found in the log.
    """
    ns_day = None
    with log_path.open("r", errors="ignore") as f:
        for line in f:
            m = PERF_RE.match(line)
            if m:
                ns_day = float(m.group(1))  # first value is (ns/day)
    if ns_day is None:
        raise ValueError(f"Could not find 'Performance:' line in {log_path}")
    return ns_day

def main():
    """Scans TARGET_ROOT for benchmark logs and prints ns/day per core count/rep."""
    if not TARGET_ROOT.exists():
        raise SystemExit(f"TARGET_ROOT does not exist: {TARGET_ROOT}")

    # Map: rep -> {np_value: performance}
    perf_by_rep: dict[int, dict[int, float]] = {}

    # Discover directories and parse logs
    for p in TARGET_ROOT.iterdir():
        if not p.is_dir():
            continue
        m = DIR_RE.match(p.name)
        if not m:
            continue

        np_val = int(m.group(1))
        rep = int(m.group(2))

        log_path = p / LOG_NAME
        if not log_path.exists():
            print(f"[WARN] Missing log: {log_path}")
            continue

        try:
            perf = extract_ns_per_day(log_path)
        except Exception as e:
            print(f"[WARN] Failed parsing {log_path}: {e}")
            continue

        perf_by_rep.setdefault(rep, {})[np_val] = perf

    if not perf_by_rep:
        raise SystemExit(f"No matching benchmark directories/logs found under {TARGET_ROOT}")

    # Union of all np values across reps, sorted
    all_np = sorted({npv for d in perf_by_rep.values() for npv in d.keys()})
    np_values = np.array(all_np, dtype=int)

    def build_rep_array(rep: int) -> np.ndarray:
        """Builds a per-np_values performance array for one rep (NaN where missing)."""
        arr = np.full(len(np_values), np.nan, dtype=float)
        rep_map = perf_by_rep.get(rep, {})
        for i, npv in enumerate(np_values):
            if npv in rep_map:
                arr[i] = rep_map[npv]
        return arr

    perf_rep1 = build_rep_array(1)
    perf_rep2 = build_rep_array(2)
    perf_rep3 = build_rep_array(3)
    perf_rep4 = build_rep_array(4)
    perf_rep5 = build_rep_array(5)
    perf_rep6 = build_rep_array(6)
    perf_rep7 = build_rep_array(7)
    perf_rep8 = build_rep_array(8)
    perf_rep9 = build_rep_array(9)
    perf_rep10 = build_rep_array(10)
    # Print results
    print("TARGET_ROOT =", str(TARGET_ROOT))
    print("cores       =", np_values)
    print("perf_rep1   =", perf_rep1)
    print("perf_rep2   =", perf_rep2)
    print("perf_rep3   =", perf_rep3)
    print("perf_rep4   =", perf_rep4)
    print("perf_rep5   =", perf_rep5)
    print("perf_rep6   =", perf_rep6)
    print("perf_rep7   =", perf_rep7)
    print("perf_rep8   =", perf_rep8)
    print("perf_rep9   =", perf_rep9)
    print("perf_rep10   =", perf_rep10)

    # Optional: save for plotting later
    # out = Path.cwd() / "benchmark_performance_seq14_1ntomp.npz"
    # np.savez(out, np_values=np_values, perf_rep1=perf_rep1, perf_rep2=perf_rep2, perf_rep3=perf_rep3)
    # print(f"Saved: {out}")

if __name__ == "__main__":
    main()
