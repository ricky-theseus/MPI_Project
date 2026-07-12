#!/usr/bin/env python3
"""
run_all.py — One-click experiment runner for dEA vs HdEA
Creates output in experiments/output/<instance>/interval_<N>/<algorithm>/

Usage:
  python run_all.py [--instances pcb442,lin318] [--algorithms dEA,HdEA_RingIndividual,...]
                    [--mode tiny|quick|paper] [--runs 10] [--procs 4]
"""

import subprocess
import os
import sys
import argparse
import time
from datetime import datetime

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_BASE = os.path.join(PROJECT_DIR, "experiments", "output")

ALGORITHMS = {
    "dEA": {
        "dir": "04_dEA_island", "exe": "TSP_dEA_island.exe",
        "extra_args": [],
    },
    "HdEA_RingIndividual": {
        "dir": "05_HdEA", "exe": "TSP_HdEA.exe",
        "extra_args": ["0"],
    },
    "HdEA_RandomIndividual": {
        "dir": "05_HdEA", "exe": "TSP_HdEA.exe",
        "extra_args": ["1"],
    },
    "HdEA_RingColony": {
        "dir": "05_HdEA", "exe": "TSP_HdEA.exe",
        "extra_args": ["2"],
    },
    "HdEA_RandomColony": {
        "dir": "05_HdEA", "exe": "TSP_HdEA.exe",
        "extra_args": ["3"],
    },
}

PAPER_INTERVALS = {
    "lin318":  [30000, 40000, 50000, 60000, 70000],
    "pcb442":  [100000, 150000, 200000, 250000, 300000],
    "u574":    [60000, 70000, 80000, 90000, 100000],
    "p654":    [50000, 60000, 70000, 80000, 90000],
    "d657":    [150000, 200000, 250000, 300000, 350000],
    "u724":    [150000, 200000, 250000, 300000, 350000],
    "rat783":  [150000, 200000, 250000, 300000, 350000],
    "dsj1000": [600000, 800000, 1000000, 1200000, 1400000],
    "pr1002":  [500000, 600000, 700000, 800000, 900000],
}

SCALE = {"tiny": 0.001, "quick": 0.01, "paper": 1.0}


def get_group_layout(procs):
    layouts = {4: (2, 2), 8: (2, 4), 9: (3, 3), 16: (4, 4),
               36: (6, 6), 64: (8, 8)}
    if procs in layouts:
        return layouts[procs]
    for m in range(int(procs ** 0.5), 0, -1):
        if procs % m == 0:
            return (m, procs // m)
    return (1, procs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instances", default="pcb442")
    parser.add_argument("--algorithms", default="all")
    parser.add_argument("--mode", default="tiny", choices=["tiny", "quick", "paper"])
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--procs", type=int, default=4)
    parser.add_argument("--intervals", default="")
    args = parser.parse_args()

    scale = SCALE[args.mode]
    M, N = get_group_layout(args.procs)

    instances = [s.strip() for s in args.instances.split(",")]
    if args.algorithms == "all":
        algs = list(ALGORITHMS.keys())
    else:
        algs = [s.strip() for s in args.algorithms.split(",")]

    tsp_dir = os.path.join(PROJECT_DIR, "tsp_instances")

    total = 0
    jobs = []
    for inst in instances:
        intervals = PAPER_INTERVALS.get(inst, [100000])
        if args.intervals:
            intervals = [int(x) for x in args.intervals.split(",")]
            scaled = [max(v, 2) for v in intervals]
        else:
            scaled = [max(int(v * scale), 2) for v in intervals]

        for alg_name in algs:
            alg = ALGORITHMS[alg_name]
            for i, mig_int in enumerate(scaled):
                out_dir = os.path.join(
                    OUTPUT_BASE, inst,
                    f"interval_{intervals[i]}",
                    alg_name)
                os.makedirs(out_dir, exist_ok=True)
                jobs.append((alg_name, alg, inst, mig_int, out_dir, M, N))
                total += 1

    print(f"{'='*60}")
    print(f"Run All Experiments")
    print(f"  Mode: {args.mode} (scale={scale})")
    print(f"  Instances: {instances}")
    print(f"  Algorithms: {algs}")
    print(f"  Procs: {args.procs} (M={M}, N={N})")
    print(f"  Runs per combo: {args.runs}")
    print(f"  Total jobs: {total}")
    print(f"  Output: {OUTPUT_BASE}")
    print(f"{'='*60}")

    for idx, (alg_name, alg, inst, mig_int, out_dir, M, N) in enumerate(jobs):
        tsp_path = os.path.join(tsp_dir, f"{inst}.tsp")
        exe_path = os.path.join(PROJECT_DIR, alg["dir"], alg["exe"])
        workdir = os.path.join(PROJECT_DIR, alg["dir"])

        cmd = ["mpiexec", "-n", str(args.procs), "-wdir", workdir, exe_path,
               tsp_path]

        if alg_name == "dEA":
            cmd += [str(args.runs),
                    str(mig_int * 20 * 100),
                    out_dir]
        else:
            cmd += [str(M), str(N), str(args.runs),
                    str(mig_int)] + alg["extra_args"] + [out_dir]

        t0 = time.time()
        tag = f"[{idx+1}/{total}] {inst}/{alg_name} mig={mig_int}"
        print(f"{tag}  [{datetime.now().strftime('%H:%M:%S')}]", flush=True)

        try:
            r = subprocess.run(cmd, cwd=workdir, capture_output=True,
                               text=True, timeout=7200)
            elapsed = time.time() - t0
            ok = r.returncode == 0
            status = "OK" if ok else f"FAIL(rc={r.returncode})"
            print(f"  {status}  {elapsed:.0f}s", flush=True)
            if not ok:
                print(f"  stderr: {r.stderr[:200]}")
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT  {time.time()-t0:.0f}s", flush=True)
        except Exception as e:
            print(f"  ERROR: {e}", flush=True)

    print(f"\nDone. Output in: {OUTPUT_BASE}")


if __name__ == "__main__":
    main()
