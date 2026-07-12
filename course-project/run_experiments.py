#!/usr/bin/env python3
"""
run_experiments.py — Automated experiment runner for dEA vs HdEA benchmarks
Based on Li & Hu 2014 parameter settings.

Usage:
  python run_experiments.py [--mode tiny|quick|paper] [--instances all|pcb442|lin318|...]
                            [--algorithms all|dEA|HdEA_0|HdEA_1|...]
                            [--procs 4] [--runs 30]

Modes:
  tiny   — migInterval scaled 0.001x paper (quick validation, ~3s/run)
  quick  — migInterval scaled 0.01x paper (faster, ~30s/run)
  paper  — full paper settings (very long, ~hours/run)
"""

import subprocess
import os
import sys
import argparse
import json
import time
from datetime import datetime

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
MPIEXEC = "mpiexec"

# Li & Hu 2014: per-instance migration intervals and optimal solutions
PAPER_CONFIG = {
    "lin318":   {"opt": 42029,  "intervals": [30000, 40000, 50000, 60000, 70000]},
    "pcb442":   {"opt": 50778,  "intervals": [100000, 150000, 200000, 250000, 300000]},
    "u574":     {"opt": 36905,  "intervals": [60000, 70000, 80000, 90000, 100000]},
    "p654":     {"opt": 34643,  "intervals": [50000, 60000, 70000, 80000, 90000]},
    "d657":     {"opt": 48912,  "intervals": [150000, 200000, 250000, 300000, 350000]},
    "u724":     {"opt": 41910,  "intervals": [150000, 200000, 250000, 300000, 350000]},
    "rat783":   {"opt": 8806,   "intervals": [150000, 200000, 250000, 300000, 350000]},
    "dsj1000":  {"opt": 18659688, "intervals": [600000, 800000, 1000000, 1200000, 1400000]},
    "pr1002":   {"opt": 259045, "intervals": [500000, 600000, 700000, 80000, 900000]},
}

ALGORITHMS = {
    "dEA": {
        "dir": "04_dEA_island",
        "exe": "TSP_dEA_island.exe",
        # dEA: totalGen = migInterval * 2000 (same as HdEA equivalent)
        "args": lambda tsp, M, N, runs, migInterval, procs: [
            tsp, str(runs), str(migInterval * 20 * 100)
        ],
    },
    "HdEA_RingIndividual": {
        "dir": "05_HdEA", "exe": "TSP_HdEA.exe",
        "args": lambda tsp, M, N, runs, migInterval, procs: [
            tsp, str(M), str(N), str(runs), str(migInterval), "0"
        ],
    },
    "HdEA_RandomIndividual": {
        "dir": "05_HdEA", "exe": "TSP_HdEA.exe",
        "args": lambda tsp, M, N, runs, migInterval, procs: [
            tsp, str(M), str(N), str(runs), str(migInterval), "1"
        ],
    },
    "HdEA_RingColony": {
        "dir": "05_HdEA", "exe": "TSP_HdEA.exe",
        "args": lambda tsp, M, N, runs, migInterval, procs: [
            tsp, str(M), str(N), str(runs), str(migInterval), "2"
        ],
    },
    "HdEA_RandomColony": {
        "dir": "05_HdEA", "exe": "TSP_HdEA.exe",
        "args": lambda tsp, M, N, runs, migInterval, procs: [
            tsp, str(M), str(N), str(runs), str(migInterval), "3"
        ],
    },
}

SCALE_FACTORS = {
    "tiny": 0.001,
    "quick": 0.01,
    "paper": 1.0,
}


def scale_interval(interval, scale):
    val = int(interval * scale)
    return max(val, 2)


def get_group_layout(procs):
    """Determine M (groups) x N (subpops/group) from process count."""
    # For dEA: all procs form a ring (no groups needed)
    # For HdEA: split into groups
    layouts = {
        4: (2, 2),    # 2 groups x 2 subpops
        8: (2, 4),    # 2 groups x 4 subpops
        9: (3, 3),    # 3 groups x 3 subpops
        16: (4, 4),   # 4 groups x 4 subpops
        36: (6, 6),   # 6 groups x 6 subpops
        64: (8, 8),   # 8 groups x 8 subpops (paper default)
    }
    if procs in layouts:
        return layouts[procs]
    M = int(procs ** 0.5)
    while procs % M != 0:
        M -= 1
    N = procs // M
    return (M, N)


def run_experiment(alg_name, alg_info, tsp_file, M, N, num_runs,
                   mig_interval, procs, workdir):
    exe_path = os.path.join(PROJECT_DIR, alg_info["dir"], alg_info["exe"])
    args = alg_info["args"](tsp_file, M, N, num_runs,
                            mig_interval, procs)
    cmd = [MPIEXEC, "-n", str(procs), "-wdir", workdir, exe_path] + args

    print(f"  [{datetime.now().strftime('%H:%M:%S')}] {alg_name} "
          f"interval={mig_interval} runs={num_runs}")
    start = time.time()
    try:
        result = subprocess.run(cmd, cwd=workdir, capture_output=True,
                                text=True, timeout=3600)
        elapsed = time.time() - start
        if result.returncode != 0:
            print(f"    FAILED (rc={result.returncode}) after {elapsed:.0f}s")
            print(f"    stderr: {result.stderr[:200]}")
            return {"status": "failed", "elapsed": elapsed}
        return {"status": "ok", "elapsed": elapsed,
                "output": result.stdout}
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        print(f"    TIMEOUT after {elapsed:.0f}s")
        return {"status": "timeout", "elapsed": elapsed}


def main():
    parser = argparse.ArgumentParser(description="dEA vs HdEA experiment runner")
    parser.add_argument("--mode", choices=["tiny", "quick", "paper"],
                        default="tiny")
    parser.add_argument("--instances", type=str, default="pcb442",
                        help="comma-separated instance names or 'all'")
    parser.add_argument("--algorithms", type=str, default="all",
                        help="comma-separated: dEA,HdEA_RingIndividual,... or 'all'")
    parser.add_argument("--procs", type=int, default=4,
                        help="number of MPI processes")
    parser.add_argument("--runs", type=int, default=30,
                        help="number of independent runs per combination")
    parser.add_argument("--intervals", type=str, default="",
                        help="comma-separated intervals (overrides paper defaults)")
    parser.add_argument("--tsp-dir", type=str,
                        default=os.path.join(PROJECT_DIR, "tsp_instances"),
                        help="directory containing .tsp files")
    args = parser.parse_args()

    scale = SCALE_FACTORS[args.mode]
    M, N = get_group_layout(args.procs)

    if args.instances == "all":
        instances = list(PAPER_CONFIG.keys())
    else:
        instances = [s.strip() for s in args.instances.split(",")]

    if args.algorithms == "all":
        algs = list(ALGORITHMS.keys())
    else:
        algs = [s.strip() for s in args.algorithms.split(",")]

    if args.intervals:
        custom_intervals = [int(x.strip())
                            for x in args.intervals.split(",")]
    else:
        custom_intervals = None

    total_jobs = 0
    skip_list = []
    for inst in instances:
        if inst not in PAPER_CONFIG:
            print(f"WARNING: Unknown instance '{inst}', skipping")
            continue
        tsp_path = os.path.join(args.tsp_dir, f"{inst}.tsp")
        if not os.path.exists(tsp_path):
            print(f"WARNING: TSP file not found: {tsp_path}, skipping {inst}")
            skip_list.append(inst)
            continue

        cfg = PAPER_CONFIG[inst]
        intervals = custom_intervals if custom_intervals else cfg["intervals"]
        if custom_intervals:
            scaled_intervals = [max(v, 2) for v in intervals]
        else:
            scaled_intervals = [scale_interval(v, scale) for v in intervals]

        for alg in algs:
            if alg not in ALGORITHMS:
                print(f"WARNING: Unknown algorithm '{alg}', skipping")
                continue
            workdir = os.path.join(PROJECT_DIR, ALGORITHMS[alg]["dir"])
            for mig_int in scaled_intervals:
                total_jobs += 1

    instances = [i for i in instances if i not in skip_list]
    if not instances:
        print("ERROR: No valid instances found. "
              f"Place .tsp files in {args.tsp_dir}")
        return 1

    print(f"{'='*70}")
    print(f"Experiment Config")
    print(f"  Mode:       {args.mode} (scale={scale})")
    print(f"  Instances:  {', '.join(instances)}")
    print(f"  Algorithms: {', '.join(algs)}")
    print(f"  Processes:  {args.procs} (M={M}, N={N})")
    print(f"  Runs/combo: {args.runs}")
    print(f"  Total jobs: {total_jobs}")
    print(f"{'='*70}")

    results = {}
    job_idx = 0
    for inst in instances:
        cfg = PAPER_CONFIG[inst]
        tsp_path = os.path.join(args.tsp_dir, f"{inst}.tsp")
        intervals = custom_intervals if custom_intervals else cfg["intervals"]
        if custom_intervals:
            scaled_intervals = [max(v, 2) for v in intervals]
        else:
            scaled_intervals = [scale_interval(v, scale) for v in intervals]
        opt = cfg["opt"]

        results[inst] = {"opt": opt, "algorithms": {}}
        for alg in algs:
            results[inst]["algorithms"][alg] = {}
            workdir = os.path.join(PROJECT_DIR, ALGORITHMS[alg]["dir"])
            for i, mig_int in enumerate(scaled_intervals):
                job_idx += 1
                print(f"\n--- Job {job_idx}/{total_jobs}: "
                      f"{inst}/{alg} interval={mig_int} "
                      f"(paper={intervals[i]}) ---")

                res = run_experiment(
                    alg, ALGORITHMS[alg], tsp_path, M, N,
                    args.runs, mig_int, args.procs, workdir)
                results[inst]["algorithms"][alg][str(mig_int)] = res

    results_path = os.path.join(PROJECT_DIR, "experiment_results.json")
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
