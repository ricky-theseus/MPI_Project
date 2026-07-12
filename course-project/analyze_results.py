#!/usr/bin/env python3
"""
analyze_results.py — Parse HdEA/dEA result files and compute t-test statistics
Based on Li & Hu 2014 methodology: 30 runs, t-test 95% confidence

Usage:
  python analyze_results.py [--dir <directory>] [--output <file>]
"""

import os
import re
import sys
import argparse
import json
from math import sqrt
from collections import defaultdict

PAPER_CONFIG = {
    "lin318":   {"opt": 42029,  "intervals": [30000, 40000, 50000, 60000, 70000]},
    "pcb442":   {"opt": 50778,  "intervals": [100000, 150000, 200000, 250000, 300000]},
    "u574":     {"opt": 36905,  "intervals": [60000, 70000, 80000, 90000, 100000]},
    "p654":     {"opt": 34643,  "intervals": [50000, 60000, 70000, 80000, 90000]},
    "d657":     {"opt": 48912,  "intervals": [150000, 200000, 250000, 300000, 350000]},
    "u724":     {"opt": 41910,  "intervals": [150000, 200000, 250000, 300000, 350000]},
    "rat783":   {"opt": 8806,   "intervals": [150000, 200000, 250000, 300000, 350000]},
    "dsj1000":  {"opt": 18659688, "intervals": [600000, 800000, 1000000, 1200000, 1400000]},
    "pr1002":   {"opt": 259045, "intervals": [500000, 600000, 700000, 800000, 900000]},
}

ALG_PATTERNS = {
    "dEA":          r"results_dEA_island",
    "RingIndividual":   r"results_HdEA_RingIndividual",
    "RandomIndividual": r"results_HdEA_RandomIndividual",
    "RingColony":       r"results_HdEA_RingColony",
    "RandomColony":     r"results_HdEA_RandomColony",
}

ALG_DIRS = {
    "dEA": "04_dEA_island",
    "RingIndividual": "05_HdEA",
    "RandomIndividual": "05_HdEA",
    "RingColony": "05_HdEA",
    "RandomColony": "05_HdEA",
}


def parse_result_file(path):
    """Parse a result file. Returns list of (run, best, elapsed, globalMigs)."""
    results = []
    try:
        with open(path, "r") as f:
            header = True
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if header and ("run\t" in line or "run " in line.split()[0]):
                    header = False
                    continue
                if header:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    parts = line.split()
                if len(parts) >= 2:
                    try:
                        run = int(parts[0])
                        best = float(parts[1])
                        elapsed = float(parts[2]) if len(parts) > 2 else 0
                        globalMigs = int(parts[3]) if len(parts) > 3 else 0
                        results.append((run, best, elapsed, globalMigs))
                    except (ValueError, IndexError):
                        continue
    except FileNotFoundError:
        pass
    return results


def t_test(samples_a, samples_b):
    """Independent two-sample t-test (Welch's). Returns t-stat, p-value approx."""
    n_a = len(samples_a)
    n_b = len(samples_b)
    if n_a < 2 or n_b < 2:
        return 0, 1.0

    mean_a = sum(samples_a) / n_a
    mean_b = sum(samples_b) / n_b

    var_a = sum((x - mean_a) ** 2 for x in samples_a) / (n_a - 1)
    var_b = sum((x - mean_b) ** 2 for x in samples_b) / (n_b - 1)

    se = sqrt(var_a / n_a + var_b / n_b)
    if se == 0:
        return 0, 1.0

    t = (mean_a - mean_b) / se

    df_num = (var_a / n_a + var_b / n_b) ** 2
    df_den = ((var_a / n_a) ** 2 / (n_a - 1) +
              (var_b / n_b) ** 2 / (n_b - 1))
    df = max(1, df_num / df_den)

    # Approximate p-value from t-distribution
    # Using a simple approximation (valid for df > 10)
    x = df / (df + t * t)
    p = _incomplete_beta(df / 2, 0.5, x) if abs(t) < 100 else 0.0

    return t, p


def _incomplete_beta(a, b, x):
    """Continued fraction approximation of regularized incomplete beta."""
    if x == 0:
        return 0.0
    if x == 1:
        return 1.0

    import math
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)

    front = math.exp(math.log(x) * a + math.log(1.0 - x) * b - lbeta) / a

    f = 1.0
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1.0)
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d

    for m in range(1, 200):
        m2 = 2 * m

        # even step
        d = 1.0 + m * (b - m) * x / ((a + m2 - 1) * (a + m2))
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + m * (a - m) * x / ((a + m2) * (a + m2 + 1))
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h = h * d * c
        d = 1.0 + (m + 1) * (b - m - 1) * x / ((a + m2 + 1) * (a + m2 + 2))
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + (m + 1) * (a - m - 1) * x / ((a + m2 + 2) * (a + m2 + 3))
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = h * d * c
        h = h * delta
        f = f * delta
        c = d

        if abs(delta - 1.0) < 1e-15:
            break

    return front * (f - 1.0)


def analyze_all(project_dir, scale=1.0):
    """Scan all result files and produce comparison tables."""
    print("=" * 80)
    print("dEA vs HdEA Experiment Analysis (Li & Hu 2014)")
    print(f"Scale factor: {scale} | t-test: 95% confidence")
    print("=" * 80)

    all_data = {}

    for inst_name, cfg in PAPER_CONFIG.items():
        all_data[inst_name] = {"opt": cfg["opt"], "algs": {}}
        for alg_name, dir_name in ALG_DIRS.items():
            workdir = os.path.join(project_dir, dir_name)
            results = parse_result_file(
                os.path.join(workdir, f"results_HdEA_{alg_name}.txt"
                             if alg_name != "dEA"
                             else "results_dEA_island.txt"))
            if results:
                all_data[inst_name]["algs"][alg_name] = results

    # Print per-instance comparison
    for inst_name, data in all_data.items():
        opt = data["opt"]
        algs = data["algs"]
        if not algs:
            continue

        print(f"\n--- {inst_name} (optimal={opt}) ---")
        print(f"{'Algorithm':<20} {'Mean':>10} {'Std':>10} "
              f"{'Best':>10} {'Gap%':>8} {'#Runs':>6}")
        print("-" * 70)

        summaries = {}
        for alg_name, results in algs.items():
            values = [r[1] for r in results]
            if not values:
                continue
            mean = sum(values) / len(values)
            std = sqrt(sum((v - mean)**2 for v in values) /
                       max(1, len(values) - 1))
            best = min(values)
            gap = (mean - opt) / opt * 100
            summaries[alg_name] = {
                "mean": mean, "std": std, "best": best,
                "gap": gap, "n": len(values), "values": values
            }
            print(f"{alg_name:<20} {mean:>10.1f} {std:>10.1f} "
                  f"{best:>10.1f} {gap:>7.2f}% {len(values):>5}")

        # t-test matrix
        alg_list = list(summaries.keys())
        if len(alg_list) > 1:
            print(f"\n{'':>20}", end="")
            for a in alg_list:
                print(f"{a[:10]:>12}", end="")
            print()

            for a in alg_list:
                print(f"{a:<20}", end="")
                for b in alg_list:
                    if a == b:
                        print(f"{'---':>12}", end="")
                    else:
                        t, p = t_test(summaries[a]["values"],
                                       summaries[b]["values"])
                        sig = "+" if p < 0.05 else " "
                        if t > 0:
                            print(f"{sig}WIN{p*100:6.1f}%", end="")
                        else:
                            print(f"{sig}LOSE{p*100:5.1f}%", end="")
                print()

        # Print details: which pairs are significant
        print("\nSignificant differences (p < 0.05):")
        found = False
        for i in range(len(alg_list)):
            for j in range(i + 1, len(alg_list)):
                a, b = alg_list[i], alg_list[j]
                t, p = t_test(summaries[a]["values"],
                               summaries[b]["values"])
                if p < 0.05:
                    winner = a if summaries[a]["mean"] < summaries[b]["mean"] else b
                    loser = b if winner == a else a
                    w_mean = summaries[winner]["mean"]
                    l_mean = summaries[loser]["mean"]
                    print(f"  {winner} > {loser} "
                          f"(means {w_mean:.1f} vs {l_mean:.1f}, "
                          f"p={p:.4f})")
                    found = True
        if not found:
            print("  None")

    return all_data


def main():
    parser = argparse.ArgumentParser(
        description="Analyze dEA vs HdEA experiment results")
    parser.add_argument("--dir", type=str,
                        default=os.path.dirname(os.path.abspath(__file__)),
                        help="project directory containing result files")
    parser.add_argument("--output", type=str, default="",
                        help="output JSON file for results")
    parser.add_argument("--scale", type=float, default=1.0,
                        help="scale factor used in experiment")
    args = parser.parse_args()

    data = analyze_all(args.dir, args.scale)

    if args.output:
        # Convert to serializable format (remove raw values)
        serializable = {}
        for inst, idata in data.items():
            serializable[inst] = {"opt": idata["opt"], "algs": {}}
            for alg, adata in idata["algs"].items():
                serializable[inst]["algs"][alg] = {
                    "mean": sum(v[1] for v in adata) / len(adata),
                    "n": len(adata)
                }
        with open(args.output, "w") as f:
            json.dump(serializable, f, indent=2)
        print(f"\nSummary saved to {args.output}")


if __name__ == "__main__":
    main()
