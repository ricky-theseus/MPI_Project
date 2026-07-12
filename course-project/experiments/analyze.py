#!/usr/bin/env python3
"""
analyze.py — Read experiment output and generate 4 figures + t-test report

Output:
  experiments/output/<instance>/analysis/
    fig1_convergence.png
    fig2_boxplot.png
    fig3_ttest.png
    fig4_runtime.png
    report.txt
"""

import os
import sys
import argparse
import glob
import json
from math import sqrt
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("WARNING: matplotlib not installed. Install: pip install matplotlib")


ALG_ORDER = ["dEA", "HdEA_RingIndividual", "HdEA_RingColony"]
ALG_LABELS = {"dEA": "dEA",
              "HdEA_RingIndividual": "HdEA RingIndiv",
              "HdEA_RingColony": "HdEA RingColony (MP)"}
COLORS = {"dEA": "#3498db",
          "HdEA_RingIndividual": "#2ecc71",
          "HdEA_RingColony": "#e74c3c"}
OPTIMAL = {
    "lin318": 42029, "pcb442": 50778, "u574": 36905, "p654": 34643,
    "d657": 48912, "u724": 41910, "rat783": 8806,
    "dsj1000": 18659688, "pr1002": 259045
}


def read_results(results_path):
    """Parse results.txt. Returns list of (run, best, elapsed)."""
    runs = []
    if not os.path.exists(results_path):
        return runs
    with open(results_path) as f:
        header = True
        for line in f:
            line = line.strip()
            if not line:
                continue
            if header and "run" in line.lower():
                header = False
                continue
            if header:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    runs.append((int(parts[0]), float(parts[1]),
                                 float(parts[2]) if len(parts) > 2 else 0))
                except ValueError:
                    continue
    return runs


def read_convergence(conv_path):
    """Parse run_N_conv.txt. Returns list of (gen, best, time)."""
    trace = []
    if not os.path.exists(conv_path):
        return trace
    with open(conv_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    trace.append((int(parts[0]), float(parts[1]),
                                  float(parts[2]) if len(parts) > 2 else 0))
                except ValueError:
                    continue
    return trace


def t_test(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return 0, 1.0
    ma = sum(a) / na
    mb = sum(b) / nb
    va = sum((x - ma) ** 2 for x in a) / (na - 1)
    vb = sum((x - mb) ** 2 for x in b) / (nb - 1)
    se = sqrt(va / na + vb / nb)
    if se == 0:
        return 0, 1.0
    t = (ma - mb) / se
    df_num = (va / na + vb / nb) ** 2
    df_den = (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    df = max(1, df_num / df_den)
    x = df / (df + t * t)
    import math
    try:
        lbeta = math.lgamma(df / 2 + 0.5) - math.lgamma(df / 2) - math.lgamma(0.5)
        p = math.exp(math.log(x) * (df / 2) + math.log(1 - x) * 0.5 - lbeta)
        p = min(p, 1 - p) * 2
    except (ValueError, OverflowError):
        p = 0.0 if abs(t) > 2 else 1.0
    return t, p


def collect_data(output_base, instance):
    """Collect all data for one instance across all intervals and algorithms."""
    inst_dir = os.path.join(output_base, instance)
    if not os.path.exists(inst_dir):
        return {}

    interval_dirs = sorted(glob.glob(os.path.join(inst_dir, "interval_*")))
    data = {}
    for intv_dir in interval_dirs:
        intv_name = os.path.basename(intv_dir)
        data[intv_name] = {}
        for alg in ALG_ORDER:
            alg_dir = os.path.join(intv_dir, alg)
            results_path = os.path.join(alg_dir, "results.txt")
            runs = read_results(results_path)
            if not runs:
                continue

            values = [r[1] for r in runs]
            times = [r[2] for r in runs]
            conv_traces = []
            for f in sorted(glob.glob(os.path.join(alg_dir, "run_*_conv.txt"))):
                conv_traces.append(read_convergence(f))

            data[intv_name][alg] = {
                "values": values, "times": times,
                "mean": sum(values) / len(values),
                "std": sqrt(sum((v - sum(values)/len(values))**2
                                for v in values) / max(1, len(values)-1)),
                "best": min(values),
                "n": len(values),
                "convergence": conv_traces,
            }
    return data


def plot_convergence(data, opt, out_dir):
    """Figure 1: Convergence curves (best vs generation)."""
    if not HAS_MPL or not data:
        return

    interval = list(data.keys())[0]
    algs_data = data[interval]

    fig, ax = plt.subplots(figsize=(10, 6))
    for alg in ALG_ORDER:
        if alg not in algs_data:
            continue
        traces = algs_data[alg]["convergence"]
        if not traces:
            continue
        max_len = max(len(t) for t in traces)
        all_gens = {t[i][0] for t in traces for i in range(len(t))
                    if t[i][0] <= 200000}
        sorted_gens = sorted(all_gens)

        avg = []
        for g in sorted_gens:
            vals = []
            for t in traces:
                for pt in t:
                    if pt[0] == g:
                        vals.append(pt[1])
                        break
            if vals:
                avg.append((g, sum(vals) / len(vals)))
        if avg:
            gs, bs = zip(*avg)
            ax.plot(gs, bs, color=COLORS.get(alg, "#333"),
                    label=ALG_LABELS.get(alg, alg), linewidth=1.5)

    if opt:
        ax.axhline(y=opt, color="gray", linestyle="--",
                   linewidth=1, label=f"Optimal ({opt})")

    ax.set_xlabel("Generation", fontsize=12)
    ax.set_ylabel("Best Tour Length", fontsize=12)
    ax.set_title(f"Convergence Curves — pcb442 (interval={interval})",
                 fontsize=14)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    path = os.path.join(out_dir, "fig1_convergence.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_boxplot(data, opt, out_dir):
    """Figure 2: Box plot of final results."""
    if not HAS_MPL or not data:
        return

    interval = list(data.keys())[0]
    algs_data = data[interval]

    all_vals = []
    labels = []
    colors = []
    for alg in ALG_ORDER:
        if alg not in algs_data:
            continue
        all_vals.append(algs_data[alg]["values"])
        labels.append(ALG_LABELS.get(alg, alg))
        colors.append(COLORS.get(alg, "#333"))

    fig, ax = plt.subplots(figsize=(9, 5))
    bp = ax.boxplot(all_vals, labels=labels, patch_artist=True,
                    widths=0.5)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)

    if opt:
        ax.axhline(y=opt, color="gray", linestyle="--",
                   linewidth=1, label=f"Optimal={opt}")
        ax.legend(fontsize=9)

    ax.set_ylabel("Best Tour Length", fontsize=12)
    ax.set_title("Final Results Distribution — pcb442 (10 runs)",
                 fontsize=14)
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    path = os.path.join(out_dir, "fig2_boxplot.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_ttest(data, out_dir):
    """Figure 3: t-test significance matrix."""
    if not HAS_MPL or not data:
        return

    interval = list(data.keys())[0]
    algs_data = data[interval]
    algs = [a for a in ALG_ORDER if a in algs_data]
    n = len(algs)

    pmat = [[1.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                _, p = t_test(algs_data[algs[i]]["values"],
                               algs_data[algs[j]]["values"])
                pmat[i][j] = p

    fig, ax = plt.subplots(figsize=(7, 6))
    cmap = plt.cm.RdYlGn
    im = ax.imshow(pmat, cmap=cmap, vmin=0, vmax=1, aspect="auto")

    labels = [ALG_LABELS.get(a, a) for a in algs]
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(labels, fontsize=9)

    for i in range(n):
        for j in range(n):
            v = pmat[i][j]
            text = f"{v:.3f}" if v < 1 else "—"
            color = "white" if v < 0.3 else "black"
            ax.text(j, i, text, ha="center", va="center",
                    fontsize=7, color=color, fontweight="bold")

    ax.set_title("t-test p-value Matrix (significant if p<0.05)", fontsize=13)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("p-value", fontsize=10)
    plt.tight_layout()
    path = os.path.join(out_dir, "fig3_ttest.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_runtime(data, out_dir):
    """Figure 4: Runtime comparison."""
    if not HAS_MPL or not data:
        return

    interval = list(data.keys())[0]
    algs_data = data[interval]

    algs = [a for a in ALG_ORDER if a in algs_data]
    means = [sum(algs_data[a]["times"]) / len(algs_data[a]["times"])
             for a in algs]
    stds = [sqrt(sum((t - m) ** 2 for t in algs_data[a]["times"])
                  / max(1, len(algs_data[a]["times"]) - 1))
            for a, m in zip(algs, means)]
    labels = [ALG_LABELS.get(a, a) for a in algs]
    colors = [COLORS.get(a, "#333") for a in algs]

    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.bar(labels, means, yerr=stds, color=colors, alpha=0.7,
                  capsize=5)
    ax.set_ylabel("Runtime (seconds)", fontsize=12)
    ax.set_title("Average Runtime — pcb442 (10 runs)", fontsize=14)
    for bar, m in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"{m:.1f}s", ha="center", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    path = os.path.join(out_dir, "fig4_runtime.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


def write_report(data, opt, out_dir):
    """Write text report with statistics and t-test results."""
    interval = list(data.keys())[0]
    algs_data = data[interval]
    algs = [a for a in ALG_ORDER if a in algs_data]

    report_path = os.path.join(out_dir, "report.txt")
    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("dEA vs HdEA Experiment Report\n")
        f.write(f"Instance: pcb442 (optimal={opt})\n")
        f.write(f"Interval: {interval}\n\n")

        f.write(f"{'Algorithm':<25} {'Mean':>10} {'Std':>10} "
                f"{'Best':>10} {'Time(s)':>10} {'#Runs':>6}\n")
        f.write("-" * 71 + "\n")
        for alg in algs:
            d = algs_data[alg]
            f.write(f"{alg:<25} {d['mean']:>10.1f} {d['std']:>10.1f} "
                    f"{d['best']:>10.1f} "
                    f"{sum(d['times'])/len(d['times']):>10.1f} "
                    f"{d['n']:>6}\n")

        f.write("\n--- t-test Results (p < 0.05 = significant) ---\n\n")
        significant = []
        for i in range(len(algs)):
            for j in range(i + 1, len(algs)):
                a, b = algs[i], algs[j]
                t, p = t_test(algs_data[a]["values"],
                               algs_data[b]["values"])
                sig = "SIGNIFICANT" if p < 0.05 else "not sig"
                if p < 0.05:
                    winner = a if algs_data[a]["mean"] < algs_data[b]["mean"] else b
                    significant.append((winner, a if winner == b else b, p))
                f.write(f"  {a} vs {b}: p={p:.4f} ({sig})\n")

        if significant:
            f.write("\nSignificant differences:\n")
            for winner, loser, p in sorted(significant, key=lambda x: x[2]):
                f.write(f"  {winner} > {loser} (p={p:.4f})\n")
        else:
            f.write("\nNo significant differences found.\n")

    print(f"  Saved: {report_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instance", default="pcb442")
    parser.add_argument("--output-base",
                        default=os.path.join(
                            os.path.dirname(os.path.dirname(__file__)),
                            "experiments", "output"))
    args = parser.parse_args()

    output_base = os.path.abspath(args.output_base)
    instance = args.instance
    opt = OPTIMAL.get(instance)

    print(f"Analyzing: {instance} (optimal={opt})")
    print(f"Data dir: {output_base}")

    data = collect_data(output_base, instance)
    if not data:
        print("ERROR: No data found. Run run_all.py first.")
        return 1

    ana_dir = os.path.join(output_base, instance, "analysis")
    os.makedirs(ana_dir, exist_ok=True)

    plot_convergence(data, opt, ana_dir)
    plot_boxplot(data, opt, ana_dir)
    plot_ttest(data, ana_dir)
    plot_runtime(data, ana_dir)
    write_report(data, opt, ana_dir)

    print(f"\nAll outputs in: {ana_dir}")


if __name__ == "__main__":
    main()
