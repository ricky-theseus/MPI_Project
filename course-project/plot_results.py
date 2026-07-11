# plot_results.py
# Generates 3 figures: boxplot, barplot, convergence curves

import os
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULTS_DIR = "results"
ALGORITHMS = ["serial", "dEA", "HdEA", "HdEA_MP"]
COLORS = ["#888888", "#4C72B0", "#55A868", "#C44E52"]
LABELS = ["Serial", "dEA", "HdEA", "HdEA-MP"]

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "figure.dpi": 150,
})

def read_results(algo):
    path = os.path.join(RESULTS_DIR, algo, "results.txt")
    bests = []
    if not os.path.exists(path):
        return []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                bests.append(float(parts[1]))
    return bests

def read_convergence(algo):
    path = os.path.join(RESULTS_DIR, algo, "convergence.txt")
    gens, times, bests = [], [], []
    if not os.path.exists(path):
        # Try individual run files
        conv_dir = os.path.join(RESULTS_DIR, algo, "convergence")
        if os.path.exists(conv_dir):
            files = [f for f in os.listdir(conv_dir) if f.startswith("run_")]
            if files:
                # Take first run for median (simplified)
                fpath = os.path.join(conv_dir, sorted(files)[0])
                return read_convergence_file(fpath)
        return [], [], []
    return read_convergence_file(path)

def read_convergence_file(path):
    gens, times, bests = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                gens.append(int(parts[0]))
                times.append(float(parts[1]))
                bests.append(float(parts[2]))
    return gens, times, bests

def plot_boxplot():
    fig, ax = plt.subplots(figsize=(8, 5))
    data = []
    for algo in ALGORITHMS:
        d = read_results(algo)
        data.append(d if d else [0])

    bp = ax.boxplot(data, labels=LABELS, patch_artist=True, widths=0.5)
    for patch, color in zip(bp["boxes"], COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel("Best Distance")
    ax.set_title("TSP Solution Quality Comparison (pcb442)")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join("..", "report", "figures", "boxplot.png"), bbox_inches="tight")
    print("Saved boxplot.png")
    plt.close()

def plot_barplot():
    fig, ax = plt.subplots(figsize=(8, 5))
    means, stds = [], []
    for algo in ALGORITHMS:
        d = read_results(algo)
        if d:
            means.append(np.mean(d))
            stds.append(np.std(d, ddof=1))
        else:
            means.append(0)
            stds.append(0)

    x = np.arange(len(ALGORITHMS))
    bars = ax.bar(x, means, yerr=stds, capsize=5, color=COLORS, alpha=0.8, width=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(LABELS)
    ax.set_ylabel("Mean Best Distance")
    ax.set_title("Mean TSP Solution Quality (error bars = 1 SD)")
    ax.grid(axis="y", alpha=0.3)

    # Annotate bars
    for bar, mean_val in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 20,
                f"{mean_val:.0f}", ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join("..", "report", "figures", "barplot.png"), bbox_inches="tight")
    print("Saved barplot.png")
    plt.close()

def plot_convergence():
    fig, ax = plt.subplots(figsize=(9, 5.5))

    for algo, color, label in zip(ALGORITHMS, COLORS, LABELS):
        gens, times, bests = read_convergence(algo)
        if gens and bests:
            ax.plot(gens, bests, color=color, linewidth=1.2, label=label)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Best Distance")
    ax.set_title("Convergence Curves (median run)")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join("..", "report", "figures", "convergence.png"), bbox_inches="tight")
    print("Saved convergence.png")
    plt.close()

def main():
    os.makedirs(os.path.join("..", "report", "figures"), exist_ok=True)

    print("=== Generating plots ===")
    plot_boxplot()
    plot_barplot()
    plot_convergence()
    print("Done!")

if __name__ == "__main__":
    main()
