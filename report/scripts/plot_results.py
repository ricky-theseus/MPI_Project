"""plot_results.py — 读 results_*.txt → 箱线图 + 柱状图 → report/figures/"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
FIGURES = os.path.join(os.path.dirname(BASE), "figures")
COURSE = os.path.join(os.path.dirname(os.path.dirname(BASE)), "course-project")

FILES = {
    "dEA (Alg1)":    os.path.join(COURSE, "results_dEA.txt"),
    "HdEA (Alg2)":   os.path.join(COURSE, "results_HdEA.txt"),
    "HdEA-MP (Alg3)": os.path.join(COURSE, "results_HdEA_MP.txt"),
}

def read_values(path):
    vals = []
    if not os.path.exists(path):
        return vals
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    vals.append(float(parts[1]))
                except ValueError:
                    continue
    return vals

def main():
    data = {}
    for name, path in FILES.items():
        v = read_values(path)
        if v:
            data[name] = np.array(v)
            print(f"{name}: {len(v)} values, mean={np.mean(v):.0f}")
        else:
            print(f"{name}: NO DATA!")

    if not data:
        print("No data to plot!")
        return

    names = list(data.keys())
    values = [data[n] for n in names]
    colors = ["#4C72B0", "#DD8452", "#55A868"]

    # --- Boxplot ---
    fig, ax = plt.subplots(figsize=(8, 5))
    bp = ax.boxplot(values, patch_artist=True, widths=0.5)
    ax.set_xticklabels(names)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    for i, v in enumerate(values):
        jitter = np.random.normal(0, 0.04, size=len(v))
        ax.scatter(np.full_like(v, i + 1) + jitter, v, color="black", alpha=0.6, s=40, zorder=3)
    ax.set_ylabel("Best Tour Distance")
    ax.set_title("TSP (pcb442) — 10-Run Comparison")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES, "boxplot.png"), dpi=200)
    print(f"Saved boxplot.png")

    # --- Bar chart ---
    fig, ax = plt.subplots(figsize=(8, 5))
    means = [np.mean(v) for v in values]
    stds  = [np.std(v, ddof=1) for v in values]
    x = np.arange(len(names))
    bars = ax.bar(x, means, yerr=stds, capsize=8, color=colors, alpha=0.7, width=0.5)
    for bar, m in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 200,
                f"{m:.0f}", ha="center", va="bottom", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Mean Best Distance")
    ax.set_title("TSP (pcb442) — Mean ± Std (10 runs)")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES, "barplot.png"), dpi=200)
    print(f"Saved barplot.png")

if __name__ == "__main__":
    main()
