"""analyze_results.py — 读 results_*.txt → 统计量 + 配对 t 检验"""
import os, sys
import numpy as np
from scipy.stats import ttest_rel

BASE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(BASE)  # report/
COURSE = os.path.join(os.path.dirname(PROJECT), "course-project")

FILES = {
    "dEA":   os.path.join(COURSE, "results_dEA.txt"),
    "HdEA":  os.path.join(COURSE, "results_HdEA.txt"),
    "HdEA_MP": os.path.join(COURSE, "results_HdEA_MP.txt"),
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
        data[name] = v
        if v:
            print(f"{name:10s}  n={len(v):2d}  mean={np.mean(v):8.0f}  std={np.std(v,ddof=1):6.0f}  min={np.min(v):8.0f}  max={np.max(v):8.0f}")
            print(f"{'':10s}  values={[f'{x:.0f}' for x in v]}")
        else:
            print(f"{name:10s}  NO DATA")

    print("\n--- Paired t-test (p<0.05 significant) ---")
    names = list(data.keys())
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            a, b = np.array(data[names[i]]), np.array(data[names[j]])
            if len(a) != len(b) or len(a) < 2:
                print(f"{names[i]:10s} vs {names[j]:10s}  insufficient data")
                continue
            t, p = ttest_rel(a, b)
            sig = "SIGNIFICANT" if p < 0.05 else "not significant"
            d = np.mean(a) - np.mean(b)
            better = f"{names[i]} better" if d < 0 else (f"{names[j]} better" if d > 0 else "tie")
            print(f"{names[i]:10s} vs {names[j]:10s}  t={t:+.4f}  p={p:.4f} ({sig})  diff={d:+.0f} ({better})")

if __name__ == "__main__":
    main()
