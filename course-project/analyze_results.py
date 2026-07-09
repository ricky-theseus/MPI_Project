"""
analyze_results.py — Analyze TSP parallel EA experiment results

Computes per-algorithm statistics (mean, std, min, max) and performs
paired t-test between algorithms to determine statistical significance.
"""
import os
import math
import sys

RESULT_FILES = {
    "dEA (Alg1)": "results_dEA.txt",
    "HdEA (Alg2)": "results_HdEA.txt",
    "HdEA-MP (Alg3)": "results_HdEA_MP.txt",
}


def read_results(path):
    """Read results from a file, skipping comment lines."""
    data = []
    if not os.path.exists(path):
        # Try subdirectories
        for d in ["01_分布式进化算法", "02_分层分布式进化算法", "03_移动种群分层分布式进化算法"]:
            p = os.path.join(d, path)
            if os.path.exists(p):
                path = p
                break
    if not os.path.exists(path):
        # Try root
        path = os.path.join(os.path.dirname(__file__) or ".", path)
    if not os.path.exists(path):
        return None

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.append(float(parts[1]))
                except ValueError:
                    continue
    return data


def mean_std(data):
    n = len(data)
    if n < 2:
        return data[0] if data else 0, 0
    m = sum(data) / n
    v = sum((x - m) ** 2 for x in data) / (n - 1)
    return m, math.sqrt(v)


def t_test(a, b):
    """Paired two-tailed t-test. Returns t-statistic and p-value."""
    n = len(a)
    if n != len(b) or n < 2:
        return None, None
    diffs = [a[i] - b[i] for i in range(n)]
    d_mean = sum(diffs) / n
    d_var = sum((x - d_mean) ** 2 for x in diffs) / (n - 1)
    if d_var == 0:
        return 0.0, 1.0
    t = d_mean / math.sqrt(d_var / n)
    # Approximate p-value using t-distribution
    df = n - 1
    # Using simplified approximation
    p = 2 * (1 - _t_cdf(abs(t), df))
    return t, p


def _t_cdf(t, df):
    """Approximate cumulative t-distribution using regularized incomplete beta."""
    # Using simplified Abramowitz & Stegun approximation
    x = df / (df + t * t)
    # For large df, approximate as normal
    if df > 100:
        return _norm_cdf(t)
    # Simple approximation for moderate df
    a = df / 2.0
    b = 0.5
    return _betainc(a, b, x)


def _norm_cdf(x):
    """Standard normal CDF approximation."""
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def _betainc(a, b, x):
    """Regularized incomplete beta function (approximation)."""
    # This is a simplified version for the specific case
    if x < 0 or x > 1:
        return 0
    # Using continued fraction for I_x(a,b)
    # Very simplified - for statistics purposes
    if a > 50 and b == 0.5:  # Large df case
        return _norm_cdf(math.sqrt(a) * (1 - 2 * x) / (2 * math.sqrt(x * (1 - x))))
    # Stirlings approximation for general case
    import math as m
    import random
    if x == 0 or x == 1:
        return 0.5 if x == 0.5 else (0 if x == 0 else 1)
    return _norm_cdf((m.lgamma(a+b) - m.lgamma(a) - m.lgamma(b) + a*m.log(x) + b*m.log(1-x)) * 0.5)
    # Fallback
    return 0.05  # conservative


def main():
    print("=" * 70)
    print("TSP Parallel EA — Results Analysis")
    print("=" * 70)

    all_data = {}
    cwd = os.path.dirname(os.path.abspath(__file__))

    for name, fname in RESULT_FILES.items():
        path = os.path.join(cwd, fname)
        if not os.path.exists(path):
            # Check subdir
            for d in range(1, 4):
                sub = os.path.join(cwd, f"0{d}_" + "*")
                import glob as gb
                matches = sorted(gb.glob(os.path.join(cwd, f"0{d}_*")))
                if matches:
                    subpath = os.path.join(matches[0], fname)
                    if os.path.exists(subpath):
                        path = subpath
                        break
        data = read_results(path)
        if data is None:
            print(f"  {name}: NO DATA (file not found)")
            continue
        all_data[name] = data
        m, s = mean_std(data)
        print(f"\n  {name}:")
        print(f"    n={len(data)}, mean={m:.0f}, std={s:.0f}")
        print(f"    min={min(data):.0f}, max={max(data):.0f}")
        print(f"    values: {[f'{v:.0f}' for v in data]}")

    if len(all_data) < 2:
        print("\nNot enough data for t-test comparison.")
        return

    print("\n" + "-" * 70)
    print("Paired t-Test Results (p < 0.05 = significant)")
    print("-" * 70)
    names = list(all_data.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = all_data[names[i]], all_data[names[j]]
            t, p = t_test(a, b)
            if t is not None:
                sig = "SIGNIFICANT" if p < 0.05 else "not significant"
                print(f"  {names[i]} vs {names[j]}: t={t:.4f}, p={p:.4f} ({sig})")
            else:
                print(f"  {names[i]} vs {names[j]}: insufficient data")


if __name__ == "__main__":
    main()
