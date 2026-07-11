# analyze_results.py
# Computes mean, std, min, max and paired t-test across algorithms

import os
import math

RESULTS_DIR = "results"
ALGORITHMS = ["serial", "dEA", "HdEA", "HdEA_MP"]

def read_results(algo):
    path = os.path.join(RESULTS_DIR, algo, "results.txt")
    runs, bests, times = [], [], []
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found")
        return [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                runs.append(int(parts[0]))
                bests.append(float(parts[1]))
                times.append(float(parts[2]))
    return runs, bests, times

def mean(data):
    return sum(data) / len(data) if data else 0

def stdev(data):
    if len(data) < 2:
        return 0
    m = mean(data)
    var = sum((x - m) ** 2 for x in data) / (len(data) - 1)
    return math.sqrt(var)

def paired_ttest(a, b):
    """Paired two-sided t-test. Returns t-statistic and p-value."""
    n = min(len(a), len(b))
    if n < 2:
        return 0, 1
    diffs = [a[i] - b[i] for i in range(n)]
    dbar = mean(diffs)
    s = stdev(diffs)
    if s == 0:
        return 0, 1.0
    t = dbar / (s / math.sqrt(n))
    # Approximate p-value using t-distribution
    df = n - 1
    p = 2 * (1 - t_cdf(abs(t), df))
    return t, p

def t_cdf(x, df):
    """Approximate Student's t CDF using regularized incomplete beta function."""
    if x <= 0:
        return 0.5
    # Using approximation for large df or simple case
    # Abramowitz and Stegun approximation
    a = df / 2.0
    b = 0.5
    z = df / (df + x * x)
    return 1 - 0.5 * betai(a, b, z)

def betai(a, b, x):
    """Incomplete beta function using continued fraction approximation."""
    if x < 0 or x > 1:
        return 0
    if x == 0 or x == 1:
        return 0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    return math.exp(a * math.log(x) + b * math.log(1 - x) - lbeta) / a * cf_cont(a, b, x)

def cf_cont(a, b, x):
    """Continued fraction for incomplete beta."""
    MAX_ITER = 200
    EPS = 3e-12
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < EPS:
        d = EPS
    d = 1.0 / d
    h = d
    for m in range(1, MAX_ITER + 1):
        m2 = 2 * m
        # Even step
        numer = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + numer * d
        if abs(d) < EPS:
            d = EPS
        c = 1.0 + numer / c
        if abs(c) < EPS:
            c = EPS
        d = 1.0 / d
        h *= d * c
        # Odd step
        numer = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + numer * d
        if abs(d) < EPS:
            d = EPS
        c = 1.0 + numer / c
        if abs(c) < EPS:
            c = EPS
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break
    return h

def print_table():
    print()
    print(f"{'Algorithm':<12} {'Mean':<10} {'Std':<10} {'Min':<10} {'Max':<10}")
    print("-" * 52)
    data = {}
    for algo in ALGORITHMS:
        runs, bests, times = read_results(algo)
        data[algo] = bests
        if bests:
            m = mean(bests)
            s = stdev(bests)
            mn = min(bests)
            mx = max(bests)
            print(f"{algo:<12} {m:<10.1f} {s:<10.1f} {mn:<10.1f} {mx:<10.1f}")
        else:
            print(f"{algo:<12} {'N/A':<10} {'N/A':<10} {'N/A':<10} {'N/A':<10}")

    print(f"\n{'Pair':<25} {'t-stat':<10} {'p-value':<10} {'Significant?':<12}")
    print("-" * 57)
    pairs = [("serial", "dEA"), ("dEA", "HdEA"), ("HdEA", "HdEA_MP"),
             ("serial", "HdEA"), ("serial", "HdEA_MP")]
    for a1, a2 in pairs:
        if a1 in data and a2 in data and data[a1] and data[a2]:
            t, p = paired_ttest(data[a1], data[a2])
            sig = "YES" if p < 0.05 else "no"
            print(f"{a1} vs {a2:<15} {t:<10.4f} {p:<10.6f} {sig:<12}")

def main():
    print("=== Experiment Results Analysis ===")
    print_table()
    print()

if __name__ == "__main__":
    main()
