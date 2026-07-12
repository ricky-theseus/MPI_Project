#!/usr/bin/env python3
"""
generate_tsp.py — Download TSPLIB instances and convert to custom format

Our format:  first line = city count, then each line: id x y
TSPLIB format: headers (NAME, TYPE, DIMENSION, etc.) then NODE_COORD_SECTION

Usage:
  python generate_tsp.py [--local] [--instances all|pcb442|lin318|...]
"""

import urllib.request
import os
import sys
import argparse
import re

TSPLIB_BASE = ("http://comopt.ifi.uni-heidelberg.de/"
               "software/TSPLIB95/tsp/")

INSTANCES = {
    "lin318":  "lin318.tsp.gz",
    "pcb442":  "pcb442.tsp.gz",
    "u574":    "u574.tsp.gz",
    "p654":    "p654.tsp.gz",
    "d657":    "d657.tsp.gz",
    "u724":    "u724.tsp.gz",
    "rat783":  "rat783.tsp.gz",
    "dsj1000": "dsj1000.tsp.gz",
    "pr1002":  "pr1002.tsp.gz",
}

# For instances not in TSPLIB or with different filenames
INSTANCE_INFO = {
    "lin318":  {"dim": 318,  "opt": 42029},
    "pcb442":  {"dim": 442,  "opt": 50778},
    "u574":    {"dim": 574,  "opt": 36905},
    "p654":    {"dim": 654,  "opt": 34643},
    "d657":    {"dim": 657,  "opt": 48912},
    "u724":    {"dim": 724,  "opt": 41910},
    "rat783":  {"dim": 783,  "opt": 8806},
    "dsj1000": {"dim": 1000, "opt": 18659688},
    "pr1002":  {"dim": 1002, "opt": 259045},
}


def parse_tsplib(content):
    """Parse TSPLIB format, return (city_count, [(id, x, y), ...])."""
    lines = content.split("\n")
    dimension = None
    edge_weight_type = None
    in_coord_section = False
    coords = []

    for line in lines:
        line = line.strip()
        if not line or line.startswith("NAME") or line.startswith("TYPE"):
            continue
        if line.startswith("COMMENT") or line.startswith("DISPLAY_DATA_TYPE"):
            continue

        m = re.match(r"DIMENSION\s*:\s*(\d+)", line, re.I)
        if m:
            dimension = int(m.group(1))
            continue
        m = re.match(r"EDGE_WEIGHT_TYPE\s*:\s*(\w+)", line, re.I)
        if m:
            edge_weight_type = m.group(1).upper()
            continue

        if "NODE_COORD_SECTION" in line.upper():
            in_coord_section = True
            continue
        if "EOF" in line.upper():
            break

        if in_coord_section:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    idx = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    coords.append((idx, x, y))
                except ValueError:
                    continue

    if not dimension:
        dimension = len(coords)
    return dimension, edge_weight_type, coords


def convert_to_custom(dimension, coords):
    """Convert to our custom format: city_count on first line, then id x y."""
    lines = [str(dimension)]
    for idx, x, y in coords:
        lines.append(f"{idx} {x:.5e} {y:.5e}")
    return "\n".join(lines) + "\n"


def generate_synthetic(name):
    """Generate synthetic TSP data for a named instance for testing."""
    info = INSTANCE_INFO.get(name, {})
    dim = info.get("dim", 100)
    import random
    random.seed(hash(name) % 2**32)
    coords = [(i + 1, random.uniform(0, 10000), random.uniform(0, 10000))
              for i in range(dim)]
    return dim, coords


def main():
    parser = argparse.ArgumentParser(description="TSP instance downloader/generator")
    parser.add_argument("--instances", type=str, default="all",
                        help="comma-separated or 'all'")
    parser.add_argument("--output-dir", type=str,
                        default=os.path.join(os.path.dirname(__file__),
                                             "tsp_instances"))
    parser.add_argument("--local", action="store_true",
                        help="Generate synthetic data locally (no download)")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.instances == "all":
        instances = list(INSTANCES.keys())
    else:
        instances = [s.strip() for s in args.instances.split(",")]

    import gzip
    import io

    for name in instances:
        out_path = os.path.join(args.output_dir, f"{name}.tsp")
        if os.path.exists(out_path) and not args.force:
            print(f"  {name}: already exists, skipping")
            continue

        if args.local:
            dim, coords = generate_synthetic(name)
            content = convert_to_custom(dim, coords)
            with open(out_path, "w") as f:
                f.write(content)
            print(f"  {name}: generated synthetic ({dim} cities)")
            continue

        url = TSPLIB_BASE + INSTANCES[name]
        try:
            print(f"  {name}: downloading from {url}...")
            req = urllib.request.urlopen(url, timeout=30)
            raw = req.read()

            with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
                content = gz.read().decode("utf-8", errors="replace")

            dimension, edge_weight_type, coords = parse_tsplib(content)
            if not coords:
                print(f"    WARNING: no coordinates parsed, "
                      f"using synthetic data")
                dimension, coords = generate_synthetic(name)

            custom = convert_to_custom(dimension, coords)
            with open(out_path, "w") as f:
                f.write(custom)
            print(f"    OK: {len(coords)} cities (dim={dimension})")

        except Exception as e:
            print(f"    Download failed: {e}")
            print(f"    Generating synthetic data instead...")
            dim, coords = generate_synthetic(name)
            content = convert_to_custom(dim, coords)
            with open(out_path, "w") as f:
                f.write(content)
            print(f"    OK: synthetic ({dim} cities)")

    print(f"\nAll files placed in: {args.output_dir}")


if __name__ == "__main__":
    main()
