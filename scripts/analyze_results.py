#!/usr/bin/env python3
"""
Analyze DiffDock_ZY docking results: parse confidence scores, classify levels,
rank predictions within and across complexes, optionally export CSV.

Works on:
  - a single result dir (containing raw_response.json / confidence_scores.txt)
  - a batch dir (containing many per-complex subdirs)

Usage:
    python scripts/analyze_results.py results/single/
    python scripts/analyze_results.py results/batch/ --export ranking.csv --top 5
"""
import argparse
import csv
import json
import os
import sys


def level(c):
    if c > 0:
        return "High"
    if c > -1.5:
        return "Moderate"
    return "Low"


def load_confidences(result_dir):
    """Return list of confidences from a single result dir, or None if not found."""
    raw = os.path.join(result_dir, "raw_response.json")
    if os.path.isfile(raw):
        with open(raw) as f:
            data = json.load(f)
        return data.get("position_confidence") or []
    # fall back to confidence_scores.txt
    txt = os.path.join(result_dir, "confidence_scores.txt")
    if os.path.isfile(txt):
        vals = []
        with open(txt) as f:
            next(f, None)  # header
            for line in f:
                parts = line.split("\t")
                if len(parts) >= 3:
                    try:
                        vals.append(float(parts[2]))
                    except ValueError:
                        pass
        return vals
    return None


def find_result_dirs(root):
    """Yield (name, dir) for each result dir under root (self or subdirs)."""
    if load_confidences(root) is not None:
        yield (os.path.basename(root.rstrip("/")) or root, root)
        return
    for entry in sorted(os.listdir(root)):
        sub = os.path.join(root, entry)
        if os.path.isdir(sub) and load_confidences(sub) is not None:
            yield (entry, sub)


def main():
    ap = argparse.ArgumentParser(description="Analyze DiffDock_ZY results")
    ap.add_argument("results_dir", help="Result dir (single) or batch dir")
    ap.add_argument("--top", type=int, default=5, help="Show top-N poses per complex")
    ap.add_argument("--threshold", type=float, help="Only show complexes with best conf >= threshold")
    ap.add_argument("--export", help="Export ranking to this CSV path")
    args = ap.parse_args()

    if not os.path.isdir(args.results_dir):
        sys.exit(f"ERROR: not a directory: {args.results_dir}")

    dirs = list(find_result_dirs(args.results_dir))
    if not dirs:
        sys.exit(f"ERROR: no docking results found under {args.results_dir}")

    rows = []
    for name, d in dirs:
        conf = load_confidences(d) or []
        if not conf:
            continue
        best = max(conf)
        if args.threshold is not None and best < args.threshold:
            continue
        rows.append((name, best, conf))

    rows.sort(key=lambda x: x[1], reverse=True)

    print("=" * 64)
    print(f"DiffDock_ZY 结果分析  ({len(rows)} 个化合物)")
    print("=" * 64)
    for rank, (name, best, conf) in enumerate(rows, 1):
        print(f"\n#{rank}  {name}   best={best:.3f} ({level(best)})")
        ranked = sorted(conf, reverse=True)[: args.top]
        for i, c in enumerate(ranked, 1):
            print(f"     pose {i}: {c:+.3f}  ({level(c)})")

    # distribution summary
    highs = sum(1 for _, b, _ in rows if b > 0)
    mods = sum(1 for _, b, _ in rows if -1.5 < b <= 0)
    lows = sum(1 for _, b, _ in rows if b <= -1.5)
    print("\n" + "-" * 64)
    print(f"分布(按最优位姿): High={highs}  Moderate={mods}  Low={lows}")

    if args.export:
        with open(args.export, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["rank", "complex_name", "best_confidence", "level"])
            for rank, (name, best, _) in enumerate(rows, 1):
                w.writerow([rank, name, f"{best:.4f}", level(best)])
        print(f"\n[ok] 排名已导出: {args.export}")


if __name__ == "__main__":
    main()
