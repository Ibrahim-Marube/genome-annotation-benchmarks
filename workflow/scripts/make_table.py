#!/usr/bin/env python3
"""Build a markdown comparison table from compare_gff.py summary files.

Usage: make_table.py results/<tool>/<genome>.summary.txt ...
"""
import re
import sys
from pathlib import Path


def parse(path):
    text = Path(path).read_text()
    pct = lambda x: f"{100 * float(x):.1f}%"
    found = re.search(r"gene found\)\s*:\s*(\d+)\s+sensitivity ([\d.]+)\s+precision ([\d.]+)", text)
    exact = re.search(r"exact matches.*?:\s*(\d+)\s+sensitivity ([\d.]+)\s+precision ([\d.]+)", text)
    return {
        "predicted": int(re.search(r"predicted genes\s*:\s*(\d+)", text).group(1)),
        "found": f"{pct(found.group(2))} / {pct(found.group(3))}",
        "exact": f"{pct(exact.group(2))} / {pct(exact.group(3))}",
        "start_diff": int(re.search(r"start disagrees.*?:\s*(\d+)", text).group(1)),
    }


def main():
    rows = []
    for p in sys.argv[1:]:
        tool, name = Path(p).parent.name, Path(p).name.replace(".summary.txt", "")
        rows.append((name, tool, parse(p)))
    print("| Genome | Tool | Predicted | Found (sens / prec) | Exact (sens / prec) | Start disagreements |")
    print("|---|---|---|---|---|---|")
    for name, tool, r in sorted(rows):
        print(f"| {name} | {tool} | {r['predicted']:,} | {r['found']} | {r['exact']} | {r['start_diff']} |")


if __name__ == "__main__":
    main()
