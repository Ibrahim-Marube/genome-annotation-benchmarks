#!/usr/bin/env python3
"""Classify why a predicted gene set disagrees with a reference annotation.

Usage: analyze_errors.py reference.gff predicted.gff
Run from anywhere; compare_gff.py must sit in the same directory.
"""
import sys
from collections import Counter

from compare_gff import parse_cds


def span(gene):
    (seqid, strand, stop), (start, _) = gene
    return seqid, strand, min(start, stop), max(start, stop)


def pseudogenes(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9 and f[2] == "pseudogene":
            out.append((f[0], f[6], int(f[3]), int(f[4])))
    return out


def overlaps(a, b):
    return a[0] == b[0] and a[2] <= b[3] and b[2] <= a[3]


def classify_extra(g, ref_spans, pseudo):
    if any(overlaps(g, p) and g[1] == p[1] for p in pseudo):
        return "overlaps reference pseudogene"
    if any(overlaps(g, r) and g[1] == r[1] for r in ref_spans):
        return "overlaps reference gene, same strand (different stop)"
    if any(overlaps(g, r) and g[1] != r[1] for r in ref_spans):
        return "overlaps reference gene on opposite strand"
    return "intergenic (no reference feature)"


def classify_missed(g, pred_spans):
    if any(overlaps(g, p) and g[1] == p[1] for p in pred_spans):
        return "overlapped by a prediction, same strand (different stop)"
    return "no prediction at all"


def length_bin(n):
    return "<300 nt" if n < 300 else "300-599 nt" if n < 600 else ">=600 nt"


def main():
    ref_path, pred_path = sys.argv[1], sys.argv[2]
    ref, pred = parse_cds(ref_path, ref=True), parse_cds(pred_path)
    ref_spans = [span(g) for g in ref.items()]
    pred_spans = [span(g) for g in pred.items()]
    pseudo = pseudogenes(ref_path)

    print("== EXTRA predictions (not in reference) ==")
    extra = [span((k, pred[k])) for k in pred.keys() - ref.keys()]
    print(f"total {len(extra)}")
    for label, n in Counter(classify_extra(g, ref_spans, pseudo) for g in extra).most_common():
        print(f"  {n:5d}  {label}")
    for label, n in sorted(Counter(length_bin(g[3] - g[2] + 1) for g in extra).items()):
        print(f"  length {label}: {n}")

    print("\n== MISSED reference genes ==")
    missed = [span((k, ref[k])) for k in ref.keys() - pred.keys()]
    print(f"total {len(missed)}")
    for label, n in Counter(classify_missed(g, pred_spans) for g in missed).most_common():
        print(f"  {n:5d}  {label}")
    for label, n in sorted(Counter(length_bin(g[3] - g[2] + 1) for g in missed).items()):
        print(f"  length {label}: {n}")

    print("\n== START-CODON disagreements (same stop, different start) ==")
    diff = [k for k in ref.keys() & pred.keys() if ref[k][0] != pred[k][0]]
    shorter = sum(1 for k in diff if abs(pred[k][0] - k[2]) < abs(ref[k][0] - k[2]))
    print(f"total {len(diff)}: Prodigal start downstream (shorter protein) {shorter}, "
          f"upstream (longer protein) {len(diff) - shorter}")


if __name__ == "__main__":
    main()
