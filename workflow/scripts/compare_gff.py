#!/usr/bin/env python3
"""Compare predicted CDS (e.g. Prodigal GFF) against a reference NCBI GFF3.

A gene is identified by (seqid, strand, stop coordinate).
  * 3'-match : same stop  -> gene found, start may differ
  * exact    : same stop and same start
Usage: compare_gff.py reference.gff predicted.gff [out_prefix]
"""
import sys
from collections import defaultdict


def parse_cds(path, ref=False):
    """Return {(seqid, strand, stop): (start, end)} using 1-based inclusive coords."""
    groups = defaultdict(list)
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if ref and attrs.get("pseudo") == "true":
            continue
        key = attrs.get("ID", f"{f[0]}:{f[3]}-{f[4]}") if ref else f"{f[0]}:{f[3]}-{f[4]}:{f[6]}"
        groups[key].append((f[0], f[6], int(f[3]), int(f[4])))
    genes = {}
    for parts in groups.values():
        seqid, strand = parts[0][0], parts[0][1]
        lo = min(p[2] for p in parts)
        hi = max(p[3] for p in parts)
        stop = hi if strand == "+" else lo
        start = lo if strand == "+" else hi
        genes[(seqid, strand, stop)] = (start, stop)
    return genes


def main():
    ref_path, pred_path = sys.argv[1], sys.argv[2]
    prefix = sys.argv[3] if len(sys.argv) > 3 else "comparison"
    ref, pred = parse_cds(ref_path, ref=True), parse_cds(pred_path)
    both = ref.keys() & pred.keys()
    exact = [k for k in both if ref[k][0] == pred[k][0]]
    n_ref, n_pred = len(ref), len(pred)
    print(f"reference genes (non-pseudo CDS): {n_ref}")
    print(f"predicted genes                 : {n_pred}")
    print(f"3'-end matches (gene found)     : {len(both)}  "
          f"sensitivity {len(both)/n_ref:.3f}  precision {len(both)/n_pred:.3f}")
    print(f"exact matches (start+stop)      : {len(exact)}  "
          f"sensitivity {len(exact)/n_ref:.3f}  precision {len(exact)/n_pred:.3f}")
    print(f"start disagrees (same stop)     : {len(both) - len(exact)}")
    for name, keys, src in (("missed", ref.keys() - pred.keys(), ref),
                            ("extra", pred.keys() - ref.keys(), pred),
                            ("start_diff", set(both) - set(exact), None)):
        with open(f"{prefix}.{name}.tsv", "w") as out:
            out.write("seqid\tstrand\tstop\tref_start\tpred_start\n")
            for k in sorted(keys):
                out.write(f"{k[0]}\t{k[1]}\t{k[2]}\t"
                          f"{ref[k][0] if k in ref else ''}\t{pred[k][0] if k in pred else ''}\n")


if __name__ == "__main__":
    main()
