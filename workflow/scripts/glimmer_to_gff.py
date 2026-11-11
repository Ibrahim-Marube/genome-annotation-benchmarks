#!/usr/bin/env python3
"""Convert a Glimmer3 .predict file to GFF3 (CDS lines, stop codon included).

Glimmer versions differ on whether the reported end coordinate includes the
stop codon, so this script checks the genome sequence itself: it counts how
many predictions have a stop codon (TAA/TAG/TGA) at the reported end versus
3 bases beyond it, picks the right convention, and prints what it found.

Usage: glimmer_to_gff.py genome.fna prefix.predict > out.gff
"""
import sys

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def read_predict(path):
    genes, seqid = [], None
    for line in open(path):
        if line.startswith(">"):
            seqid = line[1:].split()[0]
            continue
        f = line.split()
        if seqid is None or len(f) < 4 or not f[0].startswith("orf"):
            continue
        try:
            start, end = int(f[1]), int(f[2])
        except ValueError:
            continue
        strand = "-" if f[3].startswith("-") else "+"
        lo, hi = min(start, end), max(start, end)
        genes.append((seqid, f[0], strand, lo, hi))
    return genes


def codon(seqs, seqid, strand, lo, hi, offset):
    """Codon at the stop end; offset 0 = inside reported span, 3 = just beyond it."""
    s = seqs[seqid]
    if strand == "+":
        a = hi - 3 + offset
        c = s[a:a + 3]
    else:
        a = lo - 1 - offset
        c = s[a:a + 3].translate(COMP)[::-1] if a >= 0 else ""
    return c


def main():
    seqs = read_fasta(sys.argv[1])
    genes = read_predict(sys.argv[2])
    inside = sum(codon(seqs, g[0], g[2], g[3], g[4], 0) in STOPS for g in genes)
    beyond = sum(codon(seqs, g[0], g[2], g[3], g[4], 3) in STOPS for g in genes)
    n = len(genes)
    print(f"{n} predictions; stop codon at reported end: {inside}, "
          f"3 bases beyond: {beyond}", file=sys.stderr)
    if inside >= beyond and inside > 0.5 * n:
        extend = False
        print("-> coordinates already include the stop codon", file=sys.stderr)
    elif beyond > 0.5 * n:
        extend = True
        print("-> coordinates exclude the stop codon; extending by 3", file=sys.stderr)
    else:
        sys.exit("Cannot determine stop-codon convention; check input files.")
    print("##gff-version 3")
    for seqid, gid, strand, lo, hi in genes:
        if extend:
            if strand == "+":
                hi += 3
            else:
                lo -= 3
        print(f"{seqid}\tGlimmer3\tCDS\t{lo}\t{hi}\t.\t{strand}\t0\tID={gid}")


if __name__ == "__main__":
    main()
