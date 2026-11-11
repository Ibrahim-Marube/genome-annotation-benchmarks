# genome-annotation-benchmarks

Benchmarking ab initio gene prediction against curated reference annotations,
run on the University of Potsdam HPC cluster. Built while learning genome
annotation and sequence bioinformatics (MSc Bioinformatics, Uni Potsdam).

## Project 1: Prodigal vs RefSeq on *E. coli* K-12 MG1655

**Question:** how accurate is Prodigal 2.6.3 on a well-annotated bacterial
genome, and where does it go wrong?

**Data:** NCBI RefSeq GCF_000005845.2 (NC_000913.3, 4.6 Mb). The sequence is
the input and the RefSeq GFF3 is the answer key.

**Method:** predict genes with `prodigal -p single`; match genes to the
reference by strand and stop coordinate. A gene is *found* if the stop
matches and *exact* if the start matches too. Pseudogenes are excluded from
the reference gene set. Errors are then classified by `analyze_errors.py`.

### Results

| | Found (stop matches) | Exact (start + stop) |
|---|---|---|
| Sensitivity | 96.0% (4,119 / 4,290) | 90.5% (3,883) |
| Precision | 95.4% (4,119 / 4,319) | 89.9% |

**Extra predictions (200):** 137 overlap a RefSeq pseudogene, 45 are
intergenic, 12 are antisense to a real gene, 6 have a different stop.
Counting pseudogene overlaps as defensible, precision is about 98.5%.

**Missed genes (171):** 157 are shorter than 300 nt; 138 have no prediction
at all and 33 are overlapped by a prediction with a different stop.

**Start codons:** 236 genes have the right stop but a different start;
Prodigal's start is downstream in 141 cases and upstream in 95.

### Limitations
- RefSeq start sites are not all experimentally verified, so "exact match" is
  agreement with the reference, not proof of correctness.
- One genome, one tool, default settings. Next: Bakta/Prokka, and plastid and
  eukaryotic genomes.

### Reproduce (on the cluster)
~~~bash
module load bio/prodigal/2.6.3-GCCcore-10.2.0
# download genome + GFF from NCBI (see data/ecoli), then:
prodigal -i genome.fna -p single -f gff -o results/prodigal/ecoli.gff
python3 workflow/scripts/compare_gff.py reference.gff results/prodigal/ecoli.gff results/prodigal/ecoli
python3 workflow/scripts/analyze_errors.py reference.gff results/prodigal/ecoli.gff
~~~
