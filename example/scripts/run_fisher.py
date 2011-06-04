from fisher import pvalue
import sys

ct_counts = open(sys.argv[1])
ngenes = sum(1 for _ in ct_counts)
ct_counts.seek(0)
header = ct_counts.readline()

cutoff = 0.0001 / ngenes
for gene, ac, at, bc, bt in (line.rstrip().split() for line in ct_counts):
    p = pvalue(*map(int, (ac, at, bc, bt)))
    if p.two_tail > cutoff: continue
    print gene, ac, at, bc, bt, p.two_tail

