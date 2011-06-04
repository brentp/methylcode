tbl = read.delim('embryo.endosperm.counts.txt', header=TRUE, row.names=1, stringsAsFactors=TRUE)
library(DESeq)

cds = newCountDataSet(tbl,  c("embryo", "endosperm"))
cds = estimateSizeFactors(cds)
# NOTE, we dont have replication in this case.
cds = estimateVarianceFunctions(cds, method="blind")

res = nbinomTest(cds, "embryo", "endosperm")
# p-adjusted < 0.05
resvalid = res[!is.na(res$padj),]
resSig = resvalid[ resvalid$padj < .05, ]
write.table(resSig, "embryo.endosperm.sig.genes.txt", row.names=F, sep="\t",
                    quote=F)

