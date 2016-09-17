#!/usr/bin/Rscript
#
# Producing TMM-TPM normlised counts for featureCounts result table.
#
# WRITES TO STDOUT.
#
# Expects one files:
#
#  featurecounts produced count file -> here expects name featurecounts.genes.counts.gz
# 

library(edgeR)
library(methods)

fc <-  read.table(file.path('.', "featurecounts.genes.counts.gz"), header = TRUE, skip=1)

counts <- fc[,-(2:6)]
rownames(counts) <- counts[,1]
counts[,1] <- NULL

# edgeR
d <- DGEList(counts=counts)
d <- calcNormFactors(d, method="TMM")

# get table with TMM normalised TPM for all genes
cpms <- cpm(d)
colnames(cpms) <- paste(colnames(cpms),'TMMTPM',sep='_')
ccpms <- cbind(d$counts, cpms)
write.table(data.frame("Genes"=rownames(ccpms), ccpms),
            file="",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
