#!/usr/bin/Rscript
#
# Analysing featurecounts produced transcript counts for DGE.
#
#
# Expects two files:
#
#  samples.txt, sample information,
#               col1 which is used to estbalish groups
#               col3 will be used as colnames
#  e.g.
#
#    control rep1    SRR4048970
#    control rep2    SRR4048971
#    control rep3    SRR4048972
#    15d-PGJ2        rep1    SRR4048973
#    15d-PGJ2        rep2    SRR4048974
#    15d-PGJ2        rep3    SRR4048975
#
#  featurecounts produced count file as command line argument 1
#
# featureCounts can be used like this:
# e.g.
# featureCounts -a Mus_musculus.GRCm38.85.gtf -o featurecounts.genes.counts star/*.bam
# featureCounts -a Mus_musculus.GRCm38.85.gtf -o featurecounts.tx.counts -g transcript_id star/*.bam
#

# length in bases
counts2tpm <- function(counts, length)
    {
        c.rpk <- (counts/(length/1000)) # reads per kilobase
        c.rpksums <- colSums(c.rpk)/1e6 # scaling factors
        t(t(c.rpk)/c.rpksums) # scale rpks
    }

time <- format(Sys.time(), "%Y%m%d-%H%M%S")
args <- commandArgs(trailingOnly = TRUE)
if ( length(args)<1 ) {
   stop("USAGE: featureCounts_edgeR_DGE.R featurecount-file")
}

file=args[1]

fc <-  read.table(file.path('.', file), header = TRUE, skip=1)

# samples.txt
samples <- read.table(file.path('.', "samples.txt"), header = FALSE)
# some information about the groupings of samples
group <- samples[,1]
group

counts <- fc[,-(1:6)]
rownames(counts) <- fc$Geneid
colnames(counts)
tpm <- counts2tpm(counts, fc$Length)
colnames(tpm) <- paste(colnames(tpm),'TPM',sep='_')

library(edgeR)
library(methods)

# tximport type of normalising with length
# # https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# get matrix of lengths
normMat = replicate(6, fc$Length)
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(counts/normMat)) + log(colSums(counts/normMat))
d <- DGEList(counts, group=group)
d$offset <- t(t(log(normMat)) + o)
# d is now ready for estimate dispersion functions see edgeR User's Guide

d$tpm <- tpm

# create desgin matrix
#design <- model.matrix(~group)

cat("\nBEFORE filtering stats:\n")
cat("colsums:\n")
colSums(counts)
summary(colSums(counts))
cat("\nrowsums:\n")
summary(rowSums(counts))

#--------------------------------------------------------------------
# INDEPENDENT FILTERING
#--------------------------------------------------------------------
# Independent filtering to get rid of lowly expressed transcripts
# that have no chance of being DE in a stat test in the first place.
# Thus, the chance for the remainder is higher to be called DE

# Based on Nature Protocol: http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html
# based on CPM values and replicate numbers
# CRITICAL STEP In edgeR, it is recommended to remove features without
# at least 1 read per million in n of the samples, where n is the size of the smallest group of replicates
# WE ARE USING TPM
use = rowSums(tpm >1) >= min(table(d$samples$group))  # num smallest group size reps at least > 1 tpm

# apply filter
d <- d[use,]
# also to tpm table
d$tpm <- d$tpm[use,]

# reset libsizes
# this will change the TMM values for the genes as oposed to the original complete set.
# Use original CPM (from d) for non-DE related tasks.
# Not sure if strictkly speaking neceassry to adjust lib-sizes.
#d$samples$lib.size <- colSums(d$counts)

cat("\nAFTER filtering stats:\n")
cat("colsums::\n")
colSums(d$counts)
summary(colSums(d$counts))
cat("\nrowsums:\n")
summary(rowSums(d$counts))
#--------------------------------------------------------------------

# edgeR DE
# we DO not normalise again with reduced table.
# we use offsets calcualted above
#d <- calcNormFactors(d, method="TMM")
#d = estimateCommonDisp(d)
#d = estimateTagwiseDisp(d)
d = estimateDisp(d)

de.com = exactTest(d)

# print some stats
cat("\np.value<0.05:\n")
summary(decideTestsDGE(de.com, p.value=0.05))
cat("\np.value<0.01:\n")
summary(decideTestsDGE(de.com, p.value=0.01))
topTags(de.com)

#some plotting
#x = sum(p.adjust(de.com$table$PValue,method="BH") < 0.05)
#de.tags <- rownames(topTags(de.com, n=x)$table)
#plotSmear(d, de.tags=de.tags)

dge <- de.com$table
dge$FDR <- p.adjust(method="fdr",p=dge$PValue)
datasub <- cbind(dge, d$counts, d$tpm)

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(datasub), datasub),
            file=paste(time, "edgeR_DE-CONTROL-VS-TREAT.txt", sep="_"),
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

