#!/usr/bin/Rscript
#
# Analysing salmon produced transcript counts for DE genes.
# 
# Expects two files: 
#
#  tx_gene_map.txt, transcript to gene mapping e.g.
#    ENSMUST00000205326.1    ENSMUSG00000108652
#    ENSMUST00000206672.1    ENSMUSG00000108652
#    ENSMUST00000021676.11   ENSMUSG00000021252
#    ENSMUST00000124311.1    ENSMUSG00000021252
#    ...
#
#  samples.txt, sample information,
#               col1 which is used to estbalish groups 
#               col3 will be used as colnames
#               col4 which is use to load the salmon files
#  e.g.
#  
#    control rep1    SRR4048970      ./quants/SRR4048970/quant.sf
#    control rep2    SRR4048971      ./quants/SRR4048971/quant.sf
#    control rep3    SRR4048972      ./quants/SRR4048972/quant.sf
#    15d-PGJ2        rep1    SRR4048973      ./quants/SRR4048973/quant.sf
#    15d-PGJ2        rep2    SRR4048974      ./quants/SRR4048974/quant.sf
#    15d-PGJ2        rep3    SRR4048975      ./quants/SRR4048975/quant.sf
#

library(tximport)
library(edgeR)
library(methods)

# tx_gene_map.txt e.g.
t2g <- read.table(file.path('.', "tx_gene_map.txt"), header = FALSE)

# samples.txt, e.g.
samples <- read.table(file.path('.', "samples.txt"), header = FALSE)

# load salmon files specified in samples.txt col4
files <- file.path(samples[,4])
txi <- tximport(files = files, type="salmon", tx2gene = t2g )
# set colnames from samples.txt col3
colnames(txi$counts) <- samples[,3]
counts <- round(txi$counts) # round to integers

# some information about the groupings of samples from samples.txt col1
group <- samples[,1]
# create desgin matrix
design <- model.matrix(~group)

# edgeR
d <- DGEList(counts=counts,group=group)
d <- calcNormFactors(d, method="TMM")
d

# get table with TMM normalised TPM for all genes
cpms <- cpm(d)
colnames(cpms) <- paste(colnames(cpms),'TMMTPM',sep='_')
ccpms <- cbind(d$counts, cpms)
write.table(data.frame("Genes"=rownames(ccpms), ccpms),
            file="edgeR_COUNTS_TMM-TPM_samples.txt",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

#--- DE ANALYSIS ---------

cat("\nBEFORE filtering stats:\n")
cat("colsums:\n")
colSums(d$counts)
summary(colSums(d$counts))
cat("\nrowsums:\n")
summary(rowSums(d$counts))

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
use = rowSums(cpms >1) >= min(table(d$samples$group))  # num smallest group size reps at least > 1 cpm

# ALTERNATIVE FILTERING METHODS --------
# Based on rowsums alone
#use = (rowSums(d$counts)>10)

# Based on method in the DESeq viginette
# Here, we consider as a filter criterion rowsum rs, the overall
# sum of counts (irrespective of biological condition),
# and remove the genes in the lowest 40% quantile
#rs = rowSums(d$counts)
#q = quantile(rs, probs=0.4)
# if too many lowly expressed transcripts we can be more strict and
# filter on a higer tag count rowsum as these will also not be DE
#if (q<10) {
#    q=10
#}
#use = (rs > q)
#--------------------------------------

# apply filter
d1 <- d[use,]

# reset libsizes
# this will change the TMM values for the genes as oposed to the original complete set.
# Use original TPM (from d) for non-DE related tasks.
d1$samples$lib.size <- colSums(d1$counts)

cat("\nAFTER filtering stats:\n")
cat("colsums::\n")
colSums(d1$counts)
summary(colSums(d1$counts))
cat("\nrowsums:\n")
summary(rowSums(d1$counts))
#--------------------------------------------------------------------

# edgeR DE
# we normalise again with reduced table.
d1 <- calcNormFactors(d1, method="TMM")
d1 = estimateCommonDisp(d1)
d1 = estimateTagwiseDisp(d1)

de.com = exactTest(d1)

# print some stats
cat("\np.value<0.05:\n")
summary(decideTestsDGE(de.com, p.value=0.05))
cat("\np.value<0.01:\n")
summary(decideTestsDGE(de.com, p.value=0.01))
topTags(de.com)

#some plotting
#x = sum(p.adjust(de.com$table$PValue,method="BH") < 0.05)
#de.tags <- rownames(topTags(de.com, n=x)$table)
#plotSmear(d1, de.tags=de.tags)

datasub <- cbind(d1$counts, de.com$table)
datasub$FDR <- p.adjust(method="fdr",p=datasub$PValue)

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(datasub), datasub),
            file="edgeR_DE_control-vs-treat.txt",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
