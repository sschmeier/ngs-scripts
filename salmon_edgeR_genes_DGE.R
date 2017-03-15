#!/usr/bin/Rscript
#
# Analysing salmon produced transcript counts for differential gene expression (DGE)
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
time <- format(Sys.time(), "%Y%m%d-%H%M%S")
library(tximport)
library(edgeR)
library(methods)
library(readr)

## START
# tx_gene_map.txt e.g.
t2g <- read.table(file.path('.', "tx_gene_map.txt"), header = FALSE)

# samples.txt, e.g.
samples <- read.table(file.path('.', "samples.txt"), header = FALSE)

# some information about the groupings of samples from samples.txt col1
group <- samples[,1]

# create desgin matrix
#design <- model.matrix(~group)

# load salmon files specified in samples.txt col4
files <- file.path(samples[,4])
# estimate adjusted counts from TPMs
# Based on tximport:
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
txi <- tximport(files = files,
                type="salmon",
                tx2gene = t2g,
                reader=read_tsv,
                countsFromAbundance="lengthScaledTPM")

# The txi$abundance colSums might not sum up to 1e6 (as TPMs should),
# as some transcripts might be not asinged to genes 
# and are missing from the txi$abundance union values

# set colnames from samples.txt col3
colnames(txi$counts) <- samples[,3]
colnames(txi$abundance) <- paste(samples[,3],'TPM',sep='_')

# get table with tximport gene-summed TPM and estimated counts 
ctpm <- cbind(txi$counts, txi$abundance)
write.table(data.frame("Genes"=rownames(ctpm), ctpm),
            file=paste(time, "edgeR_ALLGENES-COUNTS-TPM.txt", sep="_"),
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)


d = DGEList(counts=txi$counts, group=group)
d$tpm = txi$abundance

#--- DE ANALYSIS ---------
cat("\nBEFORE filtering stats:\n")
cat("colsums:\n")
colSums(d$counts)
summary(colSums(d$counts))
cat("\nrowsums:\n")
summary(rowSums(d$counts))

FILTER=1
if ( FILTER==1 ) {
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
 # However, we use TPM for filtering
 # could use abundance/TPM as a countoff here.
 use = rowSums(d$tpm >1) >= min(table(group))  # num smallest group size reps at least > 1 tpm
 #use = rowSums(cpm(d$counts) >1) >= min(table(d$samples$group))  # num smallest group size reps at least > 1 cpm

 # apply filter
 d <- d[use,]
 d$tpm <- d$tpm[use,]

 cat("\nAFTER filtering stats:\n")
 cat("colsums::\n")
 colSums(d$counts)
 summary(colSums(d$counts))
 cat("\nrowsums:\n")
 summary(rowSums(d$counts))

 # edgeR DE	
 # We could adjsut lib-sizes as oposed to original count table, problem?
 # d$samples$lib.size = colSums(d$counts)
}
#--------------------------------------------------------------------

# we are NOT normalising with edgeR.
# we use already adjusted counts from tximport
# also not clear if we need to do TMM between-sample normalisation or if this is done already through tximport
# supposedly we don't normalise again.
# d <- calcNormFactors(d, method="TMM") 
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

