#!/usr/bin/Rscript
#
# Analysing counts for DGE.
#
#
# Expects two files:
#
#  counts.txt, header with samples, first column genes
#
#  samples.txt, sample information,
#               col1 which is used to estbalish groups
#               col3 will be used as colnames
#  e.g.
#
#    control,rep1,SRR4048970
#    control,rep2,SRR4048971
#    control,rep3,SRR4048972
#    15d-PGJ2,rep1,SRR4048973
#    15d-PGJ2,rep2,SRR4048974
#    15d-PGJ2,rep3,SRR4048975
#

args <- commandArgs(trailingOnly = TRUE)

if ( length(args)<1 ) {
   stop("USAGE:script.R counts.txt samples.txt")
}

file<-args[1]
samples<-args[2]

counts <- read.delim(file.path('.', file), sep="\t", header = TRUE, row.names=1)
samples <- read.table(file.path('.', samples), header = FALSE, sep=",")

# some information about the groupings of samples
group <- samples[,1]
group

library(edgeR)
library(methods)

colSums(counts)
summary(colSums(counts))
cat("\nrowsums:\n")
summary(rowSums(counts))

d <- DGEList(counts, group=group)

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
# WE ARE USING TPMCPM edger method
use = rowSums(cpm(counts) >1) >= min(table(d$samples$group))  # num smallest group size reps at least > 1 tpm

# apply filter
d <- d[use,]


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

design <- model.matrix(~d$samples$group)  # default if none is provided

# edgeR DE
d <- calcNormFactors(d, method="RLE")
d <- estimateDisp(d, design, robust=TRUE)

# MAKE REQUIRED TESTS

# 1
de.com <- exactTest(d, c("nonstim", "nonstim-mtb"))

# print some stats
cat("\n\n=========== Comp: nonstim ~ nonstim-mtb ===========\n")
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
dge <- dge[order(dge$FDR),]
#datasub <- cbind(dge, d$counts)
#datasub <- datasub[order(datasub$FDR),]

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(dge), dge),
            file="edgeR_DEG-nonstim_VS_nonstim-mtb.txt",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

# 2
de.com <- exactTest(d, c("nonstim", "ifng-mtb"))

# print some stats
cat("\n\n=========== Comp: nonstim ~ ifng-mtb ===========\n")
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
dge <- dge[order(dge$FDR),]
#datasub <- cbind(dge, d$counts)
#datasub <- datasub[order(datasub$FDR),]

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(dge), dge),
            file="edgeR_DEG-nonstim_VS_ifng-mtb.txt",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

# 3
de.com <- exactTest(d, c("nonstim", "il13-mtb"))
# print some stats
cat("\n\n=========== Comp: nonstim ~ il13-mtb ===========\n")
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
dge <- dge[order(dge$FDR),]
#datasub <- cbind(dge, d$counts)
#datasub <- datasub[order(datasub$FDR),]

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(dge), dge),
            file="edgeR_DEG-nonstim_VS_il13-mtb.txt",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)



# 4
de.com <- exactTest(d, c("nonstim", "il4-mtb"))
# print some stats
cat("\n\n=========== Comp: nonstim ~ il4-mtb ===========\n")
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
dge <- dge[order(dge$FDR),]
#datasub <- cbind(dge, d$counts)
#datasub <- datasub[order(datasub$FDR),]

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(dge), dge),
            file="edgeR_DEG-nonstim_VS_il4-mtb.txt",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)



# 5
de.com <- exactTest(d, c("nonstim", "il4il13-mtb"))
# print some stats
cat("\n\n=========== Comp: nonstim ~ il4il13-mtb ===========\n")
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
dge <- dge[order(dge$FDR),]
#datasub <- cbind(dge, d$counts)
#datasub <- datasub[order(datasub$FDR),]

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(dge), dge),
            file="edgeR_DEG-nonstim_VS_il4il13-mtb.txt",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
