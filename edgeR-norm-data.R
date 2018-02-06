#!/usr/bin/Rscript
#
# Transform counts to RLE or TMM normalised counts.
#
# argument
#file       - Infile, header, genes in rows, samples in columns 
#RLE/TMM 
#TRUE/FALSE - should robust method for dispiersion estimation be used
#OUTFILE

args <- commandArgs(trailingOnly = TRUE)
if ( length(args)<4 ) {
   stop("USAGE: script.R infile RLE/TMM TRUE/FALSE OUTFILE")
}

file <- args[1]
norm <- args[2]
robu <- as.logical(args[3])
out  <- args[4]

counts <- read.delim(file.path('.', file), sep="\t", header = TRUE, row.names=1)

library(edgeR)
library(methods)

d <- DGEList(counts)

# edgeR DE
d <- calcNormFactors(d, method=norm)
d <- estimateDisp(d, robust=robu)

countsNorm <- cpm(d)

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(countsNorm), countsNorm),
            file=out,
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
