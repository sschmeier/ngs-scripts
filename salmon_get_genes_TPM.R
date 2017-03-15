#!/usr/bin/Rscript
#
# Summarize salmon produced transcript counts and TPM for genes.
# Will output tximport transformed/normalised counts.
#
# WRITES TO STDOUT.
# 
# Expects two files: 
#
#  tx_gene_map.txt, transcript to gene mapping, e.g.
#    ENSMUST00000205326.1    ENSMUSG00000108652
#    ENSMUST00000206672.1    ENSMUSG00000108652
#    ENSMUST00000021676.11   ENSMUSG00000021252
#    ENSMUST00000124311.1    ENSMUSG00000021252
#    ...
#
#  samples.txt, sample information,
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
library(readr)

# tx_gene_map.txt e.g.
t2g <- read.table(file.path('.', "tx_gene_map.txt"), header = FALSE)

# samples.txt, e.g.
samples <- read.table(file.path('.', "samples.txt"), header = FALSE)

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

# set colnames from samples.txt col3
colnames(txi$counts) <- samples[,3]
colnames(txi$abundance) <- paste(samples[,3],'TPM',sep='_')

# get table with tximport gene-summed TPM and estimated counts 
ctpm <- cbind(txi$counts, txi$abundance)
write.table(data.frame("Genes"=rownames(ctpm), ctpm),
            file="",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
