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
                reader=read_tsv)

# The txi$abundance colSums might not sum up to 1e6 (as TPMs should),
# as some transcripts might be not asinged to genes 
# and are missing from the txi$abundance union values

# set colnames from samples.txt col3
colnames(txi$counts) <- samples[,3]

# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#edger
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
d <- DGEList(cts)
d$offset <- t(t(log(normMat)) + o)
# d is now ready for estimate dispersion functions see edgeR User's Guide
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

#datasub <- cbind(d$counts, cpm(d, log=TRUE))
datasub <- cpm(d, log=TRUE)

# print table
# problem always rownames column gets no header string
# we fix it with this workaround over a data.frame
write.table(data.frame("Genes"=rownames(datasub), datasub),
            file=paste(time, "salmon_edgeR_LOG2TMMTPM.txt", sep="_"),
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

