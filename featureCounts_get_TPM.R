#!/usr/bin/Rscript
#
# Producing TPM (transcripts per million) normalised counts for featureCounts result table.
#
#
# WRITES TO STDOUT.
#
# Expects one file as command-line argument:
#
#  featurecounts produced count file -> here expects name featurecounts.genes.counts.gz
#
# featureCounts can be used like this:
# e.g.
# featureCounts -a Mus_musculus.GRCm38.85.gtf -o featurecounts.genes.counts star/*.bam
# featureCounts -a Mus_musculus.GRCm38.85.gtf -o featurecounts.tx.counts -g transcript_id star/*.bam
#
#
# ATTENTION: USES THE Length column of the featureCount table. This is not perfect as
# a better value would be effectivelength -> length of tx that a read can map to
# e.g. if read_len 30 and tx_len 100, effLength = tx_len - read_len + 1 = 71
#
library(methods)

# length in bases
counts2tpm <- function(counts, length)
    {
        c.rpk <- (counts/(length/1000)) # reads per kilobase
        c.rpksums <- colSums(c.rpk)/1e6 # scaling factors
        t(t(c.rpk)/c.rpksums) # scale rpks
    }

# START
args <- commandArgs(trailingOnly = TRUE)

if ( length(args)<1 ) {
   stop("USAGE: featureCounts_get_TPM.R featurecount-file")
}

file=args[1]

fc <-  read.table(file.path('.', file), header = TRUE, skip=1)
counts <- fc[,-(1:6)]
rownames(counts) <- fc$Geneid

tpm <- counts2tpm(counts, fc$Length)
colnames(tpm) <- paste(colnames(tpm),'TPM',sep='_')

ctpm <- cbind(counts, tpm)

write.table(data.frame("ID"=rownames(ctpm), ctpm),
            file="",
            append=FALSE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
