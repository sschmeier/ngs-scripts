#!/bin/bash
# Sort reads in fq1 and fq2 and join them based on id
# (thus expectes identical ids in pe file 1 and 2).
# Reads that do not join will be removed from both files.
# Attention, tabs in the files get converted to 4 spaces upfront.
#
# However, this can probably better be achieved with tools like SolexaQA++
#
# 
# USAGE: fastq-pe-join fqR1.gz fqR2.gz
fq1=$1
fq2=$2
# create tmpfile with process id as ending
tmpfile=$(mktemp /tmp/fastq-pe-join.tmp.$$)
join -t $'\t' -1 1 -2 1 <(zcat $fq1 | sed $'s/\t/    /g' | paste -d$'\t' - - - - | sort -t$'\t' -k1,1 -S 8G) <(zcat $fq2 | sed $'s/\t/    /g' | paste -d$'\t' - - - - | sort -t$'\t' -k1,1 -S 8G) > $tmpfile
cat $tmpfile | cut -f 1,2,3,4 | tr $'\t' $'\n' | gzip > $(basename $fq1 | cut -f 1 -d '.')_1.sorted.fastq.gz
cat $tmpfile | cut -f 1,5,6,7 | tr $'\t' $'\n' | gzip > $(basename $fq2 | cut -f 1 -d '.')_2.sorted.fastq.gz
rm "$tmpfile"
