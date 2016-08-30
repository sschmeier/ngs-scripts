#!/bin/bash
# USAGE: script.sh in-dir out-dir genomesize
# recommendations taken from: http://bib.oxfordjournals.org/content/early/2016/01/12/bib.bbv110.full
#
# This script is only working for ONE treatment file, without any control files.
# Not suitable for many treatment files or control files.

dir=$1;
outdir=$2;
genomesize=$3;

now=$(date +"%Y-%m-%d_%H%M%S")

# set error and stdout file
errfile=${outdir}/macs2.${now}.stderr
outfile=${outdir}/macs2.${now}.stdout

outdir_2=${outdir}/macs2.${now};
mkdir ${outdir_2};

for i in `ls ${dir}/*.bam`; do

    echo ${i} >> ${errfile};
    echo ${i} >> ${outfile};

    # MACS2
    # -t TFILE [TFILE ...], --treatment TFILE [TFILE ...]
    # -n NAME, --name NAME  Experiment name, which will be used to generate output file names. DEFAULT: "NA"
    # -f type
    # -g GSIZE, --gsize GSIZE
    #                    Effective genome size. It can be 1.0e+9 or 1000000000,
    #                    or shortcuts:'hs' for human (2.7e9), 'mm' for mouse
    #                    (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for
    #                    fruitfly (1.2e8), Default:hs
    # -q QVALUE, --qvalue QVALUE
    #                    Minimum FDR (q-value) cutoff for peak detection.
    #                    DEFAULT: 0.05. -q, and -p are mutually exclusive.
    # --call-summits     If set, MACS will use a more sophisticated signal
    #                    processing approach to find subpeak summits in each
    #                    enriched peak region. DEFAULT: False
    macs2 callpeak -n ${outdir_2}/$(basename ${i} | sed 's/.bam//g') \
                   -t ${i} \
                   -f BAM \
                   -g ${genomesize} \
                   -q 0.05 \
	           --call-summits 2>> ${errfile} >> ${outfile};

    echo "----------" >> ${errfile};
    echo "----------" >> ${outfile};
done;
 
gzip ${outdir_2}/*.narrowPeak
gzip ${outdir_2}/*.bed
rm ${outdir_2}/*.xls
rm ${outdir_2}/*.r
