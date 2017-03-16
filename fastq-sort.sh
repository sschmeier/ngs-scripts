#!/bin/bash
# USAGE: cat fastq | fastq-sort.sh > fastq.sorted

# input from stdin if no file is given
# First I convert tabs to 4 spaces, last I convert them back to tab
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"
cat $input | sed $'s/\t/    /g' | paste -d$'\t' - - - - | sort -t$'\t' -k1,1 -S 8G | tr $'\t' $'\n' | sed $'s/    /\t/g'
