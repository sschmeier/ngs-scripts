#!/usr/bin/env python
"""
NAME: check-mates-order.py
=========

DESCRIPTION
===========

Sanity check.

Take two mate pair fastq files:
1. sort each file based on id
2. Check if ids of reads in two files match

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

v0.1   2017/03/16    Initial version.

LICENCE
=======
2017, copyright Sebastian Schmeier (s.schmeier@gmail.com), sschmeier.com

template version: 1.6 (2016/11/09)
"""

from signal import signal, SIGPIPE, SIG_DFL
from Bio.SeqIO.FastaIO import FastaIterator
from Bio import SeqIO
import itertools
import sys
import os
import os.path
import argparse
import glob
import subprocess
import gzip
import bz2
import zipfile


# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
signal(SIGPIPE, SIG_DFL)

__version__ = 'v0.1'
__date__ = '2017/01/17'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'


def parse_cmdline():
    """ Parse command-line args. """
    ## parse cmd-line -----------------------------------------------------------
    description = 'For two paired-end fastq files, iterate over entries and check if ids are the same. Report if ids do not match.'
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument(
        'str_f1',
        metavar='fastq1.gz',
        help=
        'Forward reads.')
    parser.add_argument(
        'str_f2',
        metavar='fastq2.gz',
        help=
        'Backwards reads.')

    # if no arguments supplied print help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    return args, parser
    
def main():
    """ The main funtion. """
    args, parser = parse_cmdline()
    
    f1 = args.str_f1
    f2 = args.str_f2
    
    count = 0
    countUP = 0

    f1_iter = SeqIO.parse(gzip.open(f1, 'rt'), "fastq")
    f2_iter = SeqIO.parse(gzip.open(f2, 'rt'), "fastq")
    i = 0
    for rec1, rec2 in zip(f1_iter,f2_iter):
        i += 1
        # test if ids are same
        if rec1.id != rec2.id:
            sys.stderr.write("Entry %i: Mate1!=Mate2: %s\t%s\n" % (i, rec1.id, rec2.id))
            countUP += 1
            sys.stderr.write('Stopped checking at line %i\n'%i)
            return
        count += 1
    return


if __name__ == '__main__':
    sys.exit(main())

