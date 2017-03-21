#!/usr/bin/env python
"""
NAME: salmon_get_genes_rawLOG2TPM.py
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

0.1   20170321    Initial version.

LICENCE
=======
2017, copyright Sebastian Schmeier, s.schmeier@gmail.com

template version: 1.6 (2016/11/09)
"""
from timeit import default_timer as timer
from multiprocessing import Pool
from signal import signal, SIGPIPE, SIG_DFL
import numpy as np
import sys
import os
import os.path
import argparse
import csv
import collections
import gzip
import bz2
import zipfile
import time
import glob
import fnmatch



# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
signal(SIGPIPE, SIG_DFL)

__version__ = '0.1'
__date__ = '20170321'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'

def ihs(x):
    x = np.array(x)
    return np.log(x + (x**2+1)**0.5)    

def parse_cmdline():
    """ Parse command-line args. """
    ## parse cmd-line -----------------------------------------------------------
    description = 'Read quant.genes.sf salmon files from a specified directory  and extract TPM values for genes.'
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument(
        'str_dir',
        metavar='DIR',
        help=
        'Directory.')
    parser.add_argument('-o',
                        '--out',
                        metavar='STRING',
                        dest='outfile_name',
                        default=None,
                        help='Out-file. [default: "stdout"]')
    parser.add_argument(
        '--log',
        action='store_true',
        dest='log',
        default=False,
        help='Use log2-like transform (inverse hyperbolic sine). [default: False]')

    group1 = parser.add_argument_group('Threading',
                                       'Multithreading arguments:')

    group1.add_argument(
        '-p',
        '--processes',
        metavar='INT',
        type=int,
        dest='process_number',
        default=1,
        help=
        'Number of sub-processes (workers) to use.'+\
        ' It is only logical to not give more processes'+\
        ' than cpus/cores are available. [default: 1]')
    group1.add_argument(
        '-t',
        '--time',
        action='store_true',
        dest='show_runtime',
        default=False,
        help='Time the runtime and print to stderr. [default: False]')

    # if no arguments supplied print help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    return args, parser


def load_file(filename):
    """ LOADING FILES """
    if filename in ['-', 'stdin']:
        filehandle = sys.stdin
    elif filename.split('.')[-1] == 'gz':
        filehandle = gzip.open(filename)
    elif filename.split('.')[-1] == 'bz2':
        filehandle = bz2.BZFile(filename)
    elif filename.split('.')[-1] == 'zip':
        filehandle = zipfile.Zipfile(filename)
    else:
        filehandle = open(filename)
    return filehandle


def my_func(args):
    """
    THIS IS THE ACCTUAL WORKFUNCTION THAT HAS TO BE EXECUTED MULTPLE TIMES.
    This funion will be distributed to the cores requested.
    # do stuff
    res = ...
    return (args, res)
    """
    # Do stuff and create result
    fn = load_file(args)
    name = args.split('/')[1]
    reader = csv.reader(fn, delimiter='\t')
    header = reader.next()
    d = {}
    for a in reader:
        d[a[0]] = a[4]

    res = [name, d]
    return (args, res)


def main():
    """ The main funtion. """
    args, parser = parse_cmdline()
    matches = []
    for root, dirnames, filenames in os.walk(args.str_dir):
        for filename in fnmatch.filter(filenames, 'quant.genes.sf'):
            matches.append(os.path.join(root, filename))
    

    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif args.outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(args.outfile_name, 'wb')
    else:
        outfileobj = open(args.outfile_name, 'w')

    # ------------------------------------------------------
    #  THREADING
    # ------------------------------------------------------
    # get number of subprocesses to use
    process_number = args.process_number
    if process_number < 1:
        parser.error('-p has to be > 0: EXIT.')

    # FILL ARRAY WITH PARAMETER SETS TO PROCESS
    #
    # this array contains the total amount of jobs that has to be run
    job_list = matches

    # For timing
    start_time = timer()  # very crude
    # create pool of workers ---------------------
    pool = Pool(processes=process_number)

    # "chunksize"" usually only makes a noticable performance
    # difference for very large iterables
    # Here I set it to one to get the progress bar working nicly
    # Otherwise it will not give me the correct number of processes left
    # but chunksize number.
    chunksize = 1

    result_list = pool.map_async(my_func, job_list, chunksize=chunksize)
    pool.close()  # No more work

    jobs_total = len(job_list)
    # Progress bar
    #==============================
    # This can be changed to make progressbar bigger or smaller
    progress_bar_length = 50
    #==============================
    while not result_list.ready():
        num_not_done = result_list._number_left
        num_done = jobs_total - num_not_done
        num_bar_done = num_done * progress_bar_length / jobs_total
        bar_str = ('=' * num_bar_done).ljust(progress_bar_length)
        percent = int(num_done * 100 / jobs_total)
        sys.stderr.write("JOBS (%s): [%s] (%s) %s%%\r" % (str(num_not_done).rjust(len(str(jobs_total))),
                                                          bar_str,
                                                          str(num_done).rjust(len(str(jobs_total))),
                                                          str(percent).rjust(3)))
 
        sys.stderr.flush()
        time.sleep(0.1)  # wait a bit: here we test all .1 secs
    # Finish the progress bar
    bar_str = '=' * progress_bar_length
    sys.stderr.write("JOBS (%s): [%s] (%i) 100%%\n" % ('0'.rjust(len(str(jobs_total))),
                                                       bar_str,
                                                       jobs_total))
    result_list = result_list.get()
    # --------------------------------------------

    end_time = timer()
    if args.show_runtime:
        sys.stderr.write('\nPROCESS-TIME: %.4f sec' % (end_time - start_time))
	
	
    print_time_start = timer()
    if args.show_runtime:
        sys.stderr.write('\nWRITE-RESULTS...')
        
    # Do stuff with the results
    list = []
    genes = result_list[0][1][1].keys()
    genes.sort()
    names = ['Genes']
    names = names + [t[1][0] for t in result_list]
    outfileobj.write('%s\n' % ('\t'.join(names)))
    for gene in genes:
        vals = []
        for t in result_list:
            vals.append(float(t[1][1][gene]))

        if args.log:
            vals = ihs(vals)
            
        vals = [str(f) for f in vals]
        outfileobj.write('%s\t%s\n'%(gene,'\t'.join(vals)))
            
    print_time_stop = timer()
    if args.show_runtime:
        sys.stderr.write(' %.4f sec\n' % (print_time_stop - print_time_start))		
    # ------------------------------------------------------
    outfileobj.close()
    return


if __name__ == '__main__':
    sys.exit(main())

