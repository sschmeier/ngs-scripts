#!/usr/bin/env python
"""
NAME: xxx
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

0.1.0   xxxx/xx/xx    Initial version.

LICENCE
=======
2015, copyright Sebastian Schmeier (s.schmeier@gmail.com), http://sschmeier.com

template version: 1.1 (2015/12/10)
"""
__version__='0.1.0'
__date__='xxxx/xx/xx'
__email__='s.schmeier@gmail.com'
__author__='Sebastian Schmeier'
import sys, os, os.path, argparse, csv, collections, gzip, bz2, zipfile, time
## import pandas as pd ## non standard library to be imported

# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

def parse_cmdline():
    
    ## parse cmd-line -----------------------------------------------------------
    sDescription = 'Read delimited file.' 
    sVersion='version %s, date %s' %(__version__,__date__)
    sEpilog = 'Copyright %s (%s)' %(__author__, __email__)

    oParser = argparse.ArgumentParser(description=sDescription,
                                      version=sVersion,
                                      epilog=sEpilog)
    oParser.add_argument('sFile',
                         metavar='FILE',
                         help='Delimited file. [if set to "-" or "stdin" reads from standard in]')
    oParser.add_argument('-a', '--header',
                         dest='bHead',
                         action='store_true',
                         default=False,
                         help='Header in File. [default: False]')
    oParser.add_argument('-d', '--delimiter',
                         metavar='STRING',
                         dest='sDEL',
                         default='\t',
                         help='Delimiter used in file.  [default: "tab"]')
    oParser.add_argument('-f', '--field',
                         metavar='INT',
                         type=int,
                         dest='sFIELD',
                         default=1,
                         help='Field / Column in file to use.  [default: 1]')
    oParser.add_argument('-o', '--out',
                         metavar='STRING',
                         dest='sOut',
                         default=None,
                         help='Out-file. [default: "stdout"]')
    
    # PANDAS
    ## oParser.add_argument('-l', '--labels', 
    ##                      dest='iLabel',
    ##                      metavar="INT",
    ##                      type=int,
    ##                      default=None,
    ##                      help='Column number to use as labels. [default: None]')
    # FOR CONFIG-FILE PARSING
    ## oParser.add_argument('--config',
    ##                      dest = 'sConfig',
    ##                      metavar='CONFIG-FILE',
    ##                      default='config.ini',
    ##                      help='Config-file to read. [default: config.ini]')
    
    group1 = oParser.add_argument_group('Threading', 'Multithreading arguments:')
  
    group1.add_argument('-p', '--processes',
                         metavar='INT',
                         type=int,
                         dest='iP',
                         default=1,
                         help='Number of sub-processes (workers) to use. It is only logical to not give more processes than cpus/cores are available. [default: 1]')
    group1.add_argument('-t', '--time',
                         action='store_true',
                         dest='bTIME',
                         default=False,
                         help='Time the runtime and print to stderr. [default: False]')
    
    
    oArgs = oParser.parse_args()
    return oArgs, oParser

def load_file(s):
    """ LOADING FILES """
    if s in ['-', 'stdin']:
        oF = sys.stdin
    elif s.split('.')[-1] == 'gz':
        oF = gzip.open(s)
    elif s.split('.')[-1] == 'bz2':
        oF = bz2.BZFile(s)
    elif s.split('.')[-1] == 'zip':
        oF = zipfile.Zipfile(s)
    else:
        oF = open(s)
    return oF

def my_func(args):
    """
    THIS IS THE ACCTUAL WORKFUNCTION THAT HAS TO BE EXECUTED MULTPLE TIMES.
    This funion will be distributed to the cores requested.
    # do stuff
    res = ...
    return (args, res)
    """
    # Do stuff and create result
    # EXAMPLE: Here we add up arg1 and arg2 and wait a bit.
    res = args[0] + args[1]
    time.sleep(0.2)
    return (args, res)

def main():
    oArgs, oParser = parse_cmdline()

    # get field number to use in infile
    iF = oArgs.sFIELD - 1
    if iF < 0: oParser.error('Field -f has to be an integer > 0. EXIT.')
        
    oF = load_file(oArgs.sFile)

    if not oArgs.sOut:
        oFout = sys.stdout
    elif oArgs.sOut in ['-', 'stdout']:
        oFout = sys.stdout
    elif oArgs.sOut.split('.')[-1] == 'gz':
        oFout = gzip.open(oArgs.sOut, 'wb')
    else:
        oFout = open(oArgs.sOut, 'w')
        
    # -------------------------------------------------------
    # USING a config parser
    # Read config file if it exists
    # e.g.
    # [Classes]
    # c1 = 1
    # c2 = 3
    ## from ConfigParser import SafeConfigParser  
    ## dParams = {}
    ## if os.path.isfile(oArgs.sConfig):
    ##     oConfigParser = SafeConfigParser()
    ##     oConfigParser.read(oArgs.sConfig)
    ##     for section_name in oConfigParser.sections():
    ##         for name, value in oConfigParser.items(section_name):
    ##             dParams[name] = value                   
    # -------------------------------------------------------

    # ------------------------------------------------------
    # PANDAS approach
    # Check labels
    ## if oArgs.iLabel:
    ##     if oArgs.iLabel <= 0:
    ##         oParser.error('Label column number has to be > 0. EXIT.')
    ##     iLabel = oArgs.iLabel - 1
    ## else:
    ##     iLabel = False
    ## oDF = pd.read_csv(oF, sep=oArgs.sDEL, header=header, index_col=iLabel)
    # ------------------------------------------------------
   

    # ------------------------------------------------------
    #  THREADING
    # ------------------------------------------------------
    from timeit import default_timer as timer
    from multiprocessing import Pool
    
    # get number of subprocesses to use
    iNofProcesses = oArgs.iP
    if iNofProcesses<1: oParser.error('-p has to be > 0: EXIT.')
       
    #
    # FILL ARRAY WITH PARAMETER SETS TO PROCESS
    #
    # this array contains the total amount of jobs that has to be run
    aJobs = []

    # e.g. read tasks from file
    # delimited file handler
    oR = csv.reader(oF, delimiter = oArgs.sDEL)
    if oArgs.bHead:
        aH = oR.next()
    for a in oR:
        # EXAMPLE: 2 parameters = first two cols for my_func
        # EXAMPLE: TWO INTEGERS TO ADD UP
        try: 
            aJobs.append((int(a[0]), int(a[1])))
        except:
            oParser.error('Need 2 integers in column 1 and 2 to add up.')
    oF.close()
            
    # For timing
    fStart_time = timer()  # very crude
    # create pool of workers ---------------------
    pool = Pool(processes=iNofProcesses)

    #====================================================================
    # "chunksize"" usually only makes a noticable performance
    # difference for very large iterables
    # Here I set it to one to get the progress bar working nicly
    # Otherwise it will not give me the correct number of processes left
    # but chunksize number.
    chunksize = 1
    #====================================================================
    aResults = pool.map_async(my_func, aJobs, chunksize=chunksize)
    pool.close() # No more work

    iNumJobs = len(aJobs)
    # Progress bar
    #==============================
    # This can be changed to make progressbar bigger or smaller
    iProgressBarLength = 50
    #==============================
    while not aResults.ready():
        iNumNotDone = aResults._number_left
        iNumDone = iNumJobs-iNumNotDone
        iBarDone = iNumDone*iProgressBarLength/iNumJobs
        sBar = ('=' * iBarDone).ljust(iProgressBarLength)
        iPercent = int(iNumDone*100/iNumJobs)
        sys.stderr.write("JOBS (%i): [%s] %s%%\r" % (iNumJobs,
                                                     sBar,
                                                     str(iPercent).rjust(3)))
        sys.stderr.flush()
        time.sleep(0.1)  # wait a bit: here we test all .1 secs
    # Finish the progress bar
    sBar = '=' * iProgressBarLength
    sys.stderr.write("JOBS (%i): [%s] 100%%\r" % (iNumJobs, sBar))
    aResults = aResults.get()
    # --------------------------------------------

    fEnd_time = timer()
    if oArgs.bTIME: sys.stderr.write('RUNTIME(s): %.4f\n' %(fEnd_time - fStart_time))

    # Do stuff with the results
    ## for res in aResults:
    ##     oFout.write( '%s,%s:\t%s\n'%(str(res[0][0]),
    ##                                  str(res[0][1]),
    ##                                  str(res[1]) ))
    # ------------------------------------------------------

    oFout.close()
   
    return
        
if __name__ == '__main__':
    sys.exit(main())


#!/usr/bin/env python

# Import necessary libraries:
import csv
import os
import subprocess
import zipfile

# List modules used by FastQC:
modules = ['Basic_Statistics',
           'Per_base_sequence_quality',
           'Per_tile_sequence_quality',
           'Per_sequence_quality_scores',
           'Per_base_sequence_content',
           'Per_sequence_GC_content',
           'Per_base_N_content',
           'Sequence_Length_Distribution',
           'Sequence_Duplication_Levels',
           'Overrepresented_sequences',
           'Adapter_Content',
           'Kmer_Content']

# Set dict to convert module results to integer scores:
scores = {'pass': 1,
          'warn': 0,
          'fail': -1}

# Get current working directory:
cwd = os.getcwd()

# Get list of '_fastqc.zip' files generated by FastQC:
files = [file for file in os.listdir(cwd) if file.endswith('_fastqc.zip')]

# List to collect module scores for each '_fastqc.zip' file:
all_mod_scores = []

# Read fastqc_data.txt file in each archive:
for file in files:
    archive = zipfile.ZipFile(file, 'r') # open '_fastqc.zip' file
    members = archive.namelist() # return list of archive members
    fname = [member for member in members if 'fastqc_data.txt' in member][0] # find 'fastqc_data.txt' in members
    data = archive.open(fname) # open 'fastqc_data.txt'
    
    # Get module scores for this file:
    mod_scores = [file]
    for line in data:
        text = line.decode('utf-8')
        if '>>' in text and '>>END' not in text:
            text = text.lstrip('>>').split()
            module = '_'.join(text[:-1])
            result = text[-1]
            mod_scores.append(scores[result])
            
        # Append to all module scores list:
        all_mod_scores.append(mod_scores)
        
        # close all opened files:
        data.close()
        archive.close()
            
    # Write scores out to a CSV file:
    with open('all_mod_scores.csv', 'w') as f:
        writer = csv.writer(f)
        for mod_scores in all_mod_scores:
            writer.writerow(mod_scores)
        f.close()
