#!/usr/bin/env python
"""
NAME: fastQCstat
================

DESCRIPTION
===========
Get summary statistics for many fastqc result zip-files in a directory (including all subfolders).

INSTALLATION
============
Copy file into directory.

USAGE
=====
Use `python fastQCstat.py DIR`.

VERSION HISTORY
===============

0.1.1   2016/03/04    Adjust output.
0.1.0   2016/03/03    Initial version.

LICENCE
=======
2016, copyright Sebastian Schmeier (s.schmeier@gmail.com), http://sschmeier.com

template version: 1.1 (2015/12/10)
"""
__version__='0.1.1'
__date__='2016/03/04'
__email__='s.schmeier@gmail.com'
__author__='Sebastian Schmeier'
import sys, os, os.path, argparse, csv, collections
import gzip, bz2, zipfile, time, glob, datetime

# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

def parse_cmdline():

    ## parse cmd-line -----------------------------------------------------------
    sDescription = 'Find all fastqc result zip files in directory (including sub-folders) and compile report.'
    sVersion='version %s, date %s' %(__version__,__date__)
    sEpilog = 'Copyright %s (%s)' %(__author__, __email__)

    oParser = argparse.ArgumentParser(description=sDescription,
                                      epilog=sEpilog)

    oParser.add_argument('--version', action='version', version='%s'%(sVersion))
    oParser.add_argument('sDir',
                         metavar='DIR',
                         type=str,
                         nargs='?',
                         default=os.getcwd(),
                         help='Directory to start looking for fastqc result zip-files. [DEFAULT: "."]')
    oParser.add_argument('-p', '--prefix',
                        dest='sOutPrefix',
                         metavar='PREFIX',
                         type=str,
                         nargs='?',
                         default=datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d-%H%M%S'),
                         help='Prefix to use for result-files. [DEFAULT: current date/time]')

    oArgs = oParser.parse_args()
    return oArgs, oParser

def main():
    oArgs, oParser = parse_cmdline()

    # collect all files
    aFiles = []
    for dirpath, dirnames, filenames in os.walk(oArgs.sDir):
        aFiles.extend(glob.glob(os.path.join(dirpath,'*_fastqc.zip')))
    
    if len(aFiles) <= 0:
        oParser.error('No fastqc-zip files found in %s. EXIT.'%oArgs.sDir)
    
    dModulesFiles = collections.OrderedDict()
  
    # Set dict to convert module results to integer scores:
    scores = {'pass': 1, 'warn': 0, 'fail': -1}

    # List to collect module scores for each '_fastqc.zip' file:
    dMod_scores = collections.OrderedDict()
    
    # Read fastqc_data.txt file in each archive:
    for sFile in aFiles:
        archive = zipfile.ZipFile(sFile, 'r') # open '_fastqc.zip' file
        aFiles_inarchive = archive.namelist() # return list of archive members
        
        # find 'fastqc_data.txt' in members
        fname = [member for member in aFiles_inarchive if 'fastqc_data.txt' in member][0] 
        data = archive.open(fname) # open 'fastqc_data.txt'
        
        # Get module scores for this file:
        dMod_scores[sFile] = collections.OrderedDict()
        
        for line in data:
            text = line.decode('utf-8')
            if '>>' in text and '>>END' not in text:
                aText = [s.strip() for s in text.lstrip('>>').split()]
                module = '_'.join(aText[:-1])
                result = aText[-1]
                # collect for the file the found module score
                dMod_scores[sFile][module] = scores[result]
                
                if module not in dModulesFiles:
                    dModulesFiles[module] = {'pass': [], 'warn': [], 'fail': []}
                    
                dModulesFiles[module][result].append(sFile) # safe file

        # close all opened files:
        data.close()
        archive.close()

    # Write scores out to a CSV file:
    with open('%s_fastQCstat_table.txt'%(oArgs.sOutPrefix), 'w') as f:
        oWriter = csv.writer(f, delimiter = '\t')
        bHead = 1
        for sF in dMod_scores:
            a = [sF]
            aModules = dMod_scores[sF].keys()
            assert len(aModules) == 12
            aModules.sort()
            # write header?
            if bHead:
                aHeader =  ['File'] + aModules
                oWriter.writerow(aHeader)
                bHead = 0

            for sMod in aModules:
                a.append(dMod_scores[sF][sMod])

            oWriter.writerow(a)
    f.close()

    iFiles = len(aFiles) # number of files analysed
    with open('%s_fastQCstat_summary.txt'%(oArgs.sOutPrefix), 'w') as f:
        oWriter = csv.writer(f, delimiter = '\t')
        assert len(dModulesFiles.keys()) == 12
        
        oWriter.writerow(['PASS','WARN','FAIL','MODULE'])
        for sMod in dModulesFiles:
            iLenPass = len(dModulesFiles[sMod]['pass'])
            iLenWarn = len(dModulesFiles[sMod]['warn'])
            iLenFail = len(dModulesFiles[sMod]['fail'])
            
            a = [iLenPass, iLenWarn, iLenFail, sMod]
            oWriter.writerow(a)

        f.write('>>># files analysed: %i\n'%iFiles)
        for sMod in dModulesFiles:
            if len(dModulesFiles[sMod]['fail']) > 0:
                f.write('>>>FAILED %s:\n'%(sMod))
                for sF in dModulesFiles[sMod]['fail']:
                    f.write('%s\n'%sF)
            
    f.close()
    return

if __name__ == '__main__':
    sys.exit(main())
