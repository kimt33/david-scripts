import sys
import os

cwd = os.getcwd()

def farnaz_format(dirname, filename):
    """ generates gaussian files organized by farnaz
    
    level0: molecule grouping
    level1: molecule identifier
    level2: method/basis set
    level3: output files (gaussian.com, gaussian.log, gaussian.out, gaussian.chk)

    ASSUMES:
        you are on level0
    Args:
        dirname: string that is the name of the molecule grouping directory
    """
    for dirpath, dirnames, filenames in os.walk(dirname):
        if len(dirpath.split(os.sep))==3:
            if filename in filenames:
                yield os.path.join(dirpath, filename)
            else:
                print 'NOT FOUND:',dirpath
