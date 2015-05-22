import sys
import os
import re
import glob
import read_log as log
import input_adf
import input_g09

print 'python autorun.py flag'
print 'this script is run from a base directory'
print 'flag = int.int.int... where each int specifies some flag (separator is .)'



flag = sys.argv[1]
flags = flag.split('.')
cwd = os.getcwd()
def create_dirs(list_str):
    """creates directories depending on the specified levels
    
    Args:
        list_str: list of strings that describe the directory
    Return:
        None
    """
    if len(list_str) >= 1: 
        if list_str[0] == 'help' or list_str == 'help':
            	all succeeding strings describe the name of the directory to be created, with the exception that the base is 0'
        if list_str[0] =='0':
            os.mkdir('base')
        else:
            try:
                levels = ['base']+['level'+str(i) for i in range(1,int(list_str[0])+1)]
            except ValueError:
                levels = ['base']+list_str[1:]
            for i,j in enumerate(levels):
                path = os.path.join(cwd,*levels[:i+1])
                if not os.path.isdir(path):
                    os.mkdir(path)
def get_list_files(level,ext):
    """returns list of files in provided level with extension given

    Args:
        level: string of form 'leveln' where n is an integer, or list of directories that describe the path
        ext: extension of files to find of form (.ext)
    Return:
        list of strings (paths)
    """
    if type(level) is str and re.search('^level\d+$',level):
        levels = ['base'] + ['level'+str(i) for i in range(1,int(level[5:])+1)]
    elif type(level) is list:
        levels = level
    else:
        raise AssertionError, 'given levels is bad'
    path = os.path.join(cwd,*levels)
    files = os.path.join(path,'*'+ext)
    return glob.glob(files)
def extract_files_log(program,extension,level):
    """extract the log file information from the given files (set by level and extension)

    Args:
    	program: string \in {g09, adf}
    	extension: string 
    	level: string

    """
    filepaths = get_list_files(level, ext)
    data = []
    for i in filepaths:
        data.append(log.LogData(inpdata=i,program=prog))
    return data
def make_run(program,extension,level):
    """extract the log file information from the given files (set by level and extension)

    Args:
    	program: string \in {g09, adf}
    	extension: string 
    	level: string

    """
    filepaths = get_list_files(level, ext)
    data = []
    for i in filepaths:
        data.append(log.LogData(inpdata=i,program=prog))
    return data

#create dirs
if flags[0] == '1':
    list_str = sys.argv[2].split('.')
    create_dirs(list_str)

#extract file information
if flags[0] == '2':
    prog = sys.argv[2]
    ext = sys.argv[3]
    level = sys.argv[4]
    data = extract_files_log(prog, ext, level)
    #something with data

#create input
if flag[0] == '3':
    prog = sys.argv[2]
    filename = sys.argv[3]
    template = sys.argv[4]
    if prog == 'adf':
        input = input_adf.ADF_input(template=template)
    elif prog =='g09':
        input = input_g09.g09_input(template=template)
    input.write_input(filename)

        


