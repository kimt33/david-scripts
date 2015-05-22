import sys
from read_log import LogData


if sys.argv[1] == 'logdata':
    a=LogData(inpfile=sys.argv[2],program='adf')
    print a.occ
