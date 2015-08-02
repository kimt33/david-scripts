import horton.io.smartreader as sr
import re
def extract_xyz(filename):
    regeometry = re.compile(r'^\s+Final geometry\s*$')
    recoords = re.compile(r'^\s*(\w+)\s+\d+\s+([\d\.\+\-EeDd]+)\s+([\d\.\+\-EeDd]+)\s+([\d\.\+\-EeDd]+)\s*$')
    reend = re.compile(r'^\s*Iter\s+Energy\s+Change')
    header = sr.Header(regeometry)
    coords = sr.Item(recoords, header=header, max_match=-1, reg_ender=reend)
    with open(filename, 'r') as f:
        for line in f:
            coords(line)
    temp = zip(*coords.value[0])
    atoms = list(temp[0])
    coords = zip(*temp[1:])
    coords = [[float(i) for i in j] for j in coords]
    #printable = '\n'.join('{0}{1}'.format(i,'{0:>12f}{1:>12f}{2:>12f}'.format(*j)) for i,j in zip(atoms,coords))
    return atoms, coords

