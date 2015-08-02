import re
import sys

def extractxyz_adflog(filename,coordtype='ang',start=0):
    #may need to specify starting point
    #amd does something weird with coordinatea
    infile=open(filename ,'r')
    infile.seek(start,0)
    while 1:
        line=infile.readline()
        recoord=re.search(r'^\s*Coordinates \(Cartesian\)\s*$',line)
        atoms=[]
        bohrxyz=[]
        angstromxyz=[]
        if recoord is not None:
            for i in range(6):
                line=infile.readline()
            reonexyz=re.search(r'^\s*\d+\s+a(\w+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+\d+\s+\d+\s+\d+\s*$',line)
            while reonexyz is not None:
                atoms.append(reonexyz.group(1))
                bohrxyz.append([reonexyz.group(2),reonexyz.group(3),reonexyz.group(4)])
                angstromxyz.append([reonexyz.group(5),reonexyz.group(6),reonexyz.group(7)])
    if coordtype=='ang':
        return atoms,angstromxyz
    elif coordtype=='bohr':
        return atoms,bohrxyz
            
print extractxyz_adflog(sys.argv[1])
