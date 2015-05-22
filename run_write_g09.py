#!/usr/bin/env python

import input_g09
import sys
print 'run_write_g09 atoms coords method basis route charge spinmult title filename template'
print 'atoms: atomic symbols separated by , e.g. H,O,H'
print 'coords: coordinates of each atom separated by , and / e.g. 0,0,0/0,0,1/0,0,2'
print 'things past here are supported by keyword. i.e. name=value'
print 'method: hf, cas, b3lyp, etc. default: uwb97xd'
print 'basis: sto-6g, 6-31g, etc. default: aug-cc-pvtz'
print 'route: rest of route card. default: scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry '
print 'charge: charge of molecule. default: 0'
print 'spinmult: spin multiplicity. default: 1'
print 'title: title. default: title chemformula method basis'
print 'filename: filename. default: gaussian.com'
print 'template: preset routecard: default, opt, stable, freq, pop'
print 'xyz: xyz file of atom and coord (replaces given atom and coords)'

# need to change data type for some of these
atoms = sys.argv[1].split(',')
coords = [[float(number) for number in coord.split(',')] for coord in sys.argv[2].split('/')]
others = {}
for i in sys.argv[3:]:
    input = i.split('=')
    if input[0] == 'charge':
        others[input[0]] = int(input[1])
    elif input[0] == 'spinmult':
        others[input[0]] = int(input[1])
    else:
        others[input[0]] = input[1]
input = input_g09.G09Input(atoms, coords, **others)
input.write_input()
