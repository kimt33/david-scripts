import numpy as np
import loadpaths
import name_mol
from horton.io import data_basis as db
from horton.context import context
import os
import re

class G09Input:
    '''class that contains g09 input parameters which allows creation of it
    '''
    def __init__(self, atoms, coords, method='', basis='', route='', charge=0, spinmult=1, title='', filename='', template='default'):
        '''initializes input object using default settings

        atoms
            array of chemical symbols for atoms
            in proper capitalization
        coords
            array of coordinates (array of dim 3)
        NOTE: Defaults were arbitrary and depended only on the template I used to create it
        '''
        if filename!='':
            self.filename=filename
        else:
            self.filename='gaussian.com'
        self.basename=self.filename[:-4] #NOTE: assume specified filename has three character extension
        #Link0
        self.chk=self.basename+'.chk'
        self.mem='1500MB'
        #Route section
        if method!='':
            self.method=method
        else:
            self.method='uwb97xd'
        if basis!='':
            self.basis = basis
        else:
            self.basis = 'aug-cc-pvtz'
        # this needs to be better
        if template == 'default':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry ' + route
        elif template == 'opt':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry opt=tight ' + route
        elif template == 'stable':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry stable=opt ' + route
        elif template == 'freq':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry freq ' + route
        elif template == 'pop':
            route = 'scf=(tight,xqc,fermi) integral=grid=ultrafine nosymmetry density=current pop=(chelpg,npa) IOp(6/80=1) ' + route
        route = route.split()

        self.keywords = {}
        self.routerest = []
        for i in route:
            if 'opt' == i.lower()[:3]:
                self.keywords['opt'] = i[4:]
            elif 'scf' == i.lower()[:3]:
                self.keywords['scf'] = i[4:]
            elif 'integral' == i.lower()[:8]:
                self.keywords['integral'] = i[9:]
            elif 'density' == i.lower()[:7]:
                self.keywords['density'] = i[8:]
            elif 'pop' == i.lower()[:3]:
                self.keywords['pop'] = i[4:]
            elif 'freq' == i.lower()[:4]:
                self.keywords['freq'] = i[5:]
            elif 'stable' == i.lower()[:6]:
                self.keywords['stable'] = i[7:]
            elif 'nosymmetry' == i.lower():
                self.keywords['nosymmetry'] = ''
            else:
                self.routerest.append(i)

        #title section
        chemformula = name_mol.chemicalformula(atoms)
        self.title=title+' '+chemformula+' '+self.method+'/'+self.basis
        #molecule specification (default: 0 1)
        self.charge=charge #difference in charge from neutral
        self.spinmultiplicity=spinmult #1 for singlet, 2 for doublet, ...
        #atoms
        self.atoms=atoms
        self.coords=coords

    def write_input(self,filename=''):
        if filename == '':
            filename = self.filename
        with open(filename,'w') as fp:
            fp.write('%chk='+self.chk+'\n')
            fp.write('%mem='+self.mem+'\n')
            # route
            fp.write('#p '+self.method+'/gen ')
            routerest = ''
            for keyword, val in self.keywords.iteritems():
                if val == '':
                    routerest += keyword + ' '
                else:
                    routerest += keyword + '=' + val + ' '
            routerest += ' '.join(self.routerest)
            fp.write(routerest+'\n\n')
            # title
            fp.write(self.title+'\n\n')
            # mol spec
            fp.write(str(self.charge)+' '+str(self.spinmultiplicity)+'\n')
            # coord
            for atom,coord in zip(self.atoms,self.coords):
                line = atom + '{0:>17f}{1:>17f}{2:>17f}'.format(*coord)
                fp.write(line+'\n')
            fp.write('\n')
            # basis set
            for ext in ['gbs','nwchem','davidh5']:
                if os.path.isfile(context.get_fn('basis/{0}.{1}'.format(self.basis,ext))):
                    path = context.get_fn('basis/{0}.{1}'.format(self.basis,ext))
                    break
            else:
                raise AssertionError, 'given basis {0} does not exist in horton database'.format(self.basis)
            basisinfo = db.DataBasis(path, fileformat=os.path.splitext(path)[1][1:])
            fp.write(basisinfo.write_gbs(self.atoms))
            fp.write('\n\n\n')

    def scan_geometry_atom(self,atomindex,coord1,coord2,numsteps):
        '''returns the generator for the input object where the coordinate of atom (specified by atomindex) changes incrementally following the path of coord1 to coord2

        '''
        coord1=np.array(coord1)
        coord2=np.array(coord2)
        stepvec=(coord2-coord1)/(numsteps-1)
        coord_range=[coord1+i*stepvec for i in range(numsteps)]
        for coord in coord_range:
            self.update_coord_atom(coord.tolist(),atomindex)
            yield self
'''
a=g09_input(['H','O'],[[0.0,0.0,0.0],[2.4,0.0,0.0]])
a.write_input('testest.txt')
'''
