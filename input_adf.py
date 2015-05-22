import numpy as np
class ADF_input:
    '''class that contains aim input parameters which allows creation of it
    '''
    def __init__(self, atoms=[], coords=[], template=None):
        '''initializes input object using default settings

        atoms
            array of chemical symbols for atoms
            in proper capitalization
        coords
            array of coordinates (array of dim 3)
        NOTE: Defaults were arbitrary and depended only on the template I used to create it
        '''
        #atoms
        self.atoms=atoms
        self.coords=coords
        #unrestricted
        self.unrestricted=1
        #charge
        self.charge=0
        self.spinmultiplicity=0
        #basis
        self.basis_type='DZP'
        self.basis_core='None'
        self.basis_fittype='ZORA/QZ4P'
        self.basis_createoutput='None'
        #integration
        self.integration='6.0 5.0'
        #xc
        self.exchange='GGA Becke Perdew'
        #relativistic
        self.relativistic='Scalar ZORA'
        #geometry
        self.geometry=''
        #occupation
        self.occ = []
        #save
        self.save='TAPE21 TAPE13'
        #fullscf
        self.fullscf=1
        #bader
        self.bader=1
        #bondorder
        self.bondorder=1
        #scf
        self.lshift=''
        self.adiis=0
        #exact density
        self.exactdensity=0
        #electron smearing
        self.electronsmearing=0

    def update_integration(self,*nums):
        self.integration=' '.join(str(float(num)) for num in nums)
    def update_spinmultiplicity(self,num):
        self.spinmultiplicity=num
    def update_fullscf(self,flag):
        self.fullscf=flag
    def update_bader(self,flag):
        self.bader=flag
    def update_lshift(self,num):
        self.lshift=str(float(num))
    def update_adiis(self,flag):
        self.adiis=flag
    def update_exactdensity(self,flag):
        self.exactdensity=flag
    def update_coord_atom(self,coord,atomindex):
        self.coords[atomindex]=coord
    def update_coords(self,coords):
        self.coords = coords
    def update_occupations(self, occ):
        self.occ = occ
    def update_electronsmearing(self, flag):
        self.electronsmearing = flag


    def write_input(self,filename):
        with open(filename,'w') as fp:
            #atoms
            fp.write('ATOMS\n')
            for atom,coord in zip(self.atoms,self.coords):
                line=atom+' '*(3-len(atom))+' '.join(str(float(i)) for i in coord)
                fp.write(line+'\n')
            fp.write('END\n\n')
            #unrestricted
            if self.unrestricted!=0:
                fp.write('UNRESTRICTED\n\n')
            #charge
            fp.write('CHARGE '+str(self.charge)+' '+str(self.spinmultiplicity)+'\n\n')
            #basis
            fp.write('BASIS\n')
            fp.write('type '+self.basis_type+'\n')
            fp.write('core '+self.basis_core+'\n')
            fp.write('FitType '+self.basis_fittype+'\n')
            fp.write('createoutput '+self.basis_createoutput+'\n')
            fp.write('END\n\n')
            #integration
            fp.write('INTEGRATION '+self.integration+'\n\n')
            #xc
            fp.write('XC\n')
            fp.write(self.exchange+'\n')
            fp.write('END\n\n')
            #exact density
            if self.exactdensity!=0:
                fp.write('EXACTDENSITY\n\n')
            #occupations
            if self.electronsmearing==0:
                fp.write('OCCUPATIONS Smearq=0\n\n')
            #occupations (ASSUMES C INFTY V)
            if self.occ != []:
                fp.write('OCCUPATIONS\n')
                for i, j in zip(['sigma','pi','delta','phi'],self.occ):
                    fp.write(' '+i+' '+str(j[0])+' // '+str(j[1])+'\n') 
                fp.write('END\n\n')
            #relativistic
            fp.write('RELATIVISTIC '+self.relativistic+'\n\n')
            #geometry
            if self.geometry!='':
                fp.write('GEOMETRY\n')
                fp.write(self.geometry+'\n')
                fp.write('END\n\n')
            #scf
            fp.write('SCF\n')
            if self.lshift!='':
                fp.write('lshift '+self.lshift+'\n')
            if self.adiis!=0:
                fp.write('NoADIIS\n')
            fp.write('END\n\n')


            #save
            fp.write('SAVE '+self.save+'\n\n')
            #fullscf
            if self.fullscf!=0:
                fp.write('FULLSCF\n\n')
            #bader
            if self.bader!=1:
                fp.write('BADER\n\n')
            #bond order analysis
            if self.bondorder!=1:
                fp.write('BONDORDER tol=0.05 printall\n\n')

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
a=ADF_input(['Th','O'],[[0.0,0.0,0.0],[2.4,0.0,0.0]])
for i in a.scan_geometry_1d(1,[2.4,0,0],[3,0,0],6):
    print i.coords
'''
