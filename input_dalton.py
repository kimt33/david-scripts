class DaltonInput:
    '''class that contains Dalton input parameters which allows creation of it
    '''
    def __init__(self, atoms, coords, basis='', symmetry=0, title='', filename=''):
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
            self.filename='dalton.mol'
        self.basename=self.filename[:-4] #NOTE: assume specified filename has three character extension
        
        if basis!='':
            self.basis = basis
        else:
            self.basis = 'aug-cc-pvtz'

        #title section
        chemformula = ''.join(atoms)
        self.title=title+' '+chemformula+' '+self.basis
        #atoms
        self.atoms=atoms
        self.coords = [[float(j) for j in i] for i in coords]
        # stupid things
        self.num_atoms = {i:self.atoms.count(i) for i in set(self.atoms)}
        self.atom_charge = {'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9,
                            'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17,
                            'Ar':18, 'K':19, 'Ca':20, 'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25,
                            'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32, 'As':33,
                            'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, 'Nb':41,
                            'Mo':42, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50,
                            'Sb':51, 'Te':52, 'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58,
                            'Pr':59, 'Nd':60, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67,
                            'Er':68, 'Tm':69, 'Yb':70, 'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75,
                            'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83,
                            'Th':90, 'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97,
                            'Cf':98, 'Es':99, 'Fm':100, 'Md':101, 'No':102, 'Lr':103, 'Rf':104,
                            'Db':105, 'Sg':106, 'Bh':107, 'Hs':108, 'Mt':109, 'Ds':110, 'Rg':111,
                            'Cp':112, 'Uut':113, 'Uuq':114, 'Uup':115, 'Uuh':116, 'Uus':117, 'Uuo':118}
        # symmetry
        self.symmetry = symmetry

    def write_input(self,filename=''):
        if filename == '':
            filename = self.filename
        with open(filename,'w') as fp:
            # header
            fp.write('ATOMBASIS\n')
            # title
            fp.write(self.title+'\n\n')
            # something
            fp.write('Generators={0} Atomtypes={1}\n'.format(self.symmetry, len(self.atoms)))
            # coords of each atom
            for atom in set(self.atoms):
                fp.write('Charge={0} Atoms={1} Basis={2}\n'.format(self.atom_charge[atom], self.num_atoms[atom], self.basis))
                for i, coord in enumerate(self.coords):
                    if self.atoms[i] == atom:
                        fp.write('{0} {1:>17f} {2:>17f} {3:>17f}\n'.format(atom, coord[0], coord[1], coord[2]))
            fp.write('\n\n\n')
