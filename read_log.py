import sys
import re
import horton.io.smartreader as sr

class LogData:
    def __init__(self,inpfile=None,program=''):
        '''extracts information from some log file and stores it as an object
        '''

        self.filename=inpfile
        self.jobname=''

        self.basis_type=''#wording!
        self.basisset=''

        self.num_atoms=0 #necessary? use len(atoms) instead?
        self.atoms=[]
        self.charge_atoms=[] #atomic number (should change name)
       
        self.coordinates=[] #in ANGSTROM
        self.charge_molecule=0 #notsupported by wfn reader
        
        self.pop_mulliken=[]
        self.pop_esp=[]
        self.pop_npa=[]

        self.geomcycle=[]
        self.energycycle=[]
        
        self.lda=None
        self.gradcorr=None
        self.relcorr=None
        self.coretreat=None

        self.num_elecalpha=None
        self.num_elecbeta=None
        self.num_elec=None
        self.totalspin=None
        self.occ = []  # list of dictionary that corresponds to each symmetry (one for each spin)

        self.energy_bond=None #bond energy: energy difference between molecule and fragments (fragments are atoms)

        self.warnings=[]
        self.errors=[]
        if program.lower()=='adf':
            self.extract_adf(inpfile)
        else:
            print 'ERROR: program type not recognized'
    def extract_adf(self,inpfile):
        '''extracts from log/out file of adf2013
    
        need to be careful with units, as it is not currently stored
        '''
        #find sections
        sections = {}
        with open(inpfile,'r') as infile:
            line = infile.readline()
            redivider = re.compile(r'^\s*=+\s*$')
            resection = re.compile(r'^(.+[\w\(\)\-\,\+\*\.:] +\w.*)$')
            while line != '':
                if redivider.search(line) is not None:
                    line = infile.readline()
                    possible = resection.search(line)
                    if possible is not None:
                        line = infile.readline()
                        if redivider.search(line) is not None:
                            temp_name = re.sub(r' ', '', possible.group(1))
                            temp_name = re.sub(r'\*\*\*.+$', '', temp_name)
                            temp_name = re.sub(r' +', ' ', temp_name)
                            temp_name = temp_name.lower()
                            sections[temp_name] = possible.group(1)
                line = infile.readline()

        with open(inpfile ,'r') as infile:
            def f(regexpobj,line):
                temp=regexpobj.search(line)
                if temp is not None:
                    try:
                        return temp.group(1)
                    except IndexError:
                        return True
                else:
                    return None
            line=infile.readline()
            flag_firstgeo=0 #regular expression for first geometry matches all the subsequent ones

            #reg exp objs
            rebasic=re.compile(r'^\s*M O D E L   P A R A M E T E R S\s*$')
            respin=re.compile(r'^\s*S Y M M E T R Y ,   E L E C T R O N S\s*$')
            regeometry=re.compile(r'^\s*G E O M E T R Y\s*\*')
            regeomcycle=re.compile(r'^\s*Geometry CYCLE\s+\d+\s*$')
            resymmetry = re.compile(r'^\s*Orbital Energies, all Irreps, both Spins\s*')
            reenergy=re.compile(r'^\s*B O N D I N G   E N E R G Y\s*')
            relog=re.compile(r'^\s*\(LOGFILE\)\s*$')
            while line!='':
                if f(rebasic,line):
                    reend=re.compile(r'^\s*\*+\s*$')
                    relda=re.compile(r'^\s*LDA:\s+(\w+)')
                    regradcorr=re.compile(r'^\s*Gradient Corrections:\s+(([\w\d]+ )+)\s+')
                    reunrestricted1=re.compile(r'\s*SPIN ')
                    reunrestricted2=re.compile(r'\s*OTHER ASPECTS\s*$')
                    reunrestricted3=re.compile(r'^\s*Molecule:\s+(\w+)')
                    rerelcorr=re.compile(r'^\s*Relativistic Corrections:\s+(\w+ [\(\)\w\,]+)\s+')
                    recoretreat=re.compile(r'^\s*Core Treatment:\s+(\w+ \w+)\s+')
                    while f(reend,line) is None:
                        if self.lda is None:
                            self.lda=f(relda,line)
                        if self.gradcorr is None:
                            self.gradcorr=f(regradcorr,line)
                        if f(reunrestricted1,line):
                            while f(reunrestricted2,line) is None:
                                temp=f(reunrestricted3,line)
                                if temp is not None:
                                    if temp[:2]=='UN':
                                        self.unrestricted=1
                                    else:
                                        self.unrestricted=0
                                line=infile.readline()
                        if self.relcorr is None:
                            self.relcorr=f(rerelcorr,line)
                        if self.coretreat is None:
                            self.coretreat=f(recoretreat,line)
                        line=infile.readline()
                if f(respin,line):
                    reend=re.compile(r'^\s*\*+\s*$')
                    renumelec=re.compile(r'^\s*Total:\s*([\d\.]+)[a-zA-Z\(\)\-\s\+]+([\d\.]+)')
                    while f(reend,line) is None:
                        renumelec2=renumelec.search(line)
                        if renumelec2 is not None:
                            self.num_elecalpha=float(renumelec2.group(1))
                            self.num_elecbeta=float(renumelec2.group(2))
                            self.num_elec=self.num_elecalpha+self.num_elecbeta
                            self.totalspin=abs(self.num_elecalpha-self.num_elecbeta)*0.5*2+1 #numelectron*spinofelectron*2 +1
                        line=infile.readline()

                if f(regeometry,line):
                    reend=re.compile(r'\s*FRAGMENTS\s*$')
                    recoord=re.compile(r'^\s*\d+\s+(\w+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)')
                    while f(reend,line) is None:
                        temp=recoord.search(line)
                        if temp is not None:
                            self.atoms.append(temp.group(1))
                            tempcoord=[temp.group(2),temp.group(3),temp.group(4)]
                            self.coordinates.append([float(i) for i in tempcoord])
                        line=infile.readline()
                    self.geomcycle.append(self.coordinates)
                
                if f(regeomcycle,line):
                    reend1=re.compile(r'Number of elements of the density matrix')
                    reend2=re.compile(r'B O N D I N G   E N E R G Y')
                    reenergycycle=re.compile(r'current energy\s*([\d\.\-]+)')
                    regeom1=re.compile(r'\s*Coordinates ')
                    regeom2=re.compile(r'Number of elements of the density matrix')
                    regeom3=re.compile(r'^\s*\d+\s+(\w+)\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)')
                    while f(reend1,line) is None and f(reend2,line) is None: #the flag for the end should be better handled
                        if f(reenergycycle,line) is not None:
                            self.energycycle.append(float(f(reenergycycle,line))) #NOTE: UNITS!
                        if f(regeom1,line):
                            tempcoord=[]
                            #tempatom=[]
                            while f(regeom2,line) is None:
                                regeom4=regeom3.search(line) #first three are in bohr, next are in angstrom
                                if regeom4 is not None:
                                    tempcoord.append([regeom4.group(2),regeom4.group(3),regeom4.group(4)])
                                line=infile.readline()
                            self.geomcycle.append([[float(j) for j in i] for i in tempcoord])
                        else:
                            line=infile.readline()
                
                if f(resymmetry, line):
                    reend = re.compile(r'^\s*Electron Density at Nuclei\s*$')
                    resigma = re.compile(r'^\s*SIGMA\s*')
                    repi = re.compile(r'^\s*PI\s*')
                    redelta = re.compile(r'^\s*DELTA\s*')
                    rephi = re.compile(r'^\s*PHI\s*')
                    regeneral = re.compile(r'^\s*\w+\s*')
                    reocc = re.compile(r'^\s*(\w+)\s+\d+\s+(\w+)\s+([\d\.]+)\s*')
                    temp = {'A':{'sigma':[], 'pi':[], 'delta':[], 'phi':[]}, 'B':{'sigma':[], 'pi':[], 'delta':[], 'phi':[]}}
                    while f(reend, line) is None:
                        search = reocc.search(line)
                        if search is not None:
                            if search.group(1) == 'SIGMA':
                                temp[search.group(2)]['sigma'].append(float(search.group(3)))
                            if search.group(1) == 'PI':
                                temp[search.group(2)]['pi'].append(float(search.group(3)))
                            if search.group(1) == 'DELTA':
                                temp[search.group(2)]['delta'].append(float(search.group(3)))
                            if search.group(1) == 'PHI':
                                temp[search.group(2)]['phi'].append(float(search.group(3)))
                        line = infile.readline()
                    output = [[int(sum(temp[j][i])) for j in ['A','B']] for i in ['sigma','pi','delta','phi']]
                    self.occ = output

                if f(reenergy,line):
                    reend=re.compile(r'^\s*F R A G M E N T   E N E R G Y   T E R M S\s*')
                    reenergy=re.compile(r'^\s*Total Bonding Energy:\s+([\d\.\-\+eEdD]+)')
                    while f(reend,line) is None:
                        if f(reenergy,line) is not None:
                            self.energy_bond=float(f(reenergy,line))
                        line=infile.readline()

                if f(relog,line):
                    reend=re.compile(r'^\s*(NORMAL)?\s*TERMINATION\s*$')
                    rewarning=re.compile(r'WARNING:\s*(.+)$')
                    reerror=re.compile(r'ERROR:\s*(.+)$')
                    flag_normal = 0
                    while line != '':
                        rewarning2=re.search(rewarning,line)
                        if rewarning2 is not None:
                            self.warnings.append(rewarning2.group(1))
                        reerror2=re.search(reerror,line)
                        if reerror2 is not None:
                            self.errors.append(reerror2.group(1))
                        if reend is not None:
                            flag_normal = 1
                        line=infile.readline()
                    if flag_normal == 0:
                        self.errors.append('nonnormal termination')
                line=infile.readline()
        #problem: crashes for errors
        #if self.coordinates ==[]:
        #    raise AssertionError, 'bad data file'
'''
a=LogData(inpfile=sys.argv[1],program='adf')
print a.occ
print a.atoms
print a.coordinates
print   a.lda
print   a.gradcorr
print   a.relcorr
print   a.coretreat
print a.energy_bond,' bondenerg'
print a.num_elecalpha,'alp'
print a.num_elecbeta,'bet'
print a.num_elec,'el'
print a.totalspin,'S'
print a.geomcycle,'geomcycle',len(a.geomcycle)
print a.energycycle,'energycycle',len(a.energycycle)
'''
