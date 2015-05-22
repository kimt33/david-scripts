from horton.io.smartreader import *

class XYZReader:
    """reads xyz file

    """
    def __init__(self, inpfile):
        """ 

        Args:
            inpfile: string that describes the location of xyz file
        """
        self.filename=inpfile
        self.title=''
        self.num_atoms=0
        self.atoms=[]
        self.coordinates=[]
        extract_xyz()

    def extract_xyz(self):
        """ extracts the informations from self.filename

        """
        with open(self.filename, 'r') as f:
            
