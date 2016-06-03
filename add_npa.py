import sys
import h5py
from horton.io.iodata import IOData

def add_npa(fchkfile, hdf5file):
    npa_charges = IOData.from_file(fchkfile).npa_charges
    g = h5py.File('pop_results.hdf5', 'a')
    if 'npa' in g.keys():
        del g['npa']
    g.create_dataset('/npa', data=npa_charges)

import sys
add_npa(sys.argv[1], sys.argv[2])
