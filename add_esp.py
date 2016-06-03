import sys
import h5py
from horton.io.iodata import IOData

def add_esp(fchkfile, hdf5file):
    npa_charges = IOData.from_file(fchkfile).npa_charges
    esp_charges = IOData.from_file(fchkfile).esp_charges
    mul_charges = IOData.from_file(fchkfile).mulliken_charges
    g = h5py.File('pop_results.hdf5', 'a')
    if 'npa' in g.keys():
        del g['npa']
    g.create_dataset('/npa', data=npa_charges)
    if 'esp' in g.keys():
        del g['esp']
    g.create_dataset('/esp', data=esp_charges)
    if 'mulliken' in g.keys():
        del g['mulliken']
    g.create_dataset('/mulliken', data=esp_charges)

import sys
add_esp(sys.argv[1], sys.argv[2])
