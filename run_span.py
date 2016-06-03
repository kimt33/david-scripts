import sys
import h5py
import numpy as np

from wrapper_horton import wrapper_horton
from quasi import QuasiTransformation, project
from newquasi import NewQuasiTransformation
import orthogonalization as orth
from mulliken import Mulliken
from horton.io.iodata import IOData
import orthogonalization as orth
from scipy.linalg import fractional_matrix_power

def run_extensive_span(fchkfile, hdf5file):
    pass
    f = h5py.File('quasi_results.hdf5', 'r')
    quasitype_list = ['quasi', 'newquasi']
    aao_list = ['aambs', 'ano-rcc', 'minao', 'sto-6g']
    orbtype_list = ['quambo','quao','iao 1','iao 2', 'simple']
    for aao in aao_list:
        horton_data = wrapper_horton(fchkfile, aao+'.gbs')
        print aao
        for quasitype in quasitype_list:
            for orbtype in orbtype_list:
                span=[]
                print quasitype, orbtype
                for index, coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations, basis_map in zip(range(2), *horton_data[:6]):
                    if quasitype == 'quasi':
                        quasi = QuasiTransformation(coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations.astype(bool))
                    elif quasitype == 'newquasi':
                        quasi = NewQuasiTransformation(coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations.astype(bool))

import sys
run_extensive_span(sys.argv[1], sys.argv[2])
