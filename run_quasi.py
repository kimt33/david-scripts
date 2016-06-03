import sys
import h5py
import numpy as np

from wrapper_horton import wrapper_horton
from quasi import QuasiTransformation, project
from newquasi import NewQuasiTransformation
from mulliken import Mulliken
import orthogonalization as orth

def run_extensive_quasi(fchkfile):
    g = h5py.File('quasi_results.hdf5', 'w')

    quasitype_list = ['quasi', 'newquasi']
    aao_list = ['aambs', 'ano-rcc', 'minao', 'sto-6g']
    orbtype_list = ['quambo','quao','iao 1','iao 2', 'simple']
    for aao in aao_list:
        horton_data = wrapper_horton(fchkfile, aao+'.gbs')
        print aao
        for index, coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations, basis_map in zip(range(2), *horton_data[:6]):
            for quasitype in quasitype_list:
                if quasitype == 'quasi':
                    quasi = QuasiTransformation(coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations.astype(bool))
                elif quasitype == 'newquasi':
                    quasi = NewQuasiTransformation(coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations.astype(bool))
                for orbtype in orbtype_list:
                    print quasitype, orbtype
                    quasi_func = getattr(quasi, orbtype.split()[0])
                    if 'iao' in orbtype:
                        coeff_ab_quasi = quasi_func(iaotype=int(orbtype.split()[1]), is_orth=False)
                    else:
                        coeff_ab_quasi = quasi_func(is_orth=False)
                    g.create_dataset('/'.join([aao, str(index), quasitype, orbtype, 'not_orth']), data=coeff_ab_quasi)
                    # orthogonalize
                    if orbtype == 'simple':
                        tmp_olp_ab_ab = quasi.olp_ab_ab
                        coeff_ab_omo = quasi.coeff_ab_omo
                        olp_quasi_quasi = coeff_ab_quasi.T.dot(tmp_olp_ab_ab).dot(coeff_ab_quasi)
                        olp_quasi_omo = coeff_ab_quasi.T.dot(tmp_olp_ab_ab).dot(coeff_ab_omo)
                        coeff_quasi_omo = project(olp_quasi_quasi, olp_quasi_omo)
                        occ = occupations[occupations>0]
                        continue
                    if orbtype == 'iao 2':
                        tmp_olp_ab_ab = np.zeros([quasi.num_ab+quasi.num_aao]*2)
                        tmp_olp_ab_ab[:quasi.num_ab] = np.hstack((quasi.olp_ab_ab, quasi.olp_ab_aao))
                        tmp_olp_ab_ab[quasi.num_ab:] = np.hstack((quasi.olp_aao_ab, quasi.olp_aao_aao))
                    else:
                        tmp_olp_ab_ab = quasi.olp_ab_ab
                    olp_quasi_quasi = coeff_ab_quasi.T.dot(tmp_olp_ab_ab).dot(coeff_ab_quasi)
                    coeff_ab_oquasi = coeff_ab_quasi.dot(orth.power_symmetric(olp_quasi_quasi, -0.5))
                    olp_oquasi_oquasi = coeff_ab_oquasi.T.dot(tmp_olp_ab_ab).dot(coeff_ab_oquasi)
                    coeff_ab_oquasi *= np.diag(olp_oquasi_oquasi)**(-0.5)
                    g.create_dataset('/'.join([aao, str(index), quasitype, orbtype, 'orth']), data=coeff_ab_oquasi)

import sys
run_extensive_quasi(sys.argv[1])
