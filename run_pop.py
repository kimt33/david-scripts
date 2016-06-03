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

def run_extensive_pop(fchkfile, hdf5file):
    ref_charges = IOData.from_file(fchkfile).numbers

    f = h5py.File('quasi_results.hdf5', 'r')
    g = h5py.File('pop_results.hdf5', 'w')

    quasitype_list = ['quasi', 'newquasi']
    aao_list = ['aambs', 'ano-rcc', 'minao', 'sto-6g']
    orbtype_list = ['quambo','quao','iao 1','iao 2', 'simple']
    for aao in aao_list:
        horton_data = wrapper_horton(fchkfile, aao+'.gbs')
        print aao
        for quasitype in quasitype_list:
            for orbtype in orbtype_list:
                pops = [np.zeros(ref_charges.size) for i in range(1+int(orbtype!='simple'))]
                print quasitype, orbtype
                for index, coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations, basis_map in zip(range(2), *horton_data[:6]):
                    if quasitype == 'quasi':
                        quasi = QuasiTransformation(coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations.astype(bool))
                    elif quasitype == 'newquasi':
                        quasi = NewQuasiTransformation(coeff_ab_mo, olp_ab_ab, olp_cao_ab, olp_cao_cao, occupations.astype(bool))
                    try:
                        coeff_ab_quasi = f['/'.join([aao, str(index), quasitype, orbtype, 'not_orth'])][:]
                    except KeyError:
                        print '/'.join([aao, str(index), quasitype, orbtype, 'not_orth'])
                        print f['/'.join([aao])].keys()
                        print f['/'.join([aao, str(index)])].keys()
                        print f['/'.join([aao, str(index), quasitype])].keys()
                        print f['/'.join([aao, str(index), quasitype, orbtype])].keys()
                    # orthogonalize
                    if orbtype == 'simple':
                        tmp_olp_ab_ab = quasi.olp_ab_ab
                        coeff_ab_omo = quasi.coeff_ab_omo
                        olp_quasi_quasi = coeff_ab_quasi.T.dot(tmp_olp_ab_ab).dot(coeff_ab_quasi)
                        olp_quasi_omo = coeff_ab_quasi.T.dot(tmp_olp_ab_ab).dot(coeff_ab_omo)
                        print np.linalg.matrix_rank(coeff_ab_quasi)
                        print np.linalg.matrix_rank(olp_quasi_quasi)
                        print 'xxxxxxxxxxxxxxx'
                        coeff_quasi_omo = project(olp_quasi_quasi, olp_quasi_omo)
                        print olp_quasi_omo.shape
                        print np.diag(olp_quasi_omo.T.dot(coeff_quasi_omo))
                        #print coeff_quasi_omo.T.dot(olp_quasi_quasi).dot(coeff_quasi_omo)
                        occ = occupations[occupations>0]
                        pops[0] += Mulliken(coeff_quasi_omo, occ, olp_quasi_quasi, len(set(basis_map)), basis_map).get_population()
                        continue
                    coeff_ab_oquasi = f['/'.join([aao, str(index), quasitype, orbtype, 'orth'])][:]
                    if orbtype == 'iao 2':
                        tmp_olp_ab_ab = np.zeros([quasi.num_ab+quasi.num_aao]*2)
                        tmp_olp_ab_ab[:quasi.num_ab] = np.hstack((quasi.olp_ab_ab, quasi.olp_ab_aao))
                        tmp_olp_ab_ab[quasi.num_ab:] = np.hstack((quasi.olp_aao_ab, quasi.olp_aao_aao))
                        coeff_ab_omo = np.vstack((quasi.coeff_ab_omo, np.zeros([quasi.num_aao, quasi.num_omo])))
                    else:
                        tmp_olp_ab_ab = quasi.olp_ab_ab
                        coeff_ab_omo = quasi.coeff_ab_omo
                    for i, coeff in enumerate([coeff_ab_quasi, coeff_ab_oquasi]):
                        olp_quasi_quasi = coeff.T.dot(tmp_olp_ab_ab).dot(coeff)
                        olp_quasi_omo = coeff.T.dot(tmp_olp_ab_ab).dot(coeff_ab_omo)
                        coeff_quasi_omo = project(olp_quasi_quasi, olp_quasi_omo)
                        occ = occupations[occupations>0]
                        pops[i] += Mulliken(coeff_quasi_omo, occ, olp_quasi_quasi, len(set(basis_map)), basis_map).get_population()
                print '/'.join([aao, str(index), quasitype, orbtype, 'not_orth'])
                print pops
                g.create_dataset('/'.join([aao, str(index), quasitype, orbtype, 'not_orth']), data=ref_charges-pops[0])
                if orbtype != 'simple':
                    g.create_dataset('/'.join([aao, str(index), quasitype, orbtype, 'orth']), data=ref_charges-pops[1])

import sys
run_extensive_pop(sys.argv[1], sys.argv[2])
