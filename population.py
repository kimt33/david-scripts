from horton.transformorbs.quasi import *
from horton.meanfield.project import ProjectExp, project_mol, measure_span
import horton.part.mulliken as mul
from horton.gbasis.lcgobasis import LCGOBasis
from horton.matrix.blockdiagmatrix import BlockDiagMatrix
from horton import context
import pickle
from os.path import splitext
import itertools as it
import numpy as np
from horton.log import log

log.set_level(2)

# TODO: probably smarter to use h5py instead of pickle (since i'm using dictionary of dictionaries anyways)
def run_calc(fchkfile):
    mol = Molecule.from_file(fchkfile)
    mol.orthogonalize()
    outname = splitext(fchkfile)[0]+'_pop.p'
    basis_list = ['aambs', 'anorcc', 'minao', 'sto-6g', '631g']
    quasi_list = ['quambo','quao','iao1','iao2']
    orth_list = ['orth', 'north']
    names = ['atomic basis', 'orth atomic basis']
    newbasis_flat = [LCGOBasis(mol.obasis), LCGOBasis(mol.obasis).orth_all]
    # TODO: NORMALIZE and INTRATOMICALLY ORTHOGONALIZE AB's!
    # transformations
    newbasis = {k:{j:{i:{h:None for h in orth_list} for i in ['quambo','quao','iao1','iao2','quasi']} for j in ['normal','special']} for k in basis_list}
    print 'making the basis'
    for obasis in basis_list:
        quasi = QuasiTransformation(mol, obasis2=obasis, minimal=True)
        names.append(obasis)
        newbasis_flat.append(LCGOBasis(quasi.obasis2).N)
        shell_map = quasi.obasis2.shell_map
        shell_types = quasi.obasis2.shell_types
        for quasitype in ['normal','special']:
            print obasis, quasitype
            if quasitype=='special':
                quasi = quasi.new_quasitransformation(num_oldvirtuals=100000)
            print 'not orthogonalized'
            print 'quambo',
            newbasis[obasis][quasitype]['quambo']['north'] = LCGOBasis(quasi.mol.obasis, quasi.quambo().T, index_map=shell_map, index_types=shell_types)
            print 'quao',
            newbasis[obasis][quasitype]['quao']['north'] = LCGOBasis(quasi.mol.obasis, quasi.quao().T, index_map=shell_map, index_types=shell_types)
            print 'iao1',
            newbasis[obasis][quasitype]['iao1']['north'] = LCGOBasis(quasi.mol.obasis, quasi.iao(iaotype=1).T, index_map=shell_map, index_types=shell_types)
            print 'iao2',
            newbasis[obasis][quasitype]['iao2']['north'] = LCGOBasis(quasi.mol.obasis, quasi.iao(iaotype=2).T, index_map=shell_map, index_types=shell_types)
            print 'orthogonalized'
            print 'quambo',
            newbasis[obasis][quasitype]['quambo']['orth'] = newbasis[obasis][quasitype]['quambo']['north'].orth_all
            print 'quao',
            newbasis[obasis][quasitype]['quao']['orth'] = newbasis[obasis][quasitype]['quao']['north'].orth_all
            print 'iao1',
            newbasis[obasis][quasitype]['iao1']['orth'] = newbasis[obasis][quasitype]['iao1']['north'].orth_all
            print 'iao2'
            newbasis[obasis][quasitype]['iao2']['orth'] = newbasis[obasis][quasitype]['iao2']['north'].orth_all
            newbasis[obasis][quasitype]['quasi'] = quasi
    # discard virtual mo's
    mol = mol.select_mo()
    # population
    print 'finding the populations'
    pop = mul.MullikenPopulation(mol)
    results_pop = {}
    results_pop['mulliken'] = {}
    results_pop['mulliken']['normal'] = pop.mulliken()
    pop = mul.MullikenPopulation(project_mol(mol, LCGOBasis(mol.obasis).orth_all))
    results_pop['lowdin'] = {}
    results_pop['lowdin']['normal'] = pop.mulliken()
    for i in ['quambo','quao','iao1','iao2']:
        results_pop['mulliken'][i] = {j:{k:None for k in basis_list} for j in ['normal', 'special']}
        results_pop['lowdin'][i] = {j:{k:None for k in basis_list} for j in ['normal', 'special']}
    for obasis in basis_list:
        for quasitype in ['normal', 'special']:
            for qao in quasi_list:
                for orth in orth_list:
                    print obasis, quasitype, qao, orth
                    # pop
                    quasi = newbasis[obasis][quasitype]['quasi']
                    # transform the basis
                    new_basis = newbasis[obasis][quasitype][qao][orth]
                    # project mol onto basis
                    quasi_newmol = project_mol(mol, new_basis)
                    # get pop
                    pop = mul.MullikenPopulation(quasi_newmol)
                    # weights set to 'temp' to prevent weights_mulliken function from running
                    #pop.weights = pop.weights_mulliken(quasi.obasis2)
                    # store population and name and basis
                    if orth=='north':
                        results_pop['mulliken'][qao][quasitype][obasis] = pop.mulliken()
                        names.append([obasis, quasitype, 'north', qao])
                        newbasis_flat.append(newbasis[obasis][quasitype][qao]['north'])
                    elif orth=='orth':
                        results_pop['lowdin'][qao][quasitype][obasis] = pop.mulliken()
                        names.append([obasis, quasitype, 'orth', qao])
                        newbasis_flat.append(newbasis[obasis][quasitype][qao]['orth'])
    # add occupied mo to newbasis_flat
    occ_mo_basis = mol.mo_basis
    names.append('occupied mo')
    newbasis_flat.append(occ_mo_basis)
    print 'calculating the overlaps'
    # overlap
    minmax = lambda M: np.min(np.amax(M,axis=1))
    results_overlap = np.zeros([len(names), len(names)])
    overlap_raw = []
    #results_overlap[:] = np.NAN
    for basis1, basis2 in it.combinations(enumerate(newbasis_flat),2):
        olp = ProjectExp(basis1[1], basis2[1], None).olp_12
        results_overlap[basis1[0], basis2[0]] = minmax(olp)
        overlap_raw.append(np.amax(olp, axis=1))
    print 'calculating the spans'
    # span
    results_span = np.zeros([len(names), len(names)])
    for basis1, basis2 in it.product(enumerate(newbasis_flat), repeat=2):
        # can basis 1 span basis 2
        results_span[basis1[0], basis2[0]] = measure_span(basis1[1], basis2[1])
    print 'dumping pickle'
    pickle.dump([results_pop, results_overlap, results_span, names, overlap_raw], open(outname, 'w'))
    with open(splitext(fchkfile)[0]+'.pop', 'w') as f:
        for i,j in results_pop.iteritems():
            for k,l in j.iteritems():
                if k=='normal':
                    f.write(i +' ' + k + ':\n')
                    f.write(str(l) +'\n')
                else:
                    for m,n in l.iteritems():
                        for p,q in n.iteritems():
                            f.write(i + ' ' + k + ' ' + m + ' ' + p + ':\n')
                            f.write(str(q) + '\n')
        f.write('\n')
        f.write(','.join([''.join(i) for i in names]))
        f.write('\n')
        for i in results_overlap:
            f.write(','.join([str(j) for j in i]))
            f.write('\n')

        f.write(','.join([''.join(i) for i in names]))
        f.write('\n')
        for i in results_span:
            f.write(','.join([str(j) for j in i]))
            f.write('\n')

