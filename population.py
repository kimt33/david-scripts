from horton.transformorbs.quasi import *
from horton.meanfield.project import ProjectExp, project_mol, measure_span
import horton.part.mulliken as mul
from horton.gbasis.lcgobasis import LCGOBasis
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
    newbasis_flat = [LCGOBasis(mol.obasis).N, LCGOBasis(mol.obasis).orth_all,]
    # TODO: NORMALIZE and INTRATOMICALLY ORTHOGONALIZE AB's!
    # transformations
    newbasis = {k:{j:{i:{h:None for h in orth_list} for i in ['quambo','quao','iao1','iao2','quasi']} for j in ['normal','special']} for k in basis_list} 
    for obasis in basis_list:
        quasi = QuasiTransformation(mol, obasis2=obasis, minimal=True)
        names.append(obasis)
        newbasis_flat.append(LCGOBasis(quasi.obasis2).N)
        for quasitype in ['normal','special']:
            if quasitype=='special':
                quasi = quasi.new_quasitransformation(num_oldvirtuals=0)
            newbasis[obasis][quasitype]['quambo']['north'] = LCGOBasis(quasi.mol.obasis, quasi.quambo().T)
            newbasis[obasis][quasitype]['quao']['north'] = LCGOBasis(quasi.mol.obasis, quasi.quao().T)
            newbasis[obasis][quasitype]['iao1']['north'] = LCGOBasis(quasi.mol.obasis, quasi.iao(iaotype=1).T)
            newbasis[obasis][quasitype]['iao2']['north'] = LCGOBasis(quasi.mol.obasis, quasi.iao(iaotype=2).T)
            newbasis[obasis][quasitype]['quambo']['orth'] = newbasis[obasis][quasitype]['quambo']['north'].orth_all
            newbasis[obasis][quasitype]['quao']['orth'] = newbasis[obasis][quasitype]['quao']['north'].orth_all
            newbasis[obasis][quasitype]['iao1']['orth'] = newbasis[obasis][quasitype]['iao1']['north'].orth_all
            newbasis[obasis][quasitype]['iao2']['orth'] = newbasis[obasis][quasitype]['iao2']['north'].orth_all
            newbasis[obasis][quasitype]['quasi'] = quasi
    # population
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
                    # pop
                    quasi = newbasis[obasis][quasitype]['quasi']
                    quasi_mol = quasi.mol
                    # transform the basis
                    new_basis = newbasis[obasis][quasitype][qao][orth]
                    # project mol onto basis
                    quasi_newmol = project_mol(quasi_mol, new_basis)
                    # get pop
                    pop = mul.MullikenPopulation(quasi_newmol, weights='temp')
                    # weights set to 'temp' to prevent weights_mulliken function from running
                    pop.weights = pop.weights_mulliken(quasi.obasis2)
                    # store population and name and basis
                    if orth=='north':
                        results_pop['mulliken'][qao][quasitype][obasis] = pop.mulliken()
                        names.append([obasis, quasitype, 'north', qao])
                        newbasis_flat.append(newbasis[obasis][quasitype][qao]['north'])
                    elif orth=='orth':
                        results_pop['lowdin'][qao][quasitype][obasis] = pop.mulliken()
                        names.append([obasis, quasitype, 'orth', qao])
                        newbasis_flat.append(newbasis[obasis][quasitype][qao]['orth'])
    # overlap
    minmax = lambda M: np.min(np.amax(M,axis=1))
    results_overlap = np.zeros([len(names), len(names)])
    #results_overlap[:] = np.NAN
    for basis1, basis2 in it.combinations(enumerate(newbasis_flat),2):
        results_overlap[basis1[0], basis2[0]] = minmax(ProjectExp(basis1[1], basis2[1], None).olp_12)
    # span
    results_span = np.zeros([len(names), len(names)])
    #results_span[:] = np.NAN
    for basis1, basis2 in it.combinations(enumerate(newbasis_flat),2):
        results_span[basis1[0], basis2[0]] = measure_span(basis1[1], basis2[1])
    pickle.dump([results_pop, results_overlap, results_span, names], open(outname, 'w'))
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

