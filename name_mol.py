import compactsets
def chemicalformula(atoms):
    '''returns contracted chemical formula (atom symbol + number of atoms)
    
    atoms
        list of atomic symbols (strings)
        doesn't check if each symbol is a proper symbol
    '''
    atoms_set=compactsets.compact(atoms)
    output=''
    for i in atoms_set:
        output+=str(i[0])+str(i[2])
    return output


