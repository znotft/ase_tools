"""Functions for handling OUTCAR files

   read_outcar() returns an array of Atoms() objects with postions taken 
   from an OUTCAR file. It also returns energies, forces, and the current 
   unit cell (but no stresses, yet).

   view_outcar() uses the ase gui to view the positions from an OUTCAR
   file. 

   Author: Felix Hanke    Date of last revision: 11 March 2010
"""
def read_outcar(filename='OUTCAR'):
    """Import OUTCAR type file.

    Reads unitcell, atom positions, energies, and forces from the OUTCAR file.

    CAREFUL: does not explicitly read constraints (yet?)
    """
    import os
    import numpy as np
    from ase.calculators import SinglePointCalculator
    from ase import Atoms, Atom

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename
    data    = f.readlines()
    natoms  = 0
    images  = []
    atoms   = Atoms(pbc = True)
    energy  = 0
    species = []
    species_num = []
    symbols = []
    for n,line in enumerate(data):
        if 'POSCAR:' in line:
            temp = line.split()
            species = temp[1:]
        if 'ions per type' in line:
            temp = line.split()
            for ispecies in range(len(species)):
                species_num += [int(temp[ispecies+4])]
                natoms += species_num[-1]
                for iatom in range(species_num[-1]): symbols += [species[ispecies]]
        if 'direct lattice vectors' in line:
            cell = []
            for i in range(3):
                temp = data[n+1+i].split()
                cell += [[float(temp[0]), float(temp[1]), float(temp[2])]]
            atoms.set_cell(cell)
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM' in line:
            energy = float(data[n+2].split()[4])
        if 'POSITION          ' in line:
            forces = []
            for iatom in range(natoms):
                temp    = data[n+2+iatom].split()
                atoms  += Atom(symbols[iatom],[float(temp[0]),float(temp[1]),float(temp[2])])
                forces += [[float(temp[3]),float(temp[4]),float(temp[5])]]
                atoms.set_calculator(SinglePointCalculator(energy,forces,None,None,atoms))
            images += [atoms]
            atoms = Atoms(pbc = True)
    return images

def view_outcar(filename='OUTCAR'):
    from ase import visualize
    visualize.view(read_outcar(filename))

