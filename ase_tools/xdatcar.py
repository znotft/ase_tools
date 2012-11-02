"""Functions for handling XDATCAR files

   read_xdatcar() returns an array of Atoms() objects with postions taken 
   from an XDATCAR file. 

   view_xdatcar() uses the ase gui to view the positions from an XDATCAR  
   file. 

   Author: Matthew Dyer    Date of last revision: 21 January 2009
"""
def read_xdatcar(filename='XDATCAR',poscarfile='CONTCAR'):
    """Import XDATCAR type file.

    Reads unitcell, atom positions and constraints from the POSCAR/CONTCAR
    file and tries to read atom types from POSCAR/CONTCAR header, if this fails
    the atom types are read from OUTCAR or POTCAR file.

    Then reads in positions from an XDATCAR file and returns an array of 
    Atoms() objects with the positions from this file.
    """
 
    import os
    from ase.io import vasp
    import numpy as np

    atoms=vasp.read_vasp(poscarfile)
    natoms = len(atoms)
   
    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename

    data=f.readlines()

    nimages=(len(data)-5)/(natoms+1)

    images = []
    for i in range(nimages):
        images.append(atoms.copy())
        pos=np.zeros((natoms,3),np.float)
        for j in range(natoms):
            pos[j,0]=float(data[6+i*(natoms+1)+j].split()[0])
            pos[j,1]=float(data[6+i*(natoms+1)+j].split()[1])
            pos[j,2]=float(data[6+i*(natoms+1)+j].split()[2])
        images[i].set_scaled_positions(pos)

    return images

def view_xdatcar(filename='XDATCAR',poscarfile='CONTCAR'):
    from ase import visualize
    visualize.view(read_xdatcar(filename, poscarfile))
