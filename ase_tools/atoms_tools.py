def closest_atom(atoms, atom_number, atom_type):
    """Method to find the nearest neighbour of a certain kind to an atom.

    The method returns the index of the nearest neighbour.

    atoms: ASE atoms object.
        The atoms object to work with.
    atoms_number: int
        Index of the atom of interest.
    atoms_type: str
        Atomic species of nearest neighbour to find.
    """

    d = None
    nb = None
    symbols = atoms.get_chemical_symbols()
    for n in range(len(atoms)):
        if n == atom_number:
            continue
        if atom_type == symbols[n]:
            d_test = atoms.get_distance(atom_number, n)
            if not d or atoms.get_distance(atom_number, n) < d:
                d = d_test
                nb = n
    return nb
