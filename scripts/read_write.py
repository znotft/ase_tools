#!/usr/bin/env python
"""
A script which averages a CHGCAR or LOCPOT file in one direction to make a 1D curve.
User must specify filename and direction on command line.
Depends on ase
"""

import sys
from ase.io import read, write

if len(sys.argv) != 3:
    raise IOError('Must specify input files on command line.\n For example read_write.py input.POSCAR output.geometry.in')

atoms = read(sys.argv[1])
atoms.write(sys.argv[2])
