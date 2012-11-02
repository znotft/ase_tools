#!/usr/bin/env python
"""
A script which reads POSCAR files, repeat 
the atoms and writes a new POSCAR file.

Usage:
$ repeat.py filename n1 n2 n3
where n1, n2 and n3 specifies how many times to repeat 
the unit cell in the direction of the first, second and
third unit cell vector, respectively.

To repeat only along the two first unit cell vectors it 
is sufficient to do
$ repeat.py filename n1 n2

and only along the first unit cell vector
$ repeat.py filename n1
"""


from ase.io.vasp import read_vasp, write_vasp
import sys

nargs = len(sys.argv)
if not nargs in [3, 4, 5]:
    raise SyntaxError('Must specify name of input file and \
direction on command line, e.g. \n \
$ repeat.py POSCAR 2 2 2')

file = sys.argv[1]
n1 = int(sys.argv[2])
if nargs > 3:
    n2 = int(sys.argv[3])
    if nargs > 4:
        n3 = int(sys.argv[4])
    else:
        n3 = 1
else:
    n2 = n3 = 1

atoms = read_vasp(file)
atoms.constraints = []
atoms = atoms.repeat([n1, n2, n3])
write_vasp(file+'-%dx%dx%d' % (n1, n2, n3), atoms)
