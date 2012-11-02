#!/usr/bin/env python

"""Script which generates a trajectory file from a VASP calculation.

Usage:

In the working directory of the VASP calculation, do:

$ write_traj.py filename.traj

The directory will be written to the file ``filename.traj``. Note that
this script requires that the CONTCAR, XDATCAR and OUTCAR files from
the VASP calculation exist in the working directory.
"""


from ase.calculators.vasp import Vasp, xdat2traj
import sys

outputfile = sys.argv[1].lstrip()
xdat = xdat2traj(outputfile, calc = Vasp(restart=True))
xdat.convert()
