#!/usr/bin/env python

from distutils.core import setup

setup(name='ase tools',
      version='0.2',
      description='Tools using the Atomic Simulation Environment.',
      author='Jonas Bjork, Matthew Dyer',
      author_email='jonbj@ifm.liu.se, msd30@liverpool.ac.uk',
      license='GNU GPL',
      packages=['ase_tools'],
      scripts = ['scripts/chgadd.py',
                 'scripts/chgdiff.py', 
                 'scripts/curvefit.py',
                 'scripts/read_write.py',
                 'scripts/repeat.py',
                 'scripts/slice.py', 
                 'scripts/splitparchg.py',
                 'scripts/vtotav.py',
                 'scripts/write_traj.py',]
      )

