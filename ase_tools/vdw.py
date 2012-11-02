# Copyright (2009) Jonas Bjork
# email: j.bjork@liverpool.ac.uk
# web:   pcwww.liv.ac.uk/~bjork

import pickle
from ase.parallel import rank, barrier
import os
import sys
from random import random
import time

class vdw_VASP:

    """This module defines an ASE interface to Andris Gulans's implementation
    of Langreth-Lundqvist's vdW-DF [Phys. Rev. Lett. 92, 246401 (2004)].

    To run the program the following environmental parameters MUST be set:

    .. glossary::
    
      VASP
        Points to the VASP binary to be used.

      VDW_PATH
        Points to the directory with all data files needed to run
        the vdw_nonsc program.

      VDW_NONSC
        Contains information on how the vdw_nonsc program should
        be executed. For example::
          
          $ export VDW_NONSC='mpprun vdw-nonsc-mpi'
    
    Additionally the following environmental parameter can be set:

    .. glossary::

      VASP_NONSC
        Points to additional VASP binary to be used for an extra
        VASP calculation which for example can be used to evaluate
        the full electron charge density, or calculate the
        E = E_DFT - E_c,GGA + E_c,LDA
        energy. The energy calculated in this optional step will
        be stored as 'E_2' in the pickle file. 
        Note that this step is always done NON SELF-CONSISTENTLY.

    """

    def __init__(self, name, atoms, nonsc_calculation=False):
        self.name = name
        self.atoms = atoms
        self.calc_name = atoms.get_calculator().__module__.split('.')[-1]
        if self.calc_name == 'aims':
            raise NotImplementedError('FHI-aims is not currently supported, but hopefully soon.')
        elif self.calc_name not in ['vasp']:
            raise NotImplementedError('Use a supported calculator.')
        #if not os.environ['VASP_COMMAND']:
        try:
            self.vasp = os.environ['VASP_COMMAND']
        except:
            raise KeyError('Environmental variable VASP_COMMAND, containing information on how to execute VASP, has to be set. If you are using the environmental variable VASP_SCRIPT, please change to VASP_COMMAND. For example:\n export VASP_COMMAND=\'mpirun vasp\'.')
        try:
            self.vdw_path = os.environ['VDW_PATH']
        except KeyError:
            raise KeyError('Environmental variable VDW_PATH, pointing to the directory with the vdW data files, has to be set.')
        try:
            self.vdw_nonsc = os.environ['VDW_NONSC']
        except KeyError:
            raise KeyError('Environmental variable VDW_NONSC, containing information on how to execute the vdW-DF program, has to be set.')
        if nonsc_calculation:
            try:
                self.vasp_nonsc = os.environ['VASP_NONSC']
            except:
                raise KeyError('To perform a second (non self-consistent) VASP calculation the environmental variable VASP_NONSC, containing information on how to execute vasp, has to be set.')
            self.nonsc_calculation = True
        else:
            self.nonsc_calculation = False

    def prepare(self, nonsc_calculation=None):
        """Prepare for the vdW calculation.

        The method runs the necessary VASP calculations and generates
        the necessary input files for the vdW calculation."""

        if os.path.isfile(self.name + '.pckl') or os.path.isfile(self.name + '_no-vdw.pckl'):
            self.t1 = self.t2 = int(time.time())
            return
        fd = open(self.name + '_no-vdw.pckl', 'w')
        calc = self.atoms.get_calculator()
        calc.output_template = self.name + '.' + self.calc_name
        if self.calc_name == 'aims':
            calc.set(output=['cube total_density', 
                             'cube origin 7.5 5.0 5.0'])
        sys.stdout.write('Running self-consistent DFT calculation... ')
        self.t1 = int(time.time())
        E_dft = self.atoms.get_potential_energy()
        sys.stdout.write('Done!\n\n')
        sys.stdout.flush()
        self.t2 = int(time.time())
        if nonsc_calculation is None:
            nonsc_calculation = self.nonsc_calculation
        if nonsc_calculation and self.calc_name == 'vasp':
            calc.set(istart=1,
                     icharg=11,
                     nelm=1)
            calc.output_template = self.name + '.vasp-nonsc'
            # Change command for execution of VASP
            os.environ['VASP_COMMAND'] = self.vasp_nonsc
            sys.stdout.write('Running non self-consistent DFT calculation... ')
            E_2 = self.atoms.get_potential_energy()
            sys.stdout.write('Done!\n\n')
            self.t_nonsc = time.time()
            sys.stdout.flush()
            # Change back to default command for execution of VASP
            os.environ['VASP_COMMAND'] = self.vasp
        else:
            E_2 = None
        if rank == 0:
            self.write_atoms_file()
            sys.stdout.write('Converting charge density... ')
            sys.stdout.flush()
            self.write_chg_file()
            sys.stdout.write(' Done!\n\n')
            sys.stdout.flush()
            self.t3 = int(time.time())
            if self.calc_name == 'vasp':
                os.system('cp OUTCAR %s.OUTCAR' % self.name)
            pickle.dump([{'E_dft': E_dft,
                          'E_2': E_2,
                          'E_LDA-c': None,
                          'E_LDA-x': None,
                          'E_PBE-c': None,
                          'E_PBE-x': None,
                          'E_revPBE-x': None,
                          'E_nl': None,
                          'E_vdW': None,
                          'atoms': self.atoms,
                          }],
                        fd)
            fd.close()

    def run_vdw(self, nonsc_calculation=None):
        """Starts the vdW calculation.

        The method runs only if the pickle corresponding pickle
        file does not exist.

        The vdw_nonsc program is executed in the directory set
        by the environmental variable VDW_PATH where the kernel
        and other data files must exist."""

        if os.path.isfile('%s.pckl' % self.name):
            return
        self.prepare(nonsc_calculation)
        sys.stdout.write('Running non self-consistent vdW-DF calculation... ')
        sys.stdout.flush()
        barrier()
        path = os.path.abspath('.')
        fd = open(path + '/%s.pckl' % self.name, 'w')
        data = pickle.load(open('%s_no-vdw.pckl' % self.name, 'r'))[0]
        E_dft = data['E_dft']        
        # Change path to where the vdw data files are
        os.chdir(self.vdw_path)
        # Execute vdw vdw_nonsc program
        tmp_name = path + '/%s' % self.name
        exitcode = os.system('%s %s.atoms.inp %s.chargeden > %s.vdw.out' % (self.vdw_nonsc, tmp_name, tmp_name, tmp_name))
        os.chdir(path)
        if exitcode != 0:
            raise RuntimeError('vdw_nonsc exited with exit code: %d.  ' % exitcode)
        barrier()
        sys.stdout.flush()
        # Extract and pickle the output
        if rank == 0:
            [E_LDAc, E_LDAx, E_PBEc, E_PBEx, E_revPBEx, E_nl] = self.extract()
            pickle.dump([{'E_dft': E_dft,
                          'E_2': data['E_2'],
                          'E_LDA-c': E_LDAc,
                          'E_LDA-x': E_LDAx,
                          'E_PBE-c': E_PBEc,
                          'E_PBE-x': E_PBEx,
                          'E_revPBE-x': E_revPBEx,
                          'E_nl': E_nl,
                          'E_vdW': E_dft-E_PBEc+E_LDAc+E_nl,
                          'atoms': data['atoms'],
                          }],
                        fd)
            fd.close()
            # Delete input files
            self.clean()

        sys.stdout.write('Done!\n\n')
        self.tf = int(time.time())
        print '================================================================================'
        print 'Finished vdW-DF calculation \'%s\' on %s.\n' % (self.name, time.ctime())
        print 'Timings:\n'
        print 'Self-consistent DFT calculation:        %3d hours %2d minutes %2d seconds' % \
            ((self.t2-self.t1)/3600, ((self.t2-self.t1)/60)%60, (self.t2-self.t1)%60)
        if self.nonsc_calculation:
            print 'Non self-consistent DFT calculation:    %3d hours %2d minutes %2d seconds' % \
                ((self.t_nonsc-self.t2)/3600, ((self.t_nonsc-self.t2)/60)%60, (self.t_nonsc-self.t2)%60)
            self.t2 = self.t_nonsc
        print 'Converting charge density:              %3d hours %2d minutes %2d seconds' % \
            ((self.t3-self.t2)/3600, ((self.t3-self.t2)/60)%60, (self.t3-self.t2)%60)
        print 'Non-self-consistent vdW-DF calculation: %3d hours %2d minutes %2d seconds' % \
            ((self.tf-self.t3)/3600, ((self.tf-self.t2)/60)%60, (self.tf-self.t3)%60)
        print '---------------------------------------------------------------------------'
        print 'Total time: %3d hours %2d minutes %2d seconds' % ((self.tf-self.t1)/3600, ((self.tf-self.t1)/60)%60, (self.tf-self.t1)%60)
        print '================================================================================\n'
        return

    def clean(self):
        """Deletes input files when vdW calculation has finished.

        The method deletes the charge file, the atoms positions file
        and the pickle file with the results from before the vdW
        calculation."""

        os.remove('%s_no-vdw.pckl' % self.name)
        os.remove('%s.chargeden' % self.name)
        os.remove('%s.atoms.inp' % self.name)

    def extract(self, filename=None):
        """Extracts data from the output from the vdW calculation.

        The method reads all the different (non-selfconsistent) energy
        contributions."""

        if filename is None:
            filename = '%s.vdw.out' % self.name
        for line in open(filename, 'r'):
            if line.rfind('LDA correlation energy') > -1:
                E_LDAc = float(line.split()[-2])
            elif line.rfind('LDA exchange energy') > -1:
                E_LDAx = float(line.split()[-2])
            elif line.rfind('PBE correlation energy') > -1:
                E_PBEc = float(line.split()[-2])
            elif line.rfind('revPBE exchange energy') > -1:
                E_revPBEx = float(line.split()[-2])
            elif line.rfind('PBE exchange energy') > -1:
                E_PBEx = float(line.split()[-2])
            elif line.rfind('Non-local correlation energy') > -1:
                E_nl = float(line.split()[-2])
        return [E_LDAc, E_LDAx, E_PBEc, E_PBEx, E_revPBEx, E_nl]

    def read(self):
        """Loads the data from the pickle output.

        .. highlight:: python

        To read a quantity from the pickle file do::
        
            >>> read()['quantity']

        The following quantities can be extracted from the pickle file

        .. glossary::

          E_dft
            Self-consistent DFT energy from the first
            VASP calculation.

          E_2
            Non self-consistent DFT energy from the second
            VASP calculation.

          E_LDA-c
            Non self-consistent LDA correlation energy.

          E_LDA-x
            Non self-consistent LDA exchange energy.

          E_PBE-c
            Non self-consistent (rev)PBE correlation energy.

          E_PBE-x
            Non self-consistent revPBE exchange energy.

          E_revPBE-x
            Non self-consistent revPBE exchange energy.

          E_nl
            Non self-consistent non-local correlation energy.

          E_vdW
            E_dft - E_PBE-c + E_LDA-c + E_nl

          atoms
            ASE Atoms object. From the Atoms object also information about
            the calculator can be extracted, and hence all the settings used
            for a calculation can be extracted at any point although the VASP
            input files are not stored.

        """

        self.run_vdw()
        return pickle.load(open(self.name + '.pckl', 'r'))[0]

    def get_total_energy(self):
        """Returns the total energy (including the vdW-correction).

        The method returns
        E_DFT - E_c,GGA + E_c,LDA + E_c,nl,
        thus, it uses the energy from the first VASP calculation and the
        self-consistent energies calculated by the vdw_nonsc program."""

        return self.read()['E_vdW']

    def get_dft_energy(self):
        """Returns the self-consistent DFT total energy from the first VASP calculation."""

        return self.read()['E_dft']

    def write_atoms_file(self):
        """Writes the atoms input file for the vdW calculation."""

        fd = open('%s.atoms.inp' % self.name, 'w')
        fd.write('%d\n' % len(self.atoms))
        for atom in self.atoms:
            pos = atom.position
            fd.write(' 18 %15.10f %15.10f %15.10f\n' % (pos[0], pos[1], pos[2]))
        fd.close()

    def write_chg_file(self):
        """Write the charge input file for the vdW calculation."""

        fd = open('%s.chargeden' % self.name, 'w')
        fd.write('lattice\n')
        cell = self.atoms.get_cell()
        for n in range(3):
            fd.write(' %14.7f %14.7f %14.7f\n' % (cell[n,0], cell[n,1], cell[n,2]))
        fd.write('Mesh\n')
        if self.calc_name == 'vasp':
            from ase.calculators.vasp import VaspChargeDensity
            charge = VaspChargeDensity('CHGCAR')
        elif self.calc_name == 'aims':
            from ase.io.cube import read_cube
            charge = read_cube('total_density.cube', read_data=True)
        a,b,c = charge.chg[0].shape
        charge = charge.chg[0].T.ravel()
        charge *= self.atoms.get_volume()
        fd.write(' %5d %5d %5d\n' % (a,b,c))
        fd.write('Density\n')
        for c in charge:
            fd.write('%17.15E\n' % c)
        fd.close()
