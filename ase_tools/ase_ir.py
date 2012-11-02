# -*- coding: utf-8 -*-

"""Infrared spectroscopy"""

import pickle
from math import sin, pi, sqrt, exp, log
from os import remove
from os.path import isfile

import numpy as np

import ase.units as units
from ase.io.trajectory import PickleTrajectory
from ase.parallel import rank, barrier

# Atomic mass unit to rest mass of electron
amutome = 1822.888485
# Conversion factors from atomic units to (D/Angstrom)^2/amu
# and electron x Angstrom to Debye
conv = 42055.44331
qeangtod = 4.803204272

class InfraRed:
    """Class for calculating vibrational modes and infrared intensities
    using finite difference.

    The vibrational modes are calculated from a finite difference
    approximation of the Hessian matrix and the IR intensities from
    a finite difference approximation of the gradient of the dipole
    momen. The method is described in:

      D. Porezag, M. R. Peterson:
      "Infrared intensities and Raman-scattering activities within
      density-functional theory",
      Phys. Rev. B 54, 7830 (1996)

    The calculator object (calc) linked to the Atoms object (atoms) must have the 
    attribute: calc.get_dipole_moment(atoms)

    The *summary*, *get_energies()* and *get_frequencies()*
    methods all take an optional *method* keyword.  Use
    method='Frederiksen' to use the method described in:

      T. Frederiksen, M. Paulsson, M. Brandbyge, A. P. Jauho:
      "Inelastic transport theory from first-principles: methodology
      and applications for nanoscale devices", 
      Phys. Rev. B 75, 205413 (2007) 

    atoms: Atoms object
        The atoms to work on.
    indices: list of int
        List of indices of atoms to vibrate.  Default behavior is
        to vibrate all atoms.
    name: str
        Name to use for files.
    delta: float
        Magnitude of displacements.
    nfree: int
        Number of displacements per degree of freedom.

    """
    def __init__(self, atoms, indices=None, name='ir', delta=0.01, nfree=2, idipol=4):
	if int(nfree) not in [2, 4]:
            raise ValueError('Only %i and %i supported for nfree.' % (2, 4))
        self.atoms = atoms
        if atoms.constraints:
            print "WARNING! \n Your Atoms object is constrained. Some forces may be unintended set to zero. \n"
        self.calc = atoms.get_calculator()
        if indices is None:
            indices = range(len(atoms))
        self.indices = np.asarray(indices)
        self.nfree = nfree/2
        self.ndof = 3*len(self.indices)
        self.name = name+'-d%.3f' % delta
        self.delta = delta
        self.H = None
        if idipol not in range(1,5):
            raise ValueError('"idipol" must have the value 1, 2, 3 or 4.')
        self.idipol = idipol

    def run(self):
        """Run the vibration calculations.

        This will calculate the forces for 6 displacements per atom
        ±x, ±y, ±z.  Only those calculations that are not already done
        will be started. Be aware that an interrupted calculation may
        produce an empty file (ending with .pckl), which must be deleted
        before restarting the job. Otherwise the forces will not be
        calculated for that displacement."""

        p = self.atoms.positions.copy()
        # Get forces and dipole moment at zero displacement
        filename = '%s.static.pckl' % self.name
        if not isfile(filename):
            forces = self.atoms.get_forces()
            dipole = self.calc.get_dipole_moment(self.atoms)
            if rank == 0:
                fd = open(filename, 'w')
                pickle.dump([forces, dipole], fd)
                fd.close()
        # Get forces and dipole moments for the finite displacements
        for a in self.indices:
            for i in range(3):
                for sign in [-1, 1]:
                    for n in range(1, self.nfree+1):
                        filename = '%s.%d%s%s.pckl' % (self.name, a,
                                                   'xyz'[i], n*' +-'[sign])
                        if isfile(filename):
                            continue
                        barrier()
                        self.atoms.positions[a, i] = p[a, i]+sign*n*self.delta
                        forces = self.atoms.get_forces()
                        dipole = self.calc.get_dipole_moment(self.atoms)
                        if rank == 0:
                            fd = open(filename, 'w')
                            pickle.dump([forces, dipole], fd)
                            fd.close()
                        self.atoms.positions[a, i] = p[a, i]
        self.atoms.set_positions(p)

    def clean(self):
        """Remove pckl-files."""

        name = '%s.static.pckl' % self.name
        if isfile(name):
            remove(name)
        for a in self.indices:
            for i in 'xyz':
                for sign in '-+':
                    name = '%s.%d%s%s.pckl' % (self.name, a, i, sign)
                    if isfile(name):
                        remove(name)

    def read(self, method='standard'):
        """Read files and calculate IR spectrum."""
        # Get "static" dipole moment and forces
        name = '%s.static.pckl' % self.name
        [forces_zero, dipole_zero] = pickle.load(open(name))
        self.dipole_zero = (sum(dipole_zero**2)**0.5)*qeangtod
        self.force_zero = max([sum((forces_zero[j])**2)**0.5 for j in self.indices])

        self.method = method.lower()
        assert self.method in ['standard', 'frederiksen']
        H = np.empty((self.ndof, self.ndof))
        dpdx = np.empty((self.ndof, 3))
        r = 0
        for a in self.indices:
            for i in 'xyz':
                name = '%s.%d%s' % (self.name, a, i)
                [fminus, dminus] = pickle.load(open(name + '-.pckl'))
                [fplus, dplus] = pickle.load(open(name + '+.pckl'))

                if self.nfree == 2:
                    [fminusminus, dminusminus] = pickle.load(open(name + '--.pckl'))
                    [fplusplus, dplusplus] = pickle.load(open(name + '++.pckl'))
                if self.method == 'frederiksen':
                    fminus[a] += -fminus.sum(0)
                    fplus[a] += -fplus.sum(0)
                if self.nfree == 1:
                    H[r] = (fminus - fplus)[self.indices].ravel() / (4 * self.delta)
                    dpdx[r] = -(dplus - dminus)/(2.*self.delta)
                if self.nfree == 2:
                    H[r] = (-fminusminus+8*fminus-8*fplus+fplusplus)[self.indices].ravel()/(24 * self.delta)
                    dpdx[r] = (-dplusplus + 8*dplus - 8*dminus +dminusminus) / (12*self.delta)
                if self.idipol != 4:
                    dpdx[r][(self.idipol-1-1)%3] = 0
                    dpdx[r][(self.idipol-1-2)%3] = 0
                r += 1
        # Calculate eigenfrequencies and eigenvectors
        m = self.atoms.get_masses()       
        H += H.copy().T
        self.H = H
        m = self.atoms.get_masses()
        self.im = np.repeat(m[self.indices]**-0.5, 3)
        omega2, modes = np.linalg.eigh(self.im[:, None] * H * self.im)
        self.modes = modes.T.copy()

        # Calculate intensities
        dpdq = np.array([dpdx[j]/sqrt(m[self.indices[j/3]]*amutome) for j in range(self.ndof)])
        dpdQ = np.dot(dpdq.T, modes)
        dpdQ = dpdQ.T
        intensities = np.array([sum(dpdQ[j]**2) for j in range(self.ndof)])
        # Conversion factor:
        s = units._hbar * 1e10 / sqrt(units._e * units._amu)
        self.hnu = s * omega2.astype(complex)**0.5
        self.intensities = intensities*conv

    def get_energies(self, method='standard'):
        """Get vibration energies in eV."""
        if self.H is None or method.lower() != self.method:
            self.read(method)
        return self.hnu

    def get_frequencies(self, method='standard'):
        """Get vibration frequencies in cm^-1."""
        s = 0.01 * units._e / units._c / units._hplanck
        return s * self.get_energies(method)

    def summary(self, method='standard'):
        """Print summary of eigenfrequencies and corresponding IR intensities."""

        hnu = self.get_energies(method)
        s = 0.01 * units._e / units._c / units._hplanck
        print '-------------------------------------'
        print ' Mode    Frequency        Intensity'
        print '  #    meV     cm^-1   (D/Å)^2 amu^-1'
        print '-------------------------------------'
        for n, e in enumerate(hnu):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
            print '%3d %6.1f%s  %7.1f%s  %9.4f' % (n, 1000 * e, c, s * e, c, self.intensities[n])
        print '-------------------------------------'
        print 'Zero-point energy: %.3f eV' % self.get_zero_point_energy()
        print 'Static dipole moment: %.3f D' % self.dipole_zero
        print 'Maximum "equilibrium" force on atom: %.4f eV/Å' % self.force_zero
        print

    def get_zero_point_energy(self):
        """Returns the zero point energy."""
        return 0.5 * self.hnu.real.sum()

    def get_mode(self, n):
        """Returns the eigenvector of the nth mode."""
        mode = np.zeros((len(self.atoms), 3))
        mode[self.indices] = (self.modes[n] * self.im).reshape((-1, 3))
        return mode

    def write_mode(self, n, kT=units.kB * 300, nimages=30):
        """Write mode n to trajectory file."""
        mode = self.get_mode(n) * sqrt(kT / self.hnu[n])
        p = self.atoms.positions.copy()
        n %= 3 * len(self.indices)
        traj = PickleTrajectory('%s.%d.traj' % (self.name, n), 'w')
        calc = self.atoms.get_calculator()
        self.atoms.set_calculator()
        for x in np.linspace(0, 2 * pi, nimages, endpoint=False):
            self.atoms.set_positions(p + sin(x) * mode)
            traj.write(self.atoms)
        self.atoms.set_positions(p)
        self.atoms.set_calculator(calc)
        traj.close()

    def write_spectra(self, out='ir-spectra.dat', start=800, end=4000, npts=None, width=4, type='Gaussian'):
        """Write out infrared spectrum to file.

        Start and end point, and width of the Gaussian/Lorentzian should be given in cm^-1."""
        if not npts:
            npts = (end-start)*width*5+1
        method='standard'
        frequencies = self.get_frequencies(method).real
        intensities=self.intensities
        if type.lower()=='lorentzian':
            lineshape = 1
            intensities = intensities*width*pi/2.
        elif type.lower()=='gaussian':
            lineshape = 0
            sigma = width/2./sqrt(2.*log(2.))
        else:
            raise ValueError('%s not a valid type of lineshape' % type)
        #Make array with spectrum data
        spectrum=np.zeros(npts,np.float)
        ediff = (end-start)/float(npts-1)
        for i in range(1,npts):
            energy = end - float(i)*ediff
            for j in range(len(frequencies)):
                if lineshape:
                    spectrum[i] = spectrum[i]+intensities[j]*0.5*width/pi/((energy-frequencies[j])**2+0.25*width**2)
                else:
                    spectrum[i] = spectrum[i]+intensities[j]*exp(-(energy-frequencies[j])**2/2./sigma**2)
        #Write out spectrum in file. First column is just intensities. 
        #Second column is absorbance scaled so that data runs from 1 to 0
        spectrumfile = open(out,"w")
        for i in range(1,npts):
            energy = end - float(i)*ediff
            spectrumfile.write("%f %15.5e %15.5e\n" % (energy,spectrum[i],1.-spectrum[i]/spectrum.max()))
        spectrumfile.close()
