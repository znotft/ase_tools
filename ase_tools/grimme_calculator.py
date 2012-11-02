from ase.calculators.general import Calculator
from ase.io import write
from os import system
from numpy import zeros

class GrimmeCalculator(Calculator):
    def __init__(self, grimme_command, grimme_control, grimme_geometry=None, calc_internal=None, grimme_output=None, indices=None, mask=None, hirshfeld = None):
        self.calculator     = calc_internal
        self.grimme_command = grimme_command
        self.grimme_control = grimme_control
        self.hirshfeld      = hirshfeld
        if grimme_geometry is not None:
            self.grimme_geometry = grimme_geometry
        else:
            self.grimme_geometry = "geometry.GrimmeCalculator.in"
        if grimme_output is not None:
            self.grimme_output = grimme_output
        else:
            self.grimme_output = "results.grimme"
        self.indices = indices
        self.mask = mask
        self.atoms_internal = None
        self.old_atoms = None
    
    def update(self, atoms): 
        if ((self.atoms != atoms) or
            (self.atoms != self.old_atoms)
            ):
            self.calculate(atoms)

    def calculate(self, atoms):
        if self.indices is not None:
            indices = self.indices
        elif self.mask is not None:
            indices = []
            for i in range(len(atoms)):
                if self.mask[i]:
                    indices += [i]
        else:
            indices = range(len(atoms))
        if self.calculator is not None:
            self.atoms_internal = atoms.copy()
            self.atoms_internal.set_calculator(self.calculator)
            dft_energy = self.atoms_internal.get_potential_energy()
            dft_forces = self.atoms_internal.get_forces()
        else:
            dft_energy = 0
            dft_forces = zeros((len(atoms),3))
        atoms_temp = atoms.copy()
        atoms_temp.constraints = []
        atoms_temp = atoms_temp[indices]
        atoms_temp.write(self.grimme_geometry)
        if self.hirshfeld == None:
            command = self.grimme_command+' '+self.grimme_control+' '+self.grimme_geometry+' > '+self.grimme_output
        else:
            command = self.grimme_command+' '+self.grimme_control+' '+self.grimme_geometry+'  '+self.hirshfeld+' > '+self.grimme_output 
        exitcode = system(command)
        if exitcode != 0:
            raise RuntimeError('Grimme calculation exited with nonzero exitcode, please check!!!')
        lines = open(self.grimme_output).readlines()
        grimme_forces = zeros((len(atoms),3))
        for n, line in enumerate(lines):
            if line.rfind('| Total vdw energy') > -1:
                grimme_energy = float(line.split()[-2])
            if line.rfind('Total vdw forces') > -1:
                for iatom in range(len(indices)):
                    data = lines[n+iatom+1].split()
                    for iforce in range(3):
                        grimme_forces[indices[iatom], iforce] = float(data[1+iforce])
        self.old_atoms = atoms.copy()
        self.energy = dft_energy + grimme_energy
        self.forces = dft_forces + grimme_forces 
        
    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy
    
    def get_forces(self,atoms):
        self.update(atoms)
        return self.forces
