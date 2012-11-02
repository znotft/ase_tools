from ase.lattice.cubic import *
from ase.lattice.tetragonal import *

# Define fluorite structure
class FluoriteFactory(SimpleCubicFactory):
    "A factory for creating fluorite (CaF2) lattices"
    xtal_name = 'fluorite'
    bravais_basis = [[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],
                     [0.25,0.25,0.25],[0.25,0.75,0.25,],[0.75,0.25,0.25],[0.75,0.75,0.25],
                     [0.25,0.25,0.75],[0.25,0.75,0.75,],[0.75,0.25,0.75],[0.75,0.75,0.75]]
    element_basis = (0,0,0,0,1,1,1,1,1,1,1,1)

# Define bixbyite structure
class BixbyiteFactory(BodyCenteredCubicFactory):
    "A factory for creating bixbyite (Y2O3) lattices"
    def __init__(self,u=-0.0327,x=0.3905,y=0.1518,z=0.3803):
        self.u=u
        self.x=x
        self.y=y
        self.z=z
    # Define parameters needed in class
        self.xtal_name = 'bixbyite'
        self.bravais_basis = [[1./4.,1./4.,1./4.],[1./4.,3./4.,3./4.],[3./4.,3./4.,1./4.],[3./4.,1./4.,3./4.],
                             [self.u,0,1./4.],[-self.u+1./2.,0,3./4.],[1./4.,self.u,0],[3./4.,-self.u+1./2.,0],
                             [0,1./4.,self.u],[0,3./4.,-self.u+1./2.],[-self.u,0,3./4.],[self.u+1./2.,0,1./4.],
                             [3./4.,-self.u,0],[1./4.,self.u+1./2.,0],[0,3./4.,-self.u],[0,1./4.,self.u+1./2.],
                             [self.x,self.y,self.z],[-self.x+1./2.,-self.y,self.z+1./2.],[-self.x,self.y+1./2.,-self.z+1./2.],[self.x+1./2.,-self.y+1./2.,-self.z],
                             [self.z,self.x,self.y],[self.z+1./2.,-self.x+1./2.,-self.y],[-self.z+1./2.,-self.x,self.y+1./2.],[-self.z,self.x+1./2.,-self.y+1./2.],
                             [self.y,self.z,self.x],[-self.y,self.z+1./2.,-self.x+1./2.],[self.y+1./2.,-self.z+1./2.,-self.x],[-self.y+1./2.,-self.z,self.x+1./2.],
                             [-self.x,-self.y,-self.z],[self.x+1./2.,self.y,-self.z+1./2.],[self.x,-self.y+1./2.,self.z+1./2.],[-self.x+1./2.,self.y+1./2.,self.z],
                             [-self.z,-self.x,-self.y],[-self.z+1./2.,self.x+1./2.,self.y],[self.z+1./2.,self.x,-self.y+1./2.],[self.z,-self.x+1./2.,self.y+1./2.],
                             [-self.y,-self.z,-self.x],[self.y,-self.z+1./2.,self.x+1./2.],[-self.y+1./2.,self.z+1./2.,self.x],[self.y+1./2.,self.z,-self.x+1./2.]]
        self.element_basis = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# Define cubic perovskite
class PerovskiteCubicFactory(SimpleCubicFactory):
    "A factory for creating perovskite (ABO3) lattices"
    xtal_name = 'cubic perovskite'
    bravais_basis = [[0.5,0.5,0.5],
                     [0.,0.,0.],
                     [0.5,0.,0.],[0.,0.5,0.],[0.,0.,0.5]]
    element_basis = (0,1,2,2,2)

# Define double perovskite
class DoublePerovskiteFactory(SimpleTetragonalFactory):
    "A factory for creating double perovskite ((A1,A2)BO3) lattices"
    def __init__(self,oxygen=6.):
        self.oxygen=oxygen
        self.xtal_name = 'double perovskite'
        if self.oxygen==6:
            self.bravais_basis = [[0.5,0.5,0.25],[0.5,0.5,0.75],
                                 [0.,0.,0.],[0.,0.,0.5],
                                 [0.5,0.,0.],[0.,0.5,0.],[0.,0.,0.25],
                                 [0.5,0.,0.5],[0.,0.5,0.5],[0.,0.,0.75]]
            self.element_basis = (0,1,2,2,3,3,3,3,3,3)
        elif self.oxygen==5:
            self.bravais_basis = [[0.5,0.5,0.25],[0.5,0.5,0.75],
                                 [0.,0.,0.],[0.,0.,0.5],
                                 [0.5,0.,0.],[0.,0.5,0.],
                                 [0.5,0.,0.5],[0.,0.5,0.5],[0.,0.,0.75]]
            self.element_basis = (0,1,2,2,3,3,3,3,3)
        else:
            raise ValueError("oxygen keyword only accepts values 5 or 6")
