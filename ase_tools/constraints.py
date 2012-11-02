import numpy as np
from numpy.linalg import norm, inv
#from numpy import sin, cos, arccos, arctan, pi
from math import sin, cos, acos, atan, pi
from ase.constraints import FixConstraint, slice2enlist


class FixVector(FixConstraint):
    """Constraint for keeping the direction of a vector between two atoms fixed.

    The atoms motion is non-correlated parallel to the vector, but correlated
    perpendicular to the vector. At the moment the two atoms must have the 
    same masses."""
    def __init__(self, indices=None):
        self.a = np.array(indices)

    def adjust_positions(self, old, new):
        v = old[self.a[0]] - old[self.a[1]]
        v /= norm(v)
        w1 = np.array([(v[0]*v[1]+v[1]*v[2])/v[0], -v[0], -v[1]])
        w1 /= norm(w1)
        w2 = np.cross(v, w1)
        cell = np.array([v, w1, w2])
        
        step = new - old
        mean_step = ( step[self.a[0]] + step[self.a[1]] ) / 2

        mean_step_scaled = np.dot(mean_step, inv(cell))
        step_scaled = np.dot(step, inv(cell))

        mean_step_s1 = mean_step_scaled.copy()
        mean_step_s2 = mean_step_scaled.copy()
        
        mean_step_s1[0] = step_scaled[self.a[0]][0]
        mean_step_s2[0] = step_scaled[self.a[1]][0]

        new[self.a[0]] = old[self.a[0]] + np.dot(mean_step_s1, cell)
        new[self.a[1]] = old[self.a[1]] + np.dot(mean_step_s2, cell)

    def adjust_forces(self, positions, forces):
        v = positions[self.a[0]] - positions[self.a[1]]
        v /= norm(v)
        w1 = np.array([(v[0]*v[1]+v[1]*v[2])/v[0], -v[0], -v[1]])
        w1 /= norm(w1)
        w2 = np.cross(v, w1)
        cell = np.array([v, w1, w2])

        mean_force = ( forces[self.a[0]]+forces[self.a[1]] )

        mean_force_scaled = np.dot(mean_force, inv(cell))
        forces_scaled = np.dot(forces, inv(cell))

        mean_force_s1 = mean_force_scaled.copy()
        mean_force_s2 = mean_force_scaled.copy()
        
        mean_force_s1[0] = forces_scaled[self.a[0]][0]
        mean_force_s2[0] = forces_scaled[self.a[1]][0]
        forces[self.a[0]] = np.dot(mean_force_s1, cell)
        forces[self.a[1]] = np.dot(mean_force_s2, cell)

    def index_shuffle(self, ind):
        # See docstring of superclass
        index = []
        for new, old in slice2enlist(ind):
            if old in self.a:
                index.append(new)
        if len(index) == 0:
            raise IndexError('All indices in FixVector not part of slice')
        self.a = np.asarray(index, int)

    def copy(self):
        return FixVector(self.atoms, indices=self.a.copy())
    
    def __repr__(self):
        if self.index.dtype == bool:
            return 'FixAtoms(mask=%s)' % ints2string(self.index.astype(int))
        return 'FixAtoms(indices=%s)' % ints2string(self.index)


class RigidBody(FixConstraint):
    """Treat a subset of atoms as a rigid body.

    It's possible to turn on/off translation and rotation of the rigid body
    with translate=True/False and rotate=True/False, respectively. Deafult
    is both translation and rotation turned on. (Only rotation may not be such
    a good idea in all situations.)

    masses: list of float
        The masses of the atoms object, to apply the constraint on,
        simply use atoms.get_masses().
    indices: list of int
        The indices of the atoms in the rigid body.
    translate: bool
        Allow translations of the rigid body.
    rotate: bool
        Allow rotations of the rigid body.
    """
    def __init__(self, masses, indices=None, translate=True, rotate=True):
        self.a = np.array(indices)
        self.m = masses
        self.translate = translate
        self.rotate = rotate

    def adjust_positions(self, old, new):
        step = new - old
        very_old = old.copy()

        # Calculate center of mass of rigid body
        com = np.dot(self.m[self.a], old[self.a])/self.m[self.a].sum()

        # Translation, move center of mass
        translation = np.zeros(3)
        if self.translate:
            for n in self.a:
                translation += self.m[n]*step[n]
            translation /= self.m[self.a].sum()
        old[self.a] += translation
        # New center of mass
        com_new = com+translation     

        if self.rotate:
            # Make rotations around axis [z, x, y](xyz).
            for xyz in range(3):
                com_xyz = np.array([com_new[(0+xyz)%3], com_new[(1+xyz)%3], com_new[(2+xyz)%3]])

                # Transform from cartesian to cylindrical coordinates
                new_cyl = np.empty([len(self.a), 3])
                old_cyl = np.empty([len(self.a), 3])
                for i, n in enumerate(self.a):
                    new_cyl[i][:] = cart2cyl([new[n][(0+xyz)%3], new[n][(1+xyz)%3], new[n][(2+xyz)%3]], com_xyz)
                    old_cyl[i][:] = cart2cyl([old[n][(0+xyz)%3], old[n][(1+xyz)%3], old[n][(2+xyz)%3]], com_xyz)
                step_phi = new_cyl.T[1]-old_cyl.T[1]
                # Calculate mean rotation around axis [z, x, y](xyz)
                mean_step_phi = 0
                for d_phi in step_phi:
                    if d_phi>pi:
                        d_phi -= 2*pi
                    elif d_phi<-pi:
                        d_phi += 2*pi
                    mean_step_phi += self.m[n]*d_phi
                mean_step_phi /= self.m[self.a].sum()
                step_cyl = np.array([0, mean_step_phi, 0])
                # Transform from cylindrical to cartesian coordinates
                new_cart = np.empty([len(self.a), 3])
                for i, n in enumerate(self.a):
                    ct = cyl2cart(old_cyl[i]+step_cyl, com_xyz)
                    new_cart[i] = [ct[(0-xyz)%3], ct[(1-xyz)%3], ct[(2-xyz)%3]]
                old[self.a] = new_cart

        # Update new positions
        new[self.a] = old[self.a]
        
        # Check that the atoms distance to COM has not changed.
        for n in self.a:
            a = norm(new[n]-(com_new))
            b = norm(very_old[n]-com)
            if abs(a - b) > 1E-10:
                raise RuntimeError('Distance between atom %i and COM has changed.' % n)

    def adjust_forces(self, positions, forces):
        f_tot = np.empty(len(forces))
        for i, f in enumerate(forces):
            f_tot[i] = norm(f)
        # Translational force
        mean_force = np.zeros(3)
        for n in self.a:
            mean_force += forces[n]

        forces[self.a] -= mean_force
        if not self.translate:
            mean_force = np.zeros(3)

        # Rotational force
        com = np.dot(self.m[self.a], positions[self.a])/self.m[self.a].sum()
        forces_rotate = np.zeros([len(self.a), 3])
        if self.rotate:
            for xyz in range(3):
                tot_length = 0
                com_tmp = np.array([com[(0+xyz)%3], com[(1+xyz)%3], com[(2+xyz)%3]])
                tau_phi = 0
                for i, n in enumerate(self.a):
                    pos_tmp = np.array([positions[n, (0+xyz)%3], positions[n, (1+xyz)%3], positions[n, (2+xyz)%3]])
                    unit = cyl_unit(pos_tmp, com_tmp)
                    rho = cart2cyl(pos_tmp, com_tmp)[0]
                    tau_phi += np.dot(forces[n], unit[1])*rho
                for i, n in enumerate(self.a):
                    pos_tmp = np.array([positions[n, (0+xyz)%3], positions[n, (1+xyz)%3], positions[n, (2+xyz)%3]])
                    unit = cyl_unit(pos_tmp, com_tmp)
                    rho = cart2cyl(pos_tmp, com_tmp)[0]
                    if rho>0:
                        forces_rotate[i] += tau_phi*unit[1]/rho
        forces[self.a] = forces_rotate + mean_force

    def index_shuffle(self, ind):
        # See docstring of superclass
        #index = []
        #for new, old in slice2enlist(ind):
        #    if old in self.a:
        #        index.append(new)
        #if len(index) == 0:
        #    raise IndexError('All indices in RigidBody not part of slice')
        #self.a = np.asarray(index, int)
        1

    def copy(self):
        return RigidBody(masses=self.m.copy(), indices=self.a.copy())
    
    def __repr__(self):
        return 'RigidBody(indices=%s)' % ints2string(self.a.astype(int))

def cart2cyl(coords, com):
    coords = np.array(coords)-np.array(com)
    r = norm(coords[0:2])
    if coords[0]>0 and coords[1]>=0:
        phi = atan(coords[1]/coords[0])
    elif coords[0]<0 and coords[1]<0:
        phi = atan(coords[1]/coords[0])+pi
    elif coords[0]<0 and coords[1]>=0:
        phi = pi + atan(coords[1]/coords[0])
    elif coords[0]>0 and coords[1]<0:
        phi = 2*pi + atan(coords[1]/coords[0])
    elif coords[0]==0 and coords[1]>0:
        phi = pi/2
    elif coords[0]==0 and coords[1]<0:
        phi = -pi/2
    elif coords[0]==0 and coords[1]==0:
        phi = 0
    z = coords[2]
    return np.array([r, phi, z])

def cyl2cart(coords, com):
    [r, phi, z] = coords
    x = r*cos(phi)
    y = r*sin(phi)
    z = z
    return np.array([x, y, z])+np.array(com)

def cyl_unit(coords_cart, com):
    """Returns cylindrical unit vectors in cartesian coordinates."""
    [r, phi, z] = cart2cyl(coords_cart, com)
    e_x = np.array([1, 0, 0])
    e_y = np.array([0, 1, 0])
    e_z = np.array([0, 0, 1])
    e_r = cos(phi)*e_x + sin(phi)*e_y
    e_phi = -sin(phi)*e_x + cos(phi)*e_y
    return np.array([e_z, e_phi, e_z])

