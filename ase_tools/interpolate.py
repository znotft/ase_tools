from ase.io import read
from ase.all import *
from ase_tools.conversion_tools import cartesian_to_spherical, spherical_to_cartesian, cartesian_to_cylindrical, cylindrical_to_cartesian

def linear_interpolation(atoms1, atoms2, N):
    """
    Method that makes a linear interpolation between two atoms objects, returning N+2 images (including the two initial start and end point of the interpolation.
    """
    
    pos1 = atoms1.positions
    pos2 = atoms2.positions
    images = [atoms1]
    for n in range(N):
        nn = n + 1
        atoms_tmp = atoms1.copy()
        atoms_tmp.positions = ( (N+1-nn)*pos1 + nn*pos2 ) / (N+1)
        images += [atoms_tmp]
    images += [atoms2]
    return images

def cylindrical_interpolation(atoms1, atoms2, N, n_center, z_axis=None):

    # Rotate coordinates onto new z-axis, if this is defined
    if z_axis is not None:
        z_axis = z_axis / np.linalg.norm(z_axis)
        normal = np.cross([0,0,1], z_axis)
        angle = np.arccos(np.dot(z_axis, [0, 0, 1]))
    atoms1.rotate(normal, angle, [0, 0, 0])
    atoms2.rotate(normal, angle, [0, 0, 0])

    c1 = atoms1.positions.copy()[n_center].mean(axis=0)
    c2 = atoms2.positions[n_center].mean(axis=0)
    cylindrical1 = cartesian_to_cylindrical(atoms1.positions, c1)
    cylindrical2 = cartesian_to_cylindrical(atoms2.positions, c2)

    # Make sure that the shortest rotation is performed.
    for n in range(len(cylindrical1)):
        theta1 = cylindrical1[n,1]
        theta2 = cylindrical2[n,1]
        if theta1 - theta2 > np.pi:
            cylindrical2[n,1] += 2*np.pi
        elif theta2 - theta1 > np.pi:
            cylindrical1[n,1] += 2*np.pi

    # Interpolate
    atoms1.rotate(normal, -angle, [0,0,0])
    images = [atoms1]
    for n in range(N):
        nn = n + 1
        atoms_tmp = atoms1.copy()
        center = ( (N+1-nn)*c1 + nn*c2 ) / (N+1)
        cylindrical = ( (N+1-nn)*cylindrical1 + nn*cylindrical2 ) / (N+1)
        atoms_tmp.positions = cylindrical_to_cartesian(cylindrical, center)
        atoms_tmp.rotate(normal, -angle, [0,0,0])
        images += [atoms_tmp.copy()]
    atoms2.rotate(normal, -angle, [0,0,0])
    images += [atoms2.copy()]
    return images


def multi_interpolation(atoms1, atoms2, N, indices=None, multi_type='cylindrical', center=None, z_axis=None):
    """
    Will perform a normal linear interpolation for all atoms, except those with index in indices, 
    for which a cylindrical interpolation will be performed. The center of the cylindrical coordinate
    systems is defined by the COM of these atoms subsets if not the 'center' parameter is set. The z_coordinate
    of the cylindrical coordinate system with respect to the cartesian coordinate system is defined by z_axis.
    """
    images = linear_interpolation(atoms1, atoms2, N)
    if indices is None:
        # Perform normal linear interpolationa
        return images
    atoms1_cyl = atoms1[indices]
    atoms2_cyl = atoms2[indices]
    if center is None:
        center = range(len(atoms1_cyl))
    else:
        for n, c in enumerate(center):
            for m, ind in enumerate(indices):
                if c == ind:
                    center[n] = m
    images_cyl = cylindrical_interpolation(atoms1_cyl, atoms2_cyl, N, range(len(atoms1_cyl)), z_axis)
    for n in range(1, len(images)-1):
        images[n].positions[indices] = images_cyl[n].positions
    return images



# Unsupported methods

def spherical_interpolation(atoms1, atoms2, N, n_center):
    """
    Method that makes a interpolation between two atoms in a spherical coordinate system. The center of the coordinate system is defined by an atomic index. 
    """

    pos1 = atoms1.positions.copy()
    pos2 = atoms2.positions.copy()
    c1 = pos1[n_center].mean(axis=0)
    spherical1 = cartesian_to_spherical(pos1, c1)
    c2 = pos2[n_center].mean(axis=0)
    spherical2 = cartesian_to_spherical(pos2, c2)

    # Make sure that the shortest rotation is performed.
    for n in range(len(spherical1)):
        theta1 = spherical1[n,1]
        theta2 = spherical2[n,1]
        if theta1 - theta2 > np.pi:
            spherical2[n,1] += 2*np.pi
        elif theta2 - theta1 > np.pi:
            spherical1[n,1] += 2*np.pi

    # Interpolate
    images = [atoms1]
    for n in range(N):
        nn = n + 1
        atoms_tmp = atoms1.copy()
        center = ( (N+1-nn)*c1 + nn*c2 ) / (N+1)
        spherical = ( (N+1-nn)*spherical1 + nn*spherical2 ) / (N+1)
        atoms_tmp.positions = spherical_to_cartesian(spherical, center)
        images += [atoms_tmp]
    images += [atoms2.copy()]
    return images


