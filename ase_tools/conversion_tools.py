def cartesian_to_spherical(cartesian, center=None):
    """Convert cartesian to spherical coordinates."""
    from numpy import pi, arctan, arccos, array
    from numpy.linalg import norm
    if center is None:
        cartesian_new = cartesian
    else:
        cartesian_new = cartesian - center
    spherical = []
    for n, p in enumerate(cartesian_new):
        r = norm(p)
        x, y, z = p
        if x > 0 and y > 0:
            theta = arctan(y/x)
        elif x < 0 and y > 0:
            theta = pi + arctan(y/x)
        elif x > 0 and y < 0:
            theta = 2*pi + arctan(y/x)
        elif x < 0 and y < 0:
            theta = pi + arctan(y/x)
        elif x == 0 and y > 0:
            theta = pi/2
        elif x == 0 and y < 0:
            theta = 3*pi/2
        elif x > 0 and y == 0:
            theta = 0
        elif x < 0 and y == 0:
            theta = pi
        else:
            theta = 0
        if r != 0:
            phi = arccos(z/r)
        else:
            phi = 0
        spherical.append([r, theta, phi])
    spherical = array(spherical)
    return spherical

def cartesian_to_cylindrical(cartesian, center=None):
    """Convert cartesian to cylindrical coordinates. 'center' is the cartesian coordinate that will be used as origo in the cylindrical coordinate system."""
    from numpy import pi, arctan, arccos, array
    from numpy.linalg import norm
    if center is None:
        cartesian_new = cartesian
    else:
        cartesian_new = cartesian - center
    cylindrical = []
    for n, p in enumerate(cartesian_new):
        rho = norm(p[0:2])
        x, y, z = p
        if x > 0 and y > 0:
            theta = arctan(y/x)
        elif x < 0 and y > 0:
            theta = pi + arctan(y/x)
        elif x > 0 and y < 0:
            theta = 2*pi + arctan(y/x)
        elif x < 0 and y < 0:
            theta = pi + arctan(y/x)
        elif x == 0 and y > 0:
            theta = pi/2
        elif x == 0 and y < 0:
            theta = 3*pi/2
        elif x > 0 and y == 0:
            theta = 0
        elif x < 0 and y == 0:
            theta = pi
        else:
            theta = 0
        cylindrical.append([rho, theta, z])
    cylindrical = array(cylindrical)
    return cylindrical

def spherical_to_cartesian(spherical, center=None):
    """Convert spherical to cartesian coordinates."""
    from numpy import sin, cos
    cartesian_tmp = []
    for s in spherical:
        r, theta, phi = s
        x = r*cos(theta)*sin(phi)
        y = r*sin(theta)*sin(phi)
        z = r*cos(phi)
        cartesian_tmp.append([x, y, z])
    if cartesian is None:
        cartesian = np.array(cartesian_tmp)
    else:
        cartesian = np.array(cartesian_tmp) + center
    return cartesian

def cylindrical_to_cartesian(cylindrical, center=None):

    from numpy import sin, cos, array
    cartesian_tmp = []
    for c in cylindrical:
        rho, theta, z = c
        x = rho*cos(theta)
        y = rho*sin(theta)
        cartesian_tmp.append([x, y, z])
    if center is None:
        cartesian = array(cartesian_tmp)
    else:
        cartesian = array(cartesian_tmp) + center
    return cartesian
