#!/usr/bin/env python
import numpy as np


def calculate_ellipse_radii(volume, thickness, ellipticity):
    pi = np.pi
    z_radius = thickness / 2.0
    x_radius = np.sqrt(volume * 3.0 / (4.0 * pi * ellipticity * z_radius))
    y_radius = ellipticity * x_radius

    return {'x_radius': x_radius, 'y_radius': y_radius, 'z_radius': z_radius}


def z_from_xy(x, y, x_radius, y_radius, z_radius):
    z = z_radius * np.sqrt(1.0 - (x / x_radius) ** 2 - (y / y_radius) ** 2)
    return z


def check_in_ellipsoid(x, y, z, x_radius, y_radius, z_radius):
    in_ellipsoid = False  # default to false
    coord_check = (x / x_radius) ** 2 + (y / y_radius) ** 2 + (z / z_radius) ** 2
    if coord_check < 1.0:
        in_ellipsoid = True

    return in_ellipsoid


def check_on_ellipsoid(x, y, z, x_radius, y_radius, z_radius):
    zero_tol = 1e-14
    on_ellipsoid = False  # default to false
    coord_check = (x / x_radius) ** 2 + (y / y_radius) ** 2 + (z / z_radius) ** 2
    if abs(coord_check - 1.0) < zero_tol:
        on_ellipsoid = True

    return on_ellipsoid


def angle_two_vectors(vector1, vector2):
    vector1_u = vector1 / np.linalg.norm(vector1)
    vector2_u = vector2 / np.linalg.norm(vector2)

    if (np.equal(vector1_u, vector2_u)).all():
        print('true')
        angle = 0.0
    else:
        dotprod = np.dot(vector1_u, vector2_u)
        if np.isclose(1.0, dotprod):
            angle = 0
        else:
            angle = np.arccos(dotprod)

    return angle


def element_connectivity_1D(node_loc, elems):
    # Initialise connectivity arrays
    num_elems = len(elems)
    elem_upstream = np.zeros((num_elems, 3), dtype=int)
    elem_downstream = np.zeros((num_elems, 3), dtype=int)
    num_nodes = len(node_loc)
    elems_at_node = np.zeros((num_nodes, 4), dtype=int)
    # determine elements that are associated with each node
    for ne in range(0, num_elems):
        for nn in range(1, 3):
            nnod = elems[ne][nn]
            elems_at_node[nnod][0] = elems_at_node[nnod][0] + 1
            elems_at_node[nnod][elems_at_node[nnod][0]] = ne
    # assign connectivity
    for ne in range(0, num_elems):
        nnod2 = elems[ne][2]  # second node in elem
        for noelem in range(1, elems_at_node[nnod2][0] + 1):
            ne2 = elems_at_node[nnod2][noelem]
            if ne2 != ne:
                elem_upstream[ne2][0] = elem_upstream[ne2][0] + 1
                elem_upstream[ne2][elem_upstream[ne2][0]] = ne
                elem_downstream[ne][0] = elem_downstream[ne][0] + 1
                elem_downstream[ne][elem_downstream[ne][0]] = ne2

    return {'elem_up': elem_upstream, 'elem_down': elem_downstream}
