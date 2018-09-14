#!/usr/bin/env python
import numpy as np
from scipy.spatial import Delaunay

from . import pg_utilities

"""
.. module:: generate_shapes
  :synopsis: Contains code to generate placental shapes for generic placental models.

:synopsis:Contains code to generate placental shapes for generic placental models \n
 (i.e. from literature measures without specific data from an individual
 
"""


def equispaced_data_in_ellipsoid(n, volume, thickness, ellipticity):
    """ Generates equally spaced data points in an ellipsoid.

    Inputs:
       - n: number of data points which we aim to generate
       - volume: volume of ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions

    Returns:
       - Edata: A nx3 array of datapoints, with each point being defined by its x-,y-, and z- coordinates

    A way you might want to use me is:

    >>> n = 100
    >>> volume = 10
    >>> thickness = 3
    >>> ellipticity = 1.1
    >>> equispaced_data_in_ellipsoid(n, volume, thickness, ellipticity)

   This will return 100 data points in an ellipse with z-axis thickness 3, volume 10, and with the y-axis dimension 1.1 times the x-axis dimension.

    """
    # Generates equally spaced data points in an ellipsoid with the following inputs
    # n=number of data points which we aim to generate
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
    data_spacing = (volume / n) ** (1.0 / 3.0)
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    # Aiming to generate seed points that fill a cuboid encompasing the placental volume then remove seed points that
    # are external to the ellipsoid

    num_data = 0  # zero the total number of data points

    # Calculate the number of points that should lie in each dimension in a cube
    nd_x = np.floor(2.0 * (x_radius + data_spacing) / data_spacing)
    nd_y = np.floor(2.0 * (y_radius + data_spacing) / data_spacing)
    nd_z = np.floor(2.0 * (z_radius + data_spacing) / data_spacing)
    nd_x = int(nd_x)
    nd_y = int(nd_y)
    nd_z = int(nd_z)
    # Set up edge node coordinates
    x_coord = np.linspace(-x_radius - data_spacing / 2.0, x_radius + data_spacing / 2.0, nd_x)
    y_coord = np.linspace(-y_radius - data_spacing / 2.0, y_radius + data_spacing / 2.0, nd_y)
    z_coord = np.linspace(-z_radius - data_spacing / 2.0, z_radius + data_spacing / 2.0, nd_z)

    # Use these vectors to form a unifromly spaced grid
    data_coords = np.vstack(np.meshgrid(x_coord, y_coord, z_coord)).reshape(3, -1).T

    # Store nodes that lie within ellipsoid
    Edata = np.zeros((nd_x * nd_y * nd_z, 3))
    for i in range(len(data_coords)):  # Loop through grid
        coord_check = pg_utilities.check_in_ellipsoid(data_coords[i][0], data_coords[i][1], data_coords[i][2], x_radius,
                                                      y_radius, z_radius)

        if coord_check is True:  # Has to be strictly in the ellipsoid
            Edata[num_data, :] = data_coords[i, :]  # add to data array
            num_data = num_data + 1
    Edata.resize(num_data, 3)  # resize data array to correct size

    print('Data points within ellipsoid allocated. Total = ' + str(len(Edata)))

    return Edata


def uniform_data_on_ellipsoid(n, volume, thickness, ellipticity, random_seed):
    """ Generates equally spaced data points on the positive z-surface of an ellipsoid

    Inputs:
       - n: number of data points which we aim to generate
       - volume: volume of ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions

    Returns:
       - chorion_data: A nx3 array of datapoints, with each point being defined by its x-,y-, and z- coordinates

    A way you might want to use me is:

    >>> n = 100
    >>> volume = 10
    >>> thickness = 3
    >>> ellipticity = 1.1
    >>> equispaced_data_on_ellipsoid(n, volume, thickness, ellipticity)

   This will return 100 data points on the positive z-surface ellipse with z-axis thickness 3, volume 10, and with the y-axis dimension 1.1 times the x-axis dimension.

    """
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    area_estimate = np.pi * x_radius * y_radius
    data_spacing = 0.85 * np.sqrt(area_estimate / n)

    chorion_data = np.zeros((n, 3))
    np.random.seed(random_seed)
    generated_seed = 0
    acceptable_attempts = n * 1000  # try not to have too many failures
    attempts = 0

    while generated_seed < n and attempts < acceptable_attempts:
        # generate random x-y coordinates between negative and positive radii
        new_x = np.random.uniform(-x_radius, x_radius)
        new_y = np.random.uniform(-y_radius, y_radius)
        # check if new coordinate is on the ellipse
        if ((new_x / x_radius) ** 2 + (new_y / y_radius) ** 2) < 1:  # on the surface
            if generated_seed == 0:
                generated_seed = generated_seed + 1
                new_z = pg_utilities.z_from_xy(new_x, new_y, x_radius, y_radius, z_radius)
                chorion_data[generated_seed - 1][:] = [new_x, new_y, new_z]
            else:
                reject = False
                for j in range(0, generated_seed + 1):
                    distance = (chorion_data[j - 1][0] - new_x) ** 2 + (chorion_data[j - 1][1] - new_y) ** 2
                    distance = np.sqrt(distance)
                    if distance <= data_spacing:
                        reject = True
                        break
                if reject is False:
                    generated_seed = generated_seed + 1
                    new_z = pg_utilities.z_from_xy(new_x, new_y, x_radius, y_radius, z_radius)
                    chorion_data[generated_seed - 1][:] = [new_x, new_y, new_z]

        attempts = attempts + 1
    chorion_data.resize(generated_seed, 3)  # resize data array to correct size
    print('Data points on ellipsoid allocated. Total = ' + str(len(chorion_data)))

    return chorion_data


def gen_rect_cover_ellipsoid(volume, thickness, ellipticity, x_spacing, y_spacing, z_spacing):
    # Generates equally spaced data nodes and elements and constructs a rectangular 'mesh' that covers the space that is
    # made up of an ellipsoidal placenta
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
    # X,Y,Z spacing is the number of elements required in each of the x, y z directions

    # Calculate the dimensions of the ellipsoid
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    # z height of ellipsoid is 2* zradius
    # We want number of nodes to cover height and have prescribed spaing
    nnod_x = int(np.ceil(x_radius * 2.0 / x_spacing)) + 1
    x_width = x_spacing * (nnod_x - 1)
    nnod_y = int(np.ceil(y_radius * 2.0 / y_spacing)) + 1
    y_width = y_spacing * (nnod_y - 1)
    nnod_z = int(np.ceil(z_radius * 2.0 / z_spacing)) + 1
    z_width = z_spacing * (nnod_z - 1)

    node_loc = gen_rectangular_node(x_width, y_width, z_width, nnod_x, nnod_y, nnod_z)
    # Generating the element connectivity of each cube element, 8 nodes for each 3D cube element
    elems = cube_mesh_connectivity(nnod_x, nnod_y, nnod_z)

    return {'nodes': node_loc, 'elems': elems, 'total_nodes': nnod_x * nnod_y * nnod_z,
            'total_elems': (nnod_x - 1) * (nnod_y - 1) * (nnod_z - 1)}


def gen_ellip_mesh_tet(volume, thickness, ellipticity, n):
    """ Generates ellipsoid tetrahedral mesh for 3D problems

    Inputs:
       - volume: volume of placental ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions
       - n: number of datapoints generated to create the mesh

    Returns:
       - nodes: nodes location of mesh
       - elems: element connectivity of mesh (tetrahedral element)
       - node_array: array of nodes
       - element_array: array of elements
    """

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    nodeSpacing = (n / (2 * x_radius * 2 * y_radius * 2 * z_radius)) ** (1. / 3)

    nnod_x = 2 * x_radius * nodeSpacing
    nnod_y = 2 * y_radius * nodeSpacing
    nnod_z = 2 * z_radius * nodeSpacing
    nodes = gen_rectangular_node(x_radius * 2, y_radius * 2, z_radius * 2, nnod_x, nnod_y, nnod_z)

    # nodes inside the ellipsoid
    ellipsoid_node = np.zeros((len(nodes), 3))
    count = 0
    for nnode in range(0, len(nodes)):
        coord_point = nodes[nnode][0:3]
        inside = pg_utilities.check_in_on_ellipsoid(coord_point[0], coord_point[1], coord_point[2], x_radius, y_radius,
                                                    z_radius)
        if inside:
            ellipsoid_node[count, :] = coord_point[:]
            count = count + 1
    ellipsoid_node.resize(count, 3)
    xyList = ellipsoid_node[:, [0, 1]]
    xyListUnique = np.vstack({tuple(row) for row in xyList})
    # looking for z_coordinate of surface nodes
    for xyColumn in xyListUnique:

        xyNodes = np.where(np.all(xyList == xyColumn, axis=1))[0]
        if len(xyNodes) > 1:
            x_coord = ellipsoid_node[xyNodes[0], 0]
            y_coord = ellipsoid_node[xyNodes[0], 1]
            ellipsoid_node[xyNodes[len(xyNodes) - 1], 2] = pg_utilities.z_from_xy(x_coord, y_coord, x_radius, y_radius,
                                                                                  z_radius)
            ellipsoid_node[xyNodes[0], 2] = -1 * (
                pg_utilities.z_from_xy(x_coord, y_coord, x_radius, y_radius, z_radius))

    # generate tetrahedral mesh
    pyMesh = Delaunay(ellipsoid_node)

    # Build arrays to pass into openCMISS conversion:
    node_loc = pyMesh.points
    temp_elems = pyMesh.simplices
    # CHECK ELEMENTS FOR 0 VOLUME:
    min_vol = 0.00001
    index = 0
    indexArr = []

    for element in temp_elems:
        x_coor = []
        y_coor = []
        z_coor = []
        for node in element:
            x_coor.append(node_loc[node][0])
            y_coor.append(node_loc[node][1])
            z_coor.append(node_loc[node][2])

        vmat = np.vstack((x_coor, y_coor, z_coor, [1.0, 1.0, 1.0, 1.0]))  # matrix of coor of element
        elem_volume = (1 / 6.0) * abs(np.linalg.det(vmat))  # volume of each tetrahedral element

        # if volume is not zero
        if elem_volume > min_vol:
            indexArr.append(index)
        index = index + 1

    # update arrays without 0 volume elements, to pass into openCMISS
    elems = temp_elems[indexArr, :]
    for i in range(len(elems)):
        elems[i] = [x + 1 for x in elems[i]]
    element_array = range(1, len(elems) + 1)
    node_array = range(1, len(node_loc) + 1)

    return {'nodes': node_loc, 'elems': elems, 'element_array': element_array, 'node_array': node_array,
            'nodeSpacing': nodeSpacing}


def gen_rectangular_node(x_width, y_width, z_width, nnod_x, nnod_y, nnod_z):
    # Create linspaces for x y and z coordinates
    x = np.linspace(-x_width / 2.0, x_width / 2.0, nnod_x)  # linspace for x axis
    y = np.linspace(-y_width / 2.0, y_width / 2.0, nnod_y)  # linspace for y axis
    z = np.linspace(-z_width / 2.0, z_width / 2.0, nnod_z)  # linspace for z axis
    node_loc_temp = np.vstack(np.meshgrid(y, z, x)).reshape(3, -1).T  # generate nodes for rectangular mesh
    node_loc = np.zeros((len(node_loc_temp), 3))
    for i in range(0, len(node_loc)):
        node_loc[i][0] = node_loc_temp[i][2]
        node_loc[i][1] = node_loc_temp[i][0]
        node_loc[i][2] = node_loc_temp[i][1]

    return node_loc


def gen_rectangular_mesh2(nel_x, nel_y, nel_z, xdim, ydim, zdim, element_type):
    # generates a rectangular mesh of defined dimenions using either linear or quadratic elements
    if element_type == 1:  # linear element
        nnod_x = int(nel_x + 1)
        nnod_y = int(nel_y + 1)
        nnod_z = int(nel_z + 1)
    elif element_type == 2:  # quadratic element
        nnod_x = int((nel_x * 2) + 1)
        nnod_y = int((nel_y * 2) + 1)
        nnod_z = int((nel_z * 2) + 1)

    node = gen_rectangular_node(xdim, ydim, zdim, nnod_x, nnod_y, nnod_z)  # getting nodes

    if element_type == 1:  # linear element
        elems = cube_mesh_connectivity(nnod_x, nnod_y, nnod_z)  # getting elem connectivity
    elif element_type == 2:  # quadratic element
        elems = cube_mesh_connectivity_quadratic(nel_x, nel_y, nel_z, nnod_x, nnod_y,
                                                 nnod_z)  # getting element connectivity

    element_array = range(1, len(elems) + 1)
    node_array = range(1, len(node) + 1)
    if element_type == 2:
        surfacenodes = identify_surface_node_quad(nel_x, nel_y, nel_z)
    else:
        print("This element type has no implemented surface node definition")
        surfacenodes = 0

    return {'nodes': node, 'elems': elems, 'element_array': element_array,
            'node_array': node_array, 'surface_nodes': surfacenodes}


def gen_3d_ellipsoid(nel_x, nel_y, nel_z, volume, thickness, ellipticity, element_type):
    """ Generates ellipsoid placental mesh to solve 3D problems (note this is not a quality structured mesh)

    Inputs:
       - nel: number of element in x,y,z axis , the more nel, the rounder the mesh
       - volume: volume of placental ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions
       
    Returns:
       - placental_node_coor: nodes location of mesh
       - placental_el_con: element connectivity of mesh (tetrahedral element)
       - node_array: array of nodes
       - element_array: array of elements
    """
    # creating cube between -1 and 1 with n number of element
    # cubelength=2
    if element_type == 1:  # linear element
        nnod_x = int(nel_x + 1)
        nnod_y = int(nel_y + 1)
        nnod_z = int(nel_z + 1)
    elif element_type == 2:  # quadratic element
        nnod_x = int((nel_x * 2) + 1)
        nnod_y = int((nel_y * 2) + 1)
        nnod_z = int((nel_z * 2) + 1)

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    cube_node = gen_rectangular_node(2 * x_radius, 2 * y_radius, 2 * z_radius, nnod_x, nnod_y, nnod_z)
    if element_type == 1:  # linear element
        cube_elems = cube_mesh_connectivity(nnod_x, nnod_y, nnod_z)  # getting elem connectivity
    elif element_type == 2:  # quadratic element
        cube_elems = cube_mesh_connectivity_quadratic(nel_x, nel_y, nel_z, nnod_x, nnod_y,
                                                      nnod_z)  # getting element connectivity

    ellipsoid_coor = np.zeros((len(cube_node), 3))

    for ii in range(0, len(cube_node)):
        ellipsoid_coor[ii, 0] = cube_node[ii, 0] * np.sqrt(1.0 - cube_node[ii, 1] ** 2 / (2.0 * y_radius ** 2) -
                                                           cube_node[ii, 2] ** 2 / (2.0 * z_radius ** 2) + cube_node[
                                                               ii, 1] ** 2 *
                                                           cube_node[ii, 2] ** 2 / (
                                                                       3.0 * y_radius ** 2 * z_radius ** 2))  # for  x_coor
        ellipsoid_coor[ii, 1] = cube_node[ii, 1] * np.sqrt(1.0 - cube_node[ii, 0] ** 2 / (2.0 * x_radius ** 2) -
                                                           cube_node[ii, 2] ** 2 / (2.0 * z_radius ** 2) + cube_node[
                                                               ii, 0] ** 2 * cube_node[ii, 2] ** 2
                                                           / (3.0 * x_radius ** 2 * z_radius ** 2))  # for  y_coor
        ellipsoid_coor[ii, 2] = cube_node[ii, 2] * np.sqrt(1.0 - cube_node[ii, 1] ** 2 / (2.0 * y_radius ** 2) -
                                                           cube_node[ii, 0] ** 2 / (2.0 * x_radius ** 2) + cube_node[
                                                               ii, 1] ** 2 * cube_node[ii, 0] ** 2
                                                           / (3.0 * y_radius ** 2 * x_radius ** 2))  # for  z_coor

    element_array = range(1, len(cube_elems) + 1)
    node_array = range(1, len(ellipsoid_coor) + 1)
    if element_type == 2:
        surfacenodes = identify_surface_node_quad(nel_x, nel_y, nel_z)
    else:
        print("This element type has no implemented surface node definition")
        surfacenodes = 0

    return {'placental_node_coor': ellipsoid_coor, 'placental_el_con': cube_elems, 'element_array': element_array,
            'node_array': node_array, 'surface_nodes': surfacenodes}


def cube_mesh_connectivity(nnod_x, nnod_y, nnod_z):
    """Generates element connectivity in cube mesh
      
       Inputs:
         - nnod_x:number of node in x axis
         - nnod_y:number of node in y axis
         - nnod_z:number of node in z axis

       Outputs:
         - elems: array of element connectivity

    """
    num_elems = (nnod_x - 1) * (nnod_y - 1) * (nnod_z - 1)
    elems = np.zeros((num_elems, 9),
                     dtype=int)  # this stores first element number and then the nodes of each mesh element
    element_number = 0
    ne = 0
    # loop through elements
    for k in range(1, nnod_z):
        for j in range(1, nnod_y):
            for i in range(1, nnod_x):
                elems[ne][0] = ne  # store element number
                elems[ne][1] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1)  # lowest coordinates
                elems[ne][2] = elems[ne][1] + 1  # add one in x
                elems[ne][3] = elems[ne][1] + nnod_x  # go through x and find first in y
                elems[ne][4] = elems[ne][3] + 1  # add one in y
                elems[ne][5] = elems[ne][1] + nnod_x * nnod_y  # same as 1 -4 but at higher z -coord
                elems[ne][6] = elems[ne][2] + nnod_x * nnod_y
                elems[ne][7] = elems[ne][3] + nnod_x * nnod_y
                elems[ne][8] = elems[ne][4] + nnod_x * nnod_y
                ne = ne + 1

    return elems


def cube_mesh_connectivity_quadratic(nel_x, nel_y, nel_z, nnod_x, nnod_y, nnod_z):
    """Generates element connectivity in quadratic cube mesh
      
       Inputs:
         - nnod_x:number of node in x axis
         - nnod_y:number of node in y axis
         - nnod_z:number of node in z axis

       Outputs:
         - elems: array of element connectivity in quadratic

    """

    num_elems = nel_x * nel_y * nel_z
    elems = np.zeros((num_elems, 28), dtype=int)

    element_number = 0
    ne = 0
    # Got the element
    for k in range(1, nnod_z, 2):
        for j in range(1, nnod_y, 2):
            for i in range(1, nnod_x, 2):
                # 1st layer
                elems[ne][0] = ne
                elems[ne][1] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1)  # 1st node
                elems[ne][2] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1) + 1  # right subsequent node
                elems[ne][3] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1) + 2  # right subsequent node
                elems[ne][4] = elems[ne][1] + nnod_x  # 1st node in another y layer
                elems[ne][5] = elems[ne][1] + nnod_x + 1  # right subsequent node
                elems[ne][6] = elems[ne][1] + nnod_x + 2  # right subsequent node
                elems[ne][7] = elems[ne][1] + 2 * (nnod_x)  # 1st node in another y layer
                elems[ne][8] = elems[ne][1] + 2 * (nnod_x) + 1  # right subsequent node
                elems[ne][9] = elems[ne][1] + 2 * (nnod_x) + 2  # right subsequent node

                # 2nd layer
                elems[ne][10] = elems[ne][1] + nnod_x * nnod_y  # same in one z layer
                elems[ne][11] = elems[ne][2] + nnod_x * nnod_y
                elems[ne][12] = elems[ne][3] + nnod_x * nnod_y
                elems[ne][13] = elems[ne][4] + nnod_x * nnod_y
                elems[ne][14] = elems[ne][5] + nnod_x * nnod_y
                elems[ne][15] = elems[ne][6] + nnod_x * nnod_y
                elems[ne][16] = elems[ne][7] + nnod_x * nnod_y
                elems[ne][17] = elems[ne][8] + nnod_x * nnod_y
                elems[ne][18] = elems[ne][9] + nnod_x * nnod_y

                # thrid layer

                elems[ne][19] = elems[ne][1] + nnod_x * nnod_y * 2  # same in another z layer
                elems[ne][20] = elems[ne][2] + nnod_x * nnod_y * 2
                elems[ne][21] = elems[ne][3] + nnod_x * nnod_y * 2
                elems[ne][22] = elems[ne][4] + nnod_x * nnod_y * 2
                elems[ne][23] = elems[ne][5] + nnod_x * nnod_y * 2
                elems[ne][24] = elems[ne][6] + nnod_x * nnod_y * 2
                elems[ne][25] = elems[ne][7] + nnod_x * nnod_y * 2
                elems[ne][26] = elems[ne][8] + nnod_x * nnod_y * 2
                elems[ne][27] = elems[ne][9] + nnod_x * nnod_y * 2

                ne = ne + 1

    return elems


def identify_surface_node_quad(nel_x, nel_y, nel_z):
    """Generates collection of nodes that are on the surface of in quadratic placental mesh
      
       Inputs:
         - nel_x:number of elem in x axis
         - nel_y:number of elem in y axis
         - nel_z:number of elem in z axis

       Outputs:
         - surfacenode: collection of nodes on the surface of placental mesh

    """
    nnod_x = int((nel_x * 2) + 1)  # number of nodes in x axis
    nnod_y = int((nel_y * 2) + 1)  # number of nodes in y axis
    nnod_z = int((nel_z * 2) + 1)  # number of nodes in z axis
    # For left and right surface
    sIEN = np.zeros((9, nel_y * nel_z), dtype=int)  # to store surface indiviaul element nodes (sIEN)
    e = 0
    for k in range(1, nnod_x * nnod_y * (nnod_z - 1), (nnod_x * nnod_y) * 2):  # go up
        for j in range(1, nnod_x * (nnod_y - 1), 2 * nnod_x):  # go left

            sIEN[0, e] = j + (k - 1)  # 1st node
            sIEN[1, e] = sIEN[0, e] + (nnod_x) * (nnod_y)  # 2nd node
            sIEN[2, e] = sIEN[1, e] + (nnod_x) * (nnod_y)  # 3rd node
            sIEN[3, e] = sIEN[0, e] + nnod_x  # 4th node
            sIEN[4, e] = sIEN[1, e] + nnod_x  # 5th node
            sIEN[5, e] = sIEN[2, e] + nnod_x  # 6th node
            sIEN[6, e] = sIEN[3, e] + nnod_x  # 7th node
            sIEN[7, e] = sIEN[4, e] + nnod_x  # 8th node
            sIEN[8, e] = sIEN[5, e] + nnod_x  # 9th node
            e = e + 1

        left = np.unique(sIEN)  # collection of nodes of left surface
        right = np.unique(sIEN.T + (nnod_x - 1))  # collection of nodes on right surface

    # For front and back surface
    sIEN = np.zeros((9, nel_x * nel_z), dtype=int)
    e = 0
    for k in range(1, nnod_x * nnod_y * (nnod_z - 2), (nnod_x * nnod_y) * 2):  # go up
        for i in range(1, nnod_x - 1, 2):  # go right
            sIEN[0, e] = i + (k - 1)
            sIEN[1, e] = sIEN[0, e] + 1
            sIEN[2, e] = sIEN[0, e] + 2
            sIEN[3, e] = sIEN[0, e] + (nnod_x * nnod_y)
            sIEN[4, e] = sIEN[3, e] + 1
            sIEN[5, e] = sIEN[3, e] + 2
            sIEN[6, e] = sIEN[3, e] + (nnod_x * nnod_y)
            sIEN[7, e] = sIEN[6, e] + 1
            sIEN[8, e] = sIEN[6, e] + 2
            e = e + 1

        front = np.unique(sIEN)  # collection of nodes on front surface
        back = np.unique(sIEN.T + (nnod_x * (nnod_y - 1)))  # collection of nodes on back surface

    # For top and bottom surface
    sIEN = np.zeros((9, nel_x * nel_y), dtype=int)
    e = 0
    for j in range(1, nnod_x * (nnod_y - 1), nnod_x * 2):  # go up
        for i in range(1, nnod_x - 1, 2):  # go back
            sIEN[0, e] = i + (j - 1)
            sIEN[1, e] = sIEN[0, e] + 1
            sIEN[2, e] = sIEN[0, e] + 2
            sIEN[3, e] = sIEN[0, e] + nnod_x
            sIEN[4, e] = sIEN[3, e] + 1
            sIEN[5, e] = sIEN[3, e] + 2
            sIEN[6, e] = sIEN[3, e] + nnod_x
            sIEN[7, e] = sIEN[6, e] + 1
            sIEN[8, e] = sIEN[6, e] + 2
            e = e + 1

        bottom = np.unique(sIEN)  # collection of nodes on bottom surface
        top = np.unique(sIEN.T + (nnod_x * nnod_y) * (nnod_z - 1))  # collection of nodes on top surface

    surfacenode = np.hstack((front, back, left, right, bottom, top))
    surfacenode = np.unique(surfacenode)  # collection of surface nodes from all surface
    return surfacenode


def identify_node_from_coord(nodes, filename):
    # reading in the node location
    xyz = open(filename, 'r')
    xyz_coor = xyz.readlines()  # readlines
    startLines = range(0, len(xyz_coor))

    for i in range(len(xyz_coor)):
        xyz_coor[i] = xyz_coor[i].split()
    xyzList = []
    for i in startLines:
        targetpoint = []
        targetpoint.append(float(xyz_coor[i][0]))  # x coor
        targetpoint.append((float(xyz_coor[i][1])))  # y coor
        targetpoint.append((float(xyz_coor[i][1])))  # y coor
        xyzList.append(targetpoint)
    xyz.close()

    node_list = np.zeros(len(xyzList))

    mindist = 100000
    for i in range(0, len(xyzList)):
        for j in range(0, len(nodes)):
            print(xyzList[i][0], nodes[j][0])

    return i


def identify_vessel_node(ellipsoid_coor, surfacenode, stem_file, volume,thickness, ellipticity):
    """Generates array of spiral artery nodes and decidual vein nodes. Spiral artery nodes are mapped with stem villi.
      
       Inputs:
         - ellipsoid_coor:coordinate of nodes of placental mesh
         - surfacenode:array of surface nodes
         - stem_file:txt file that described stem villi locations

       Outputs:
         - spiral_array: array of spiral artery nodes 
         - decidual_array: array of decidual artery nodes
         - vesselnode: array of both spiral and decidual nodes
         - surfnode_ex_vessel: array of surface node excluding vessel nodes
    """
    sa_radius = 3.7/2.0
    dv_radius = sa_radius

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    xyList = np.zeros((len(surfacenode), 4))
    count = 0
    for i in range(0, len(surfacenode)):  # taking only x and y coordinates
        if ellipsoid_coor[surfacenode[i] - 1, 3] > 0:  # take if upper surface nodes as this is where vessele reside
            # location from upper surface nodes only
            xyList[count, 0] = ellipsoid_coor[surfacenode[i] - 1, 0] #node number
            xyList[count, 1] = ellipsoid_coor[surfacenode[i] - 1, 1] #x-coordinate
            xyList[count, 2] = ellipsoid_coor[surfacenode[i] - 1, 2]  #y-coordinate
            xyList[count, 3] = ellipsoid_coor[surfacenode[i] - 1, 3]   #z-coordinate

            count = count + 1

    xyList = xyList[0:count, :]


    surfnode_ex_vessel = np.copy(surfacenode)
    vesselnode_temp = np.vstack({tuple(row) for row in xyList})  #nodes that might be vessels

    # reading in the stem vessel to map the spiral artery location
    stem_xy = open(stem_file, 'r')
    stem_coor = stem_xy.readlines()  # readlines
    startLines = range(0, len(stem_coor))

    for i in range(len(stem_coor)):
        stem_coor[i] = stem_coor[i].split()
    stem_xyList = []
    for i in startLines:
        node = []
        node.append(float(stem_coor[i][0]))  # x coor of stem villi
        node.append((float(stem_coor[i][1])))  # y coor of stem villi
        stem_xyList.append(node)
    stem_xy.close()

    print('Total stem read = '+ str(len(stem_xyList)))

    vessel_mapped_stem = stem_xyList  # this is the x,y location where we want to put spiral artery
    spiral_array = np.zeros((len(xyList)), dtype=int)  # store the node nuber of spiral artery
    decidual_array = np.zeros((len(xyList)), dtype=int)  # store the node number of decidual vein

    check = ellipsoid_coor[:, 0:2]
    np.random.seed(0)
    sa_nodes = 0
    dv_nodes = 0
    for i in range(0, len(vessel_mapped_stem)):  # for each blood vessel,Cycle through to find closest nodes

        closest_node = 0
        for nodeX in vesselnode_temp:
            distance=np.sqrt((vessel_mapped_stem[i][0] - nodeX[1]) ** 2 + (
                        vessel_mapped_stem[i][1] - nodeX[2]) ** 2 )  # distance from the nodes
            if(distance < sa_radius):
                #print('SA Node', int(nodeX[0]),nodeX[1],nodeX[2],vessel_mapped_stem[i][0],vessel_mapped_stem[i][1])
                arterynode = nodeX[0]
                A = np.where(vesselnode_temp == arterynode)
                vesselnode_temp = np.delete(vesselnode_temp, A[0], axis=0)
                A2 = np.where(surfnode_ex_vessel == int(arterynode))
                surfnode_ex_vessel = np.delete(surfnode_ex_vessel, A2)
                spiral_array[sa_nodes] = arterynode
                sa_nodes = sa_nodes +1

        #print(closest_node[0])
        #arterynode = closest_node[0]
        #A = np.where(vesselnode_temp == arterynode)
        #vesselnode_temp =  np.delete(vesselnode_temp, A[0], axis=0)
        #A2 = np.where(surfnode_ex_vessel == int(arterynode))
        #surfnode_ex_vessel = np.delete(surfnode_ex_vessel, A2)
        #spiral_array[i] = arterynode
        #sa_nodes = sa_nodes +1

    #Doing decidual veins after arteries to make sure we dont take up any spots that arteries would have otherwise beein
    for i in range(0, len(vessel_mapped_stem)): #need same number of arteries as veins
        V = np.random.choice(len(vesselnode_temp))  # choosing random , won't repeat arteries as they are already
        vessel_location = vesselnode_temp[V]
        for nodeX in vesselnode_temp:
            distance=np.sqrt((vessel_location[1] - nodeX[1]) ** 2 + (
                        vessel_location[2] - nodeX[2]) ** 2 )  # distance from the nodes
            dv_from_centre = np.sqrt(nodeX[1] ** 2 +  nodeX[2] ** 2 )
            if(distance < dv_radius and dv_from_centre < 0.9*x_radius):
                #print('DV Node', int(nodeX[0]))
                veinnode = nodeX[0]
                V = np.where(vesselnode_temp == veinnode)
                vesselnode_temp = np.delete(vesselnode_temp, V[0], axis=0)
                V2 = np.where(surfnode_ex_vessel == int(veinnode))
                surfnode_ex_vessel = np.delete(surfnode_ex_vessel, V2)
                decidual_array[dv_nodes] = veinnode
                dv_nodes = dv_nodes +1


        #veinnode = vesselnode_temp[V][0]
        #vesselnode_temp = np.delete(vesselnode_temp, V, axis=0)
        #V2 = np.where(surfnode_ex_vessel == int(veinnode))
        #surfnode_ex_vessel = np.delete(surfnode_ex_vessel, V2)
        #decidual_array[i] = veinnode
        #dv_nodes = dv_nodes+1

    spiral_array = np.resize(spiral_array,sa_nodes)
    print('SAs found = ' + str(sa_nodes))
    decidual_array = np.resize(decidual_array, dv_nodes)
    #print('dec',decidual_array)

    return {'spiral_array': spiral_array, 'decidual_array': decidual_array, 'surfnode_ex_vessel': surfnode_ex_vessel}


def identify_vessel_node_test_mesh(ellipsoid_coor, surfacenode,volume, thickness, ellipticity):
    """Generates array of spiral artery nodes and decidual vein nodes. Spiral artery nodes are mapped with stem villi.

       Inputs:
         - ellipsoid_coor:coordinate of nodes of placental mesh
         - surfacenode:array of surface nodes
         - stem_file:txt file that described stem villi locations

       Outputs:
         - spiral_array: array of spiral artery nodes
         - decidual_array: array of decidual artery nodes
         - vesselnode: array of both spiral and decidual nodes
         - surfnode_ex_vessel: array of surface node excluding vessel nodes
    """
    sa_radius = 3.7 / 2.0
    dv_radius = sa_radius

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    xyList = np.zeros((len(surfacenode), 4))
    count = 0
    for i in range(0, len(surfacenode)):  # taking only x and y coordinates
        if ellipsoid_coor[surfacenode[i] - 1, 3] > 0:  # take if upper surface nodes as this is where vessele reside
            # location from upper surface nodes only
            xyList[count, 0] = ellipsoid_coor[surfacenode[i] - 1, 0]  # node number
            xyList[count, 1] = ellipsoid_coor[surfacenode[i] - 1, 1]  # x-coordinate
            xyList[count, 2] = ellipsoid_coor[surfacenode[i] - 1, 2]  # y-coordinate
            xyList[count, 3] = ellipsoid_coor[surfacenode[i] - 1, 3]  # z-coordinate

            count = count + 1

    xyList = xyList[0:count, :]

    surfnode_ex_vessel = np.copy(surfacenode)
    vesselnode_temp = np.vstack({tuple(row) for row in xyList})  # nodes that might be vessels

    # reading in the stem vessel to map the spiral artery location
    vessel_mapped_stem = [-9.822741e+00, 1.550285e+01]
    vessel_mapped_stem_v = [1.155144e+01, 1.435972e+01]

    spiral_array = np.zeros((len(xyList)), dtype=int)  # store the node nuber of spiral artery
    decidual_array = np.zeros((len(xyList)), dtype=int)  # store the node number of decidual vein

    check = ellipsoid_coor[:, 0:2]
    np.random.seed(0)
    sa_nodes = 0
    dv_nodes = 0
    for i in range(0, len(vessel_mapped_stem)):  # for each blood vessel,Cycle through to find closest nodes
        for nodeX in vesselnode_temp:
            distance = np.sqrt((vessel_mapped_stem[0] - nodeX[1]) ** 2 + (
                    vessel_mapped_stem[1] - nodeX[2]) ** 2)  # distance from the nodes
            if (distance < sa_radius):
                #print('SA Node', int(nodeX[0]))
                arterynode = nodeX[0]
                A = np.where(vesselnode_temp == arterynode)
                vesselnode_temp = np.delete(vesselnode_temp, A[0], axis=0)
                A2 = np.where(surfnode_ex_vessel == int(arterynode))
                surfnode_ex_vessel = np.delete(surfnode_ex_vessel, A2)
                spiral_array[sa_nodes] = arterynode
                sa_nodes = sa_nodes + 1

    # Doing decidual veins after arteries to make sure we dont take up any spots that arteries would have otherwise beein
    for i in range(0, len(vessel_mapped_stem_v)):  # need same number of arteries as veins
        for nodeX in vesselnode_temp:
            distance = np.sqrt((vessel_mapped_stem_v[0] - nodeX[1]) ** 2 + (
                    vessel_mapped_stem_v[1] - nodeX[2]) ** 2)  # distance from the nodes
            if (distance < dv_radius):
                #print('DV Node', int(nodeX[0]))
                veinnode = nodeX[0]
                V = np.where(vesselnode_temp == veinnode)
                vesselnode_temp = np.delete(vesselnode_temp, V[0], axis=0)
                V2 = np.where(surfnode_ex_vessel == int(veinnode))
                surfnode_ex_vessel = np.delete(surfnode_ex_vessel, V2)
                decidual_array[dv_nodes] = veinnode
                dv_nodes = dv_nodes + 1


    spiral_array = np.resize(spiral_array, sa_nodes)
    decidual_array = np.resize(decidual_array, dv_nodes)
    #print('dec', decidual_array)

    return {'spiral_array': spiral_array, 'decidual_array': decidual_array, 'surfnode_ex_vessel': surfnode_ex_vessel}


def gen_3d_ellipsoid_structured(size_el, volume, thickness, ellipticity, squareSizeRatio, circle_prop, el_type, debug):
    """ Generates a structured ellipsoid mesh to solve 3D problems. The aim is for a quality computational mesh that
    has as regular elements as possible, within the constraints of typical dimensions of ellipsoids representing the
    volume of the placenta. This code is derived from an openCMISS example written by Chris Bradley, which is used to
    simulate fluid structure interactions in a cylinder. Note that this hasn't been tested on linear elements


    Inputs:
       - size_el: approximate dimension of an element in each axis that we are aiming for
       - volume: volume of placental ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimension
       - squareSizeRatio: ratio of square in mesh cross-section to radius
       - circle_prop: proportion of ellipse in x-y that is made up by 'plate' of nodes and elements
       - debug (True or False) allows you to print certain statements to screen

    Returns:
       - placental_node_coor: nodes location of mesh
       - placental_el_con: element connectivity of mesh (tetrahedral element)
       - node_array: array of nodes
       - element_array: array of elements
    """
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    ellipsoidRadius_z = radii['z_radius']
    ellipsoidRadius_x = radii['x_radius']
    ellipsoidRadius_y = radii['y_radius']

    if (debug):
        print('Solving a model with x-radius: ' + str(ellipsoidRadius_x) + ' y-radius: ' + str(
            ellipsoidRadius_y) + 'z-radius: ' + str(ellipsoidRadius_z))

    nel_x = int(np.floor((ellipsoidRadius_x * 2) / size_el))  # number of elems needed in x aixs in mesh
    nel_y = int(np.floor((ellipsoidRadius_y * 2) / size_el))  # number of elems needed in  in y axis in mesh (need to
    #  implement having different x,y element numbers
    nel_z = int(np.floor((ellipsoidRadius_z * 2) / size_el))  # number of elems needed in  in z axis in  mesh
    # total number of elements in x,y are number in square plus 2* number in arm
    # If square takes up half the radius need even numbers in arm and square at one third of total each
    # If square takes up squareSizeRatio of the total, then need the square part to be multiplied by that proportion
    total_square_arm = 2.0 * nel_x / 3.0
    numberOfSquareElements = int(np.floor(squareSizeRatio * total_square_arm))
    numberOfArmElements = int(np.floor((1 - squareSizeRatio) * total_square_arm))
    # In future for cross-sections that deviate a lot from circular will need different number of elements in square
    # and arm in x- and y-
    numberOfZElements = nel_z

    if (el_type == 1):  # linear
        numberOfNodesXi = 2
    elif (el_type == 2):  # quadratic
        numberOfNodesXi = 3

    numberOfLocalNodes = numberOfNodesXi * numberOfNodesXi * numberOfNodesXi
    numberOfLocalInterfaceNodes = numberOfNodesXi * numberOfNodesXi
    localNodeIdx000 = 0
    localNodeIdx100 = numberOfNodesXi - 1
    localNodeIdx010 = numberOfNodesXi * (numberOfNodesXi - 1)
    localNodeIdx110 = numberOfNodesXi * numberOfNodesXi - 1
    localNodeIdx001 = numberOfNodesXi * numberOfNodesXi * (numberOfNodesXi - 1)
    localNodeIdx101 = numberOfNodesXi - 1 + numberOfNodesXi * numberOfNodesXi * (numberOfNodesXi - 1)
    localNodeIdx011 = numberOfNodesXi * (numberOfNodesXi - 1) + numberOfNodesXi * numberOfNodesXi * (
            numberOfNodesXi - 1)
    localNodeIdx111 = numberOfLocalNodes - 1

    numberOfNodesPerBlock = numberOfSquareElements * (numberOfNodesXi - 1) * (
                numberOfArmElements * (numberOfNodesXi - 1) + 1)
    numberOfElementsPerBlock = numberOfSquareElements * numberOfArmElements
    numberOfNodesPerLength = 4 * numberOfNodesPerBlock + \
                             (numberOfSquareElements * (numberOfNodesXi - 1) - 1) * (
                                         numberOfSquareElements * (numberOfNodesXi - 1) - 1)
    numberOfElementsPerLength = 4 * numberOfElementsPerBlock + numberOfSquareElements * numberOfSquareElements
    numberOfNodes = numberOfNodesPerLength * (numberOfZElements * (numberOfNodesXi - 1) + 1)
    numberOfElements = numberOfElementsPerLength * numberOfZElements

    if debug:
        print('  Mesh Parameters:')
        print('    numberOfSquareElements: %d' % (numberOfSquareElements))
        print('    numberOfArmElements: %d' % (numberOfArmElements))
        print('    numberOfZElements: %d' % (numberOfZElements))
        print('    numberOfNodesXi: %d' % (numberOfNodesXi))
        print('    numberOfNodesPerBlock: %d' % (numberOfNodesPerBlock))
        print('    numberOfElementPerBlock: %d' % (numberOfElementsPerBlock))
        print('    numberOfNodesPerLength: %d' % (numberOfNodesPerLength))
        print('    numberOfElementsPerLength: %d' % (numberOfElementsPerLength))
        print('    numberOfNodes: %d' % (numberOfNodes))
        print('    numberOfElements: %d' % (numberOfElements))
        print('    numberOfLocalNodes: %d' % (numberOfLocalNodes))


    elems = np.zeros((numberOfElements, numberOfLocalNodes+1), dtype='int32')
    node_array = np.zeros((numberOfNodes,4))
    nodelist = [0]*numberOfNodes
    surface_nodes = [0] * numberOfNodes
    num_surface_nodes = 0

    for zElementIdx in range(1, max(numberOfZElements + 1, 2)):
        # Handle the arm blocks first
        previousBlock = 4
        for blockIdx in range(1, 5):  # generating arm blocks
            # DEFINING NODES AND ELEMENTS WITH CONNECTIVITY
            for yElementIdx in range(1, numberOfArmElements + 1):
                for xElementIdx in range(1, numberOfSquareElements + 1):
                    localNodes = [0] * numberOfLocalNodes  # Nodes local to this arm block
                    elementNumber = xElementIdx + (yElementIdx - 1) * numberOfSquareElements + (
                            blockIdx - 1) * numberOfSquareElements * numberOfArmElements + \
                                    (zElementIdx - 1) * numberOfElementsPerLength
                    if (xElementIdx == 1):
                        localNodes[localNodeIdx000] = (
                                                              previousBlock - 1) * numberOfNodesPerBlock + numberOfSquareElements * (
                                                              numberOfNodesXi - 1) + \
                                                      (yElementIdx - 1) * (
                                                              numberOfNodesXi - 1) * numberOfSquareElements * (
                                                              numberOfNodesXi - 1) + \
                                                      (zElementIdx - 1) * numberOfNodesPerLength * (numberOfNodesXi - 1)
                        localNodes[localNodeIdx100] = (blockIdx - 1) * numberOfNodesPerBlock + numberOfNodesXi - 1 + \
                                                      (yElementIdx - 1) * (
                                                              numberOfNodesXi - 1) * numberOfSquareElements * (
                                                              numberOfNodesXi - 1) + \
                                                      (zElementIdx - 1) * numberOfNodesPerLength * (numberOfNodesXi - 1)
                    else:
                        localNodes[localNodeIdx000] = (blockIdx - 1) * numberOfNodesPerBlock + (xElementIdx - 2) * (
                                numberOfNodesXi - 1) + (numberOfNodesXi - 2) + 1 + \
                                                      (yElementIdx - 1) * (numberOfNodesXi - 1) * (
                                                              numberOfSquareElements * (numberOfNodesXi - 1)) + \
                                                      (zElementIdx - 1) * numberOfNodesPerLength * (numberOfNodesXi - 1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000] + numberOfNodesXi - 1
                    localNodes[localNodeIdx010] = localNodes[localNodeIdx000] + numberOfSquareElements * (
                            numberOfNodesXi - 1) * (numberOfNodesXi - 1)
                    localNodes[localNodeIdx110] = localNodes[localNodeIdx100] + numberOfSquareElements * (
                            numberOfNodesXi - 1) * (numberOfNodesXi - 1)
                    localNodes[localNodeIdx001] = localNodes[localNodeIdx000] + numberOfNodesPerLength * (
                            numberOfNodesXi - 1)
                    localNodes[localNodeIdx101] = localNodes[localNodeIdx100] + numberOfNodesPerLength * (
                            numberOfNodesXi - 1)
                    localNodes[localNodeIdx011] = localNodes[localNodeIdx010] + numberOfNodesPerLength * (
                            numberOfNodesXi - 1)
                    localNodes[localNodeIdx111] = localNodes[localNodeIdx110] + numberOfNodesPerLength * (
                            numberOfNodesXi - 1)
                    if (el_type == 2):
                        localNodes[1] = localNodes[localNodeIdx100] - 1
                        localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements * (numberOfNodesXi - 1)
                        localNodes[4] = localNodes[1] + numberOfSquareElements * (numberOfNodesXi - 1)
                        localNodes[5] = localNodes[4] + 1
                        localNodes[7] = localNodes[localNodeIdx110] - 1
                        localNodes[9] = localNodes[0] + numberOfNodesPerLength
                        localNodes[10] = localNodes[1] + numberOfNodesPerLength
                        localNodes[11] = localNodes[2] + numberOfNodesPerLength
                        localNodes[12] = localNodes[3] + numberOfNodesPerLength
                        localNodes[13] = localNodes[4] + numberOfNodesPerLength
                        localNodes[14] = localNodes[5] + numberOfNodesPerLength
                        localNodes[15] = localNodes[6] + numberOfNodesPerLength
                        localNodes[16] = localNodes[7] + numberOfNodesPerLength
                        localNodes[17] = localNodes[8] + numberOfNodesPerLength
                        localNodes[19] = localNodes[10] + numberOfNodesPerLength
                        localNodes[21] = localNodes[12] + numberOfNodesPerLength
                        localNodes[22] = localNodes[13] + numberOfNodesPerLength
                        localNodes[23] = localNodes[14] + numberOfNodesPerLength
                        localNodes[25] = localNodes[16] + numberOfNodesPerLength
                    linearNodes = [localNodes[localNodeIdx000], localNodes[localNodeIdx100],
                                   localNodes[localNodeIdx010], localNodes[localNodeIdx110], \
                                   localNodes[localNodeIdx001], localNodes[localNodeIdx101],
                                   localNodes[localNodeIdx011], localNodes[localNodeIdx111]]
                    if (debug):
                        print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (elementNumber, linearNodes[0], linearNodes[1], linearNodes[2], linearNodes[3],
                               linearNodes[4], linearNodes[5], linearNodes[6], linearNodes[7]))
                        if (el_type == 2):
                            print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[0], localNodes[1], localNodes[2], localNodes[3], localNodes[4],
                                   localNodes[5], localNodes[6], localNodes[7], localNodes[8]))
                            print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[9], localNodes[10], localNodes[11], localNodes[12], localNodes[13],
                                   localNodes[14], localNodes[15], localNodes[16], localNodes[17]))
                            print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[18], localNodes[19], localNodes[20], localNodes[21], localNodes[22],
                                   localNodes[23], localNodes[24], localNodes[25], localNodes[26]))

                    if (el_type == 1):
                        elems[elementNumber-1][0] = elementNumber
                        elems[elementNumber-1][1:numberOfLocalNodes+1] = linearNodes
                    if (el_type == 2):
                        elems[elementNumber - 1][0] = elementNumber
                        elems[elementNumber - 1][1:numberOfLocalNodes + 1] = localNodes
            previousBlock = blockIdx
            # Handle the square block
        if (numberOfSquareElements == 1):
            elementNumber = elementNumber + 1
            localNodes[localNodeIdx000] = 3 * numberOfNodesPerBlock + \
                                          (zElementIdx - 1) * numberOfNodesPerLength * (numberOfNodesXi - 1)
            localNodes[localNodeIdx100] = 4 * numberOfNodesPerBlock + \
                                          (zElementIdx - 1) * numberOfNodesPerLength * (numberOfNodesXi - 1)
            localNodes[localNodeIdx010] = 2 * numberOfNodesPerBlock + \
                                          (zElementIdx - 1) * numberOfNodesPerLength * (numberOfNodesXi - 1)
            localNodes[localNodeIdx110] = numberOfNodesPerBlock + \
                                          (zElementIdx - 1) * numberOfNodesPerLength * (numberOfNodesXi - 1)
            if (el_type == 2):
                localNodes[1] = localNodes[localNodeIdx100] - 1
                localNodes[3] = localNodes[localNodeIdx000] - 1
                localNodes[4] = localNodes[localNodeIdx100] + 1
                localNodes[5] = localNodes[localNodeIdx110] - 1
                localNodes[7] = localNodes[localNodeIdx010] - 1
            localNodes[localNodeIdx001] = localNodes[localNodeIdx000] + numberOfNodesPerLength * (numberOfNodesXi - 1)
            localNodes[localNodeIdx101] = localNodes[localNodeIdx100] + numberOfNodesPerLength * (numberOfNodesXi - 1)
            localNodes[localNodeIdx011] = localNodes[localNodeIdx010] + numberOfNodesPerLength * (numberOfNodesXi - 1)
            localNodes[localNodeIdx111] = localNodes[localNodeIdx110] + numberOfNodesPerLength * (numberOfNodesXi - 1)
            linearNodes = [localNodes[localNodeIdx000], localNodes[localNodeIdx100], localNodes[localNodeIdx010],
                           localNodes[localNodeIdx110], \
                           localNodes[localNodeIdx001], localNodes[localNodeIdx101], localNodes[localNodeIdx011],
                           localNodes[localNodeIdx111]]
            if (el_type == 2):
                localNodes[9] = localNodes[0] + numberOfNodesPerLength
                localNodes[10] = localNodes[1] + numberOfNodesPerLength
                localNodes[11] = localNodes[2] + numberOfNodesPerLength
                localNodes[12] = localNodes[3] + numberOfNodesPerLength
                localNodes[13] = localNodes[4] + numberOfNodesPerLength
                localNodes[14] = localNodes[5] + numberOfNodesPerLength
                localNodes[15] = localNodes[6] + numberOfNodesPerLength
                localNodes[16] = localNodes[7] + numberOfNodesPerLength
                localNodes[17] = localNodes[8] + numberOfNodesPerLength
                localNodes[19] = localNodes[10] + numberOfNodesPerLength
                localNodes[21] = localNodes[12] + numberOfNodesPerLength
                localNodes[22] = localNodes[13] + numberOfNodesPerLength
                localNodes[23] = localNodes[14] + numberOfNodesPerLength
                localNodes[25] = localNodes[16] + numberOfNodesPerLength
            if (debug):
                print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (elementNumber, linearNodes[0], linearNodes[1], linearNodes[2], linearNodes[3], linearNodes[4],
                       linearNodes[5], linearNodes[6], linearNodes[7]))
                if (el_type == 2):
                    print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (localNodes[0], localNodes[1], localNodes[2], localNodes[3], localNodes[4], localNodes[5],
                           localNodes[6], localNodes[7], localNodes[8]))
                    print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (
                          localNodes[9], localNodes[10], localNodes[11], localNodes[12], localNodes[13], localNodes[14],
                          localNodes[15], localNodes[16], localNodes[17]))
                    print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (localNodes[18], localNodes[19], localNodes[20], localNodes[21], localNodes[22],
                           localNodes[23], localNodes[24], localNodes[25], localNodes[26]))

            if (el_type == 1):
                elems[elementNumber - 1][0] = elementNumber
                elems[elementNumber - 1][1:numberOfLocalNodes + 1] = linearNodes
            if (el_type == 2):
                elems[elementNumber - 1][0] = elementNumber
                elems[elementNumber - 1][1:numberOfLocalNodes + 1] = localNodes
        else:
            for yElementIdx in range(1, numberOfSquareElements + 1):
                for xElementIdx in range(1, numberOfSquareElements + 1):
                    localNodes = [0] * numberOfLocalNodes
                    elementNumber = 4 * numberOfElementsPerBlock + xElementIdx + (
                            yElementIdx - 1) * numberOfSquareElements + \
                                    (zElementIdx - 1) * numberOfElementsPerLength
                    if (yElementIdx == 1):
                        if (xElementIdx == 1):
                            # Bottom-left
                            localNodes[localNodeIdx000] = 3 * numberOfNodesPerBlock + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = 3 * numberOfNodesPerBlock + numberOfArmElements * (
                                        numberOfNodesXi - 1) * \
                                                          numberOfSquareElements * (numberOfNodesXi - 1) + (
                                                                      numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = 3 * numberOfNodesPerBlock - (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          (numberOfNodesXi - 2) + (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx100] - 1
                                localNodes[3] = localNodes[localNodeIdx000] - 1
                                localNodes[4] = localNodes[localNodeIdx110] - numberOfSquareElements * (
                                            numberOfNodesXi - 1)
                                localNodes[5] = localNodes[4] + 1
                                localNodes[7] = localNodes[localNodeIdx110] - 1
                        elif (xElementIdx == numberOfSquareElements):
                            # Bottom-right
                            localNodes[localNodeIdx000] = 4 * numberOfNodesPerBlock - (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = 4 * numberOfNodesPerBlock + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          (numberOfNodesXi - 1) - (numberOfNodesXi - 2) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = numberOfSquareElements * (numberOfNodesXi - 1) * \
                                                          numberOfArmElements * (numberOfNodesXi - 1) + (
                                                                      numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx000] + 1
                                localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements * (
                                            numberOfNodesXi - 1) + 1
                                localNodes[4] = localNodes[3] + 1
                                localNodes[5] = localNodes[localNodeIdx110] - 1
                                localNodes[7] = localNodes[localNodeIdx010] + 1
                        else:
                            # Bottom
                            localNodes[localNodeIdx000] = 3 * numberOfNodesPerBlock + numberOfSquareElements * (
                                        numberOfNodesXi - 1) * \
                                                          numberOfArmElements * (numberOfNodesXi - 1) + (
                                                                      xElementIdx - 1) * (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = localNodes[localNodeIdx000] + (numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          (numberOfNodesXi - 2) + (xElementIdx - 1) * (
                                                                      numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = localNodes[localNodeIdx010] + (numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx000] + 1
                                localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements * (
                                            numberOfNodesXi - 1) + 1
                                localNodes[4] = localNodes[3] + 1
                                localNodes[5] = localNodes[4] + 1
                                localNodes[7] = localNodes[localNodeIdx110] - 1
                    elif (yElementIdx == numberOfSquareElements):
                        if (xElementIdx == 1):
                            # Top-left
                            localNodes[localNodeIdx000] = 2 * numberOfNodesPerBlock + numberOfSquareElements * (
                                        numberOfNodesXi - 1) * \
                                                          numberOfArmElements * (numberOfNodesXi - 1) + (
                                                                      numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          ((numberOfSquareElements - 1) * (numberOfNodesXi - 1) - 1) + (
                                                                      numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = 2 * numberOfNodesPerBlock + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = 2 * numberOfNodesPerBlock - (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx100] - 1
                                localNodes[3] = localNodes[localNodeIdx000] - 1
                                localNodes[4] = localNodes[1] + numberOfSquareElements * (numberOfNodesXi - 1) - 1
                                localNodes[5] = localNodes[4] + 1
                                localNodes[7] = localNodes[localNodeIdx110] + 1
                        elif (xElementIdx == numberOfSquareElements):
                            # Top-right
                            localNodes[localNodeIdx000] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          ((numberOfSquareElements - 1) * (numberOfNodesXi - 1) - 1) + \
                                                          (numberOfSquareElements - 1) * (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = numberOfSquareElements * (numberOfNodesXi - 1) * \
                                                          numberOfArmElements * (numberOfNodesXi - 1) + \
                                                          (numberOfSquareElements - 1) * (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = numberOfNodesPerBlock + numberOfSquareElements * (
                                        numberOfNodesXi - 1) * \
                                                          numberOfArmElements * (numberOfNodesXi - 1) + (
                                                                      numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = numberOfNodesPerBlock + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx000] + 1
                                localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements * (
                                            numberOfNodesXi - 1) - 1
                                localNodes[4] = localNodes[3] + 1
                                localNodes[5] = localNodes[localNodeIdx110] - 1
                                localNodes[7] = localNodes[localNodeIdx010] - 1
                        else:
                            # Top
                            localNodes[localNodeIdx000] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          ((numberOfSquareElements - 1) * (numberOfNodesXi - 1) - 1) + \
                                                          (xElementIdx - 1) * (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = localNodes[localNodeIdx000] + (numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = 2 * numberOfNodesPerBlock - (xElementIdx - 1) * (
                                        numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = localNodes[localNodeIdx010] - (numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx000] + 1
                                localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements * (
                                            numberOfNodesXi - 1) - 1
                                localNodes[4] = localNodes[3] + 1
                                localNodes[5] = localNodes[4] + 1
                                localNodes[7] = localNodes[localNodeIdx010] - 1
                    else:
                        if (xElementIdx == 1):
                            # Left
                            localNodes[localNodeIdx000] = 3 * numberOfNodesPerBlock - (yElementIdx - 1) * (
                                        numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          ((yElementIdx - 1) * (numberOfNodesXi - 1) - 1) + (
                                                                      numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = localNodes[localNodeIdx000] - (numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = localNodes[localNodeIdx100] + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          (numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx100] - 1
                                localNodes[3] = localNodes[localNodeIdx000] - 1
                                localNodes[4] = localNodes[localNodeIdx110] - numberOfSquareElements * (
                                            numberOfNodesXi - 1)
                                localNodes[5] = localNodes[4] + 1
                                localNodes[7] = localNodes[localNodeIdx110] - 1
                        elif (xElementIdx == numberOfSquareElements):
                            # Right
                            localNodes[localNodeIdx000] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          ((yElementIdx - 1) * (numberOfNodesXi - 1) - 1) + (
                                                                      numberOfSquareElements - 1) * (
                                                                      numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = numberOfSquareElements * (
                                        numberOfNodesXi - 1) * numberOfArmElements * (numberOfNodesXi - 1) + \
                                                          (yElementIdx - 1) * (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = localNodes[localNodeIdx000] + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          (numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = localNodes[localNodeIdx100] + (numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx000] + 1
                                localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements * (
                                            numberOfNodesXi - 1) + 1
                                localNodes[4] = localNodes[3] + 1
                                localNodes[5] = localNodes[localNodeIdx100] + 1
                                localNodes[7] = localNodes[localNodeIdx010] + 1
                        else:
                            # Middle
                            localNodes[localNodeIdx000] = 4 * numberOfNodesPerBlock + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          ((yElementIdx - 1) * (numberOfNodesXi - 1) - 1) + (
                                                                      xElementIdx - 1) * (numberOfNodesXi - 1) + \
                                                          (zElementIdx - 1) * numberOfNodesPerLength * (
                                                                      numberOfNodesXi - 1)
                            localNodes[localNodeIdx100] = localNodes[localNodeIdx000] + (numberOfNodesXi - 1)
                            localNodes[localNodeIdx010] = localNodes[localNodeIdx000] + (
                                        numberOfSquareElements * (numberOfNodesXi - 1) - 1) * \
                                                          (numberOfNodesXi - 1)
                            localNodes[localNodeIdx110] = localNodes[localNodeIdx010] + (numberOfNodesXi - 1)
                            if (el_type == 2):
                                localNodes[1] = localNodes[localNodeIdx000] + 1
                                localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements * (
                                            numberOfNodesXi - 1) - 1
                                localNodes[4] = localNodes[3] + 1
                                localNodes[5] = localNodes[4] + 1
                                localNodes[7] = localNodes[localNodeIdx010] + 1
                    localNodes[localNodeIdx001] = localNodes[localNodeIdx000] + numberOfNodesPerLength * (
                                numberOfNodesXi - 1)
                    localNodes[localNodeIdx101] = localNodes[localNodeIdx100] + numberOfNodesPerLength * (
                                numberOfNodesXi - 1)
                    localNodes[localNodeIdx011] = localNodes[localNodeIdx010] + numberOfNodesPerLength * (
                                numberOfNodesXi - 1)
                    localNodes[localNodeIdx111] = localNodes[localNodeIdx110] + numberOfNodesPerLength * (
                                numberOfNodesXi - 1)
                    linearNodes = [localNodes[localNodeIdx000], localNodes[localNodeIdx100],
                                   localNodes[localNodeIdx010], localNodes[localNodeIdx110], \
                                   localNodes[localNodeIdx001], localNodes[localNodeIdx101],
                                   localNodes[localNodeIdx011], localNodes[localNodeIdx111]]
                    if (el_type == 2):
                        localNodes[9] = localNodes[0] + numberOfNodesPerLength
                        localNodes[10] = localNodes[1] + numberOfNodesPerLength
                        localNodes[11] = localNodes[2] + numberOfNodesPerLength
                        localNodes[12] = localNodes[3] + numberOfNodesPerLength
                        localNodes[13] = localNodes[4] + numberOfNodesPerLength
                        localNodes[14] = localNodes[5] + numberOfNodesPerLength
                        localNodes[15] = localNodes[6] + numberOfNodesPerLength
                        localNodes[16] = localNodes[7] + numberOfNodesPerLength
                        localNodes[17] = localNodes[8] + numberOfNodesPerLength
                        localNodes[19] = localNodes[10] + numberOfNodesPerLength
                        localNodes[21] = localNodes[12] + numberOfNodesPerLength
                        localNodes[22] = localNodes[13] + numberOfNodesPerLength
                        localNodes[23] = localNodes[14] + numberOfNodesPerLength
                        localNodes[25] = localNodes[16] + numberOfNodesPerLength
                    if (debug):
                        print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (elementNumber, linearNodes[0], linearNodes[1], linearNodes[2], linearNodes[3],
                               linearNodes[4], linearNodes[5], linearNodes[6], linearNodes[7]))
                        if (el_type == 2):
                            print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[0], localNodes[1], localNodes[2], localNodes[3], localNodes[4],
                                   localNodes[5], localNodes[6], localNodes[7], localNodes[8]))
                            print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[9], localNodes[10], localNodes[11], localNodes[12], localNodes[13],
                                   localNodes[14], localNodes[15], localNodes[16], localNodes[17]))
                            print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[18], localNodes[19], localNodes[20], localNodes[21], localNodes[22],
                                   localNodes[23], localNodes[24], localNodes[25], localNodes[26]))
                    if (el_type == 1):
                        elems[elementNumber-1][0] = elementNumber
                        elems[elementNumber-1][1:numberOfLocalNodes+1] = linearNodes
                    if (el_type == 2):
                        elems[elementNumber - 1][0] = elementNumber
                        elems[elementNumber - 1][1:numberOfLocalNodes + 1] = localNodes

    if (debug):
        print('   Nodes:')
    for zNodeIdx in range(1, numberOfZElements * (numberOfNodesXi - 1) + 2):
        prop = 1 - (1 - circle_prop) * abs(2.0 * (zNodeIdx - 1) / float(numberOfZElements * (numberOfNodesXi - 1)) -
                                           1.0)
        sign = np.sign(2.0 * (zNodeIdx - 1) / float(numberOfZElements * (numberOfNodesXi - 1)) - 1.0)
        zPosition = sign * ellipsoidRadius_z * np.sqrt(1 - prop ** 2)

        new_x_radius = ellipsoidRadius_x * np.sqrt(ellipsoidRadius_z ** 2.0 - zPosition ** 2.0) \
                       / (ellipsoidRadius_z)
        new_y_radius = ellipsoidRadius_y * np.sqrt(ellipsoidRadius_z ** 2.0 - zPosition ** 2.0) \
                       / (ellipsoidRadius_z)

        angle_theta = np.arctan(new_y_radius / new_x_radius)

        squareSize_x = squareSizeRatio * new_x_radius * np.cos(angle_theta)
        squareSize_y = squareSizeRatio * new_y_radius * np.sin(angle_theta)

        # Handle the arm blocks first
        previousBlock = 4
        for blockIdx in range(1, 5):
            # print('Block which ' + str(blockIdx) + ' ' + str(zNodeIdx))
            for yNodeIdx in range(1, numberOfArmElements * (numberOfNodesXi - 1) + 2):
                for xNodeIdx in range(1, numberOfSquareElements * (numberOfNodesXi - 1) + 1):
                    nodeNumber = (blockIdx - 1) * numberOfNodesPerBlock + xNodeIdx + (
                                yNodeIdx - 1) * numberOfSquareElements * (numberOfNodesXi - 1) + \
                                 (zNodeIdx - 1) * numberOfNodesPerLength
                    #nodeDomain = decomposition.NodeDomainGet(nodeNumber, 1)
                    if(nodeNumber > numberOfNodes + 1):
                        print(nodeNumber)
                    else:
                        if (yNodeIdx == numberOfArmElements * (numberOfNodesXi - 1) + 1):
                            # On the square
                            # print('On the square', xNodeIdx,yNodeIdx)

                            if (blockIdx == 1):
                                xPosition = squareSize_x - 2.0 * xNodeIdx * squareSize_x / (
                                            numberOfSquareElements * (numberOfNodesXi - 1))
                                yPosition = squareSize_y
                            elif (blockIdx == 2):
                                xPosition = -squareSize_x
                                yPosition = squareSize_y - 2.0 * xNodeIdx * squareSize_y / (
                                            numberOfSquareElements * (numberOfNodesXi - 1))
                            elif (blockIdx == 3):
                                xPosition = -squareSize_x + 2.0 * xNodeIdx * squareSize_x / (
                                            numberOfSquareElements * (numberOfNodesXi - 1))
                                yPosition = -squareSize_y
                            elif (blockIdx == 4):
                                xPosition = squareSize_x
                                yPosition = -squareSize_y + 2.0 * xNodeIdx * squareSize_y / (
                                            numberOfSquareElements * (numberOfNodesXi - 1))
                        else:
                            # In the arm
                            # Work out the r, theta position of each point equally spread on the block
                            if (blockIdx == 1):
                                start_theta = np.arctan(new_y_radius / new_x_radius)
                                end_theta = np.pi - start_theta
                                theta = start_theta + xNodeIdx * (end_theta - start_theta) / (numberOfSquareElements * (
                                        numberOfNodesXi - 1))
                                # theta is the angle from the centre of the mesh to the surface of the ellipsoid
                                sq_x = squareSize_x - 2.0 * xNodeIdx * squareSize_x / (
                                        numberOfSquareElements * (numberOfNodesXi - 1))
                                sq_y = squareSize_y
                            elif (blockIdx == 2):
                                start_theta = np.pi - np.arctan(new_y_radius / new_x_radius)
                                end_theta = np.pi + np.arctan(new_y_radius / new_x_radius)
                                theta = start_theta + xNodeIdx * (end_theta - start_theta) / (numberOfSquareElements * (
                                        numberOfNodesXi - 1))
                                sq_x = -squareSize_x
                                sq_y = squareSize_y - 2.0 * xNodeIdx * squareSize_y / (
                                        numberOfSquareElements * (numberOfNodesXi - 1))
                            elif (blockIdx == 3):
                                start_theta = np.pi + np.arctan(new_y_radius / new_x_radius)
                                end_theta = 2.0 * np.pi - np.arctan(new_y_radius / new_x_radius)
                                theta = start_theta + xNodeIdx * (end_theta - start_theta) / (numberOfSquareElements * (
                                        numberOfNodesXi - 1))
                                sq_x = -squareSize_x + 2.0 * xNodeIdx * squareSize_x / (
                                            numberOfSquareElements * (numberOfNodesXi - 1))
                                sq_y = -squareSize_y
                            elif (blockIdx == 4):
                                start_theta = 2.0 * np.pi - np.arctan(new_y_radius / new_x_radius)
                                end_theta = 2.0 * np.pi + np.arctan(new_y_radius / new_x_radius)
                                theta = start_theta + xNodeIdx * (end_theta - start_theta) / (numberOfSquareElements * (
                                        numberOfNodesXi - 1))
                                sq_x = squareSize_x
                                sq_y = -squareSize_y + 2.0 * xNodeIdx * squareSize_y / (
                                        numberOfSquareElements * (numberOfNodesXi - 1))

                            armRadius = new_y_radius * new_x_radius / np.sqrt(
                                new_x_radius ** 2 * np.sin(theta) ** 2 + new_y_radius ** 2 * np.cos(theta) ** 2)
                            arm_x = armRadius * np.cos(theta)
                            arm_y = armRadius * np.sin(theta)

                            arm_no = (yNodeIdx - 1) / (numberOfArmElements * (numberOfNodesXi - 1) + 1.0)

                            xPosition = arm_x - arm_no * (arm_x - sq_x)
                            yPosition = arm_y - arm_no * (arm_y - sq_y)

                        if (zNodeIdx == 1):
                            zPosition = -ellipsoidRadius_z * np.sqrt(
                                1 - xPosition ** 2 / ellipsoidRadius_x ** 2 - yPosition ** 2 / ellipsoidRadius_y ** 2)
                            num_surface_nodes = num_surface_nodes + 1
                            surface_nodes[num_surface_nodes] = nodeNumber
                        elif (zNodeIdx == numberOfZElements * (numberOfNodesXi - 1) + 1):
                            zPosition = ellipsoidRadius_z * np.sqrt(
                                1 - xPosition ** 2 / ellipsoidRadius_x ** 2 - yPosition ** 2 / ellipsoidRadius_y ** 2)
                            num_surface_nodes = num_surface_nodes + 1
                            surface_nodes[num_surface_nodes] = nodeNumber
                        elif(yNodeIdx == 1):
                            num_surface_nodes = num_surface_nodes + 1
                            surface_nodes[num_surface_nodes] = nodeNumber
                        nodelist[nodeNumber-1] = nodeNumber
                        node_array[nodeNumber-1][:]= [nodeNumber,xPosition,yPosition,zPosition]

                        if (debug):
                            print('      Node        %d:' % (nodeNumber))
                            print('         Position         = [ %.2f, %.2f, %.2f ]' % (
                            xPosition, yPosition, zPosition))

        # Now handle square
        for yNodeIdx in range(2, numberOfSquareElements * (numberOfNodesXi - 1) + 1):
            for xNodeIdx in range(2, numberOfSquareElements * (numberOfNodesXi - 1) + 1):
                nodeNumber = 4 * numberOfNodesPerBlock + (xNodeIdx - 1) + (yNodeIdx - 2) * (
                            numberOfSquareElements * (numberOfNodesXi - 1) - 1) + \
                             (zNodeIdx - 1) * numberOfNodesPerLength
                if (nodeNumber > numberOfNodes + 1):
                    print(nodeNumber)
                else:
                    xPosition = squareSize_x - squareSize_x * (yNodeIdx - 1) / (
                                numberOfSquareElements * (numberOfNodesXi - 1)) * 2.0
                    yPosition = -squareSize_y + squareSize_y * (xNodeIdx - 1) / (
                                numberOfSquareElements * (numberOfNodesXi - 1)) * 2.0

                    if (zNodeIdx == 1):
                        zPosition = -ellipsoidRadius_z * np.sqrt(
                            1 - xPosition ** 2 / ellipsoidRadius_x ** 2 - yPosition ** 2 / ellipsoidRadius_y ** 2)
                        num_surface_nodes = num_surface_nodes + 1
                        surface_nodes[num_surface_nodes] = nodeNumber
                    elif (zNodeIdx == numberOfZElements * (numberOfNodesXi - 1) + 1):
                        zPosition = ellipsoidRadius_z * np.sqrt(
                            1 - xPosition ** 2 / ellipsoidRadius_x ** 2 - yPosition ** 2 / ellipsoidRadius_y ** 2)
                        num_surface_nodes = num_surface_nodes + 1
                        surface_nodes[num_surface_nodes] = nodeNumber
                    nodelist[nodeNumber - 1] = nodeNumber
                    node_array[nodeNumber - 1][:] = [nodeNumber, xPosition, yPosition, zPosition]

                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition, yPosition, zPosition))


    surface_nodes = np.unique(surface_nodes)
    nzid = np.nonzero(surface_nodes)
    surface_nodes = surface_nodes[nzid]
    return{'nodes':node_array,'elems':elems,'surface_nodes':surface_nodes, 'node_list':nodelist}