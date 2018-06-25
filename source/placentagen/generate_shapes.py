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


def gen_rectangular_mesh(volume, thickness, ellipticity, x_spacing, y_spacing, z_spacing):
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
    
    node_loc=gen_rectangular_node(x_width, y_width, z_width, nnod_x, nnod_y, nnod_z)
    # Generating the element connectivity of each cube element, 8 nodes for each 3D cube element
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
                elems[ne][1] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1) #lowest coordinates
                elems[ne][2] = elems[ne][1] + 1 #add one in x
                elems[ne][3] = elems[ne][1] + nnod_x #go through x and find first in y
                elems[ne][4] = elems[ne][3] + 1 #add one in y
                elems[ne][5] = elems[ne][1] + nnod_x * nnod_y #same as 1 -4 but at higher z -coord
                elems[ne][6] = elems[ne][2] + nnod_x * nnod_y
                elems[ne][7] = elems[ne][3] + nnod_x * nnod_y
                elems[ne][8] = elems[ne][4] + nnod_x * nnod_y
                ne = ne + 1

    return {'nodes': node_loc, 'elems': elems, 'total_nodes': nnod_x * nnod_y * nnod_z,
            'total_elems': (nnod_x - 1) * (nnod_y - 1) * (nnod_z - 1)}


def gen_mesh_darcy(volume,thickness,ellipticity,n):
    """ Generates ellipsoid placental mesh to solve darcy flow

    Inputs:
       - volume: volume of placental ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions
       - n: number of datapoints (optimal is 682000)

    Returns:
       - nodes: nodes location of mesh
       - elems: element connectivity of mesh (tetrahedral element)
       - node_array: array of nodes
       - element_array: array of elements"""

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    
    nodeSpacing = (n/(2*x_radius*2*y_radius*2*z_radius)) **(1./3)
        
    nnod_x=2*x_radius*nodeSpacing
    nnod_y=2*y_radius*nodeSpacing
    nnod_z=2*z_radius*nodeSpacing
    nodes=gen_rectangular_node(x_radius*2, y_radius*2, z_radius*2, nnod_x, nnod_y, nnod_z)
    
    #nodes inside the ellipsoid    
    ellipsoid_node=np.zeros((len(nodes),3))
    count=0
    for nnode in range (0, len(nodes)):
       coord_point = nodes[nnode][0:3]
       inside=pg_utilities.check_in_on_ellipsoid(coord_point[0], coord_point[1], coord_point[2], x_radius, y_radius, z_radius)
       if inside:
          ellipsoid_node[count,:]=coord_point[:]
          count=count+1
    ellipsoid_node.resize(count,3)     
    xyList = ellipsoid_node[:,[0,1]]
    xyListUnique = np.vstack({tuple(row) for row in xyList})
    #looking for z_coordinate of surface nodes
    for xyColumn in xyListUnique:
        
        xyNodes = np.where(np.all(xyList == xyColumn, axis = 1))[0]
        if len(xyNodes) > 1:
           x_coord=ellipsoid_node[xyNodes[0],0]
           y_coord=ellipsoid_node[xyNodes[0],1]
           ellipsoid_node[xyNodes[len(xyNodes) - 1],2] = pg_utilities.z_from_xy(x_coord, y_coord, x_radius, y_radius, z_radius)   
           ellipsoid_node[xyNodes[0],2] =-1*( pg_utilities.z_from_xy(x_coord, y_coord, x_radius, y_radius, z_radius))

    #generate tetrahedral mesh
    pyMesh = Delaunay(ellipsoid_node)

    #Build arrays to pass into openCMISS conversion:
    node_loc = pyMesh.points
    temp_elems = pyMesh.simplices
    #CHECK ELEMENTS FOR 0 VOLUME:
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
          
     vmat = np.vstack((x_coor,y_coor,z_coor,[1.0,1.0,1.0,1.0]))#matrix of coor of element
     elem_volume = (1/6.0) * abs(np.linalg.det(vmat))#volume of each tetrahedral element
     
     #if volume is not zero
     if elem_volume > min_vol:
       
       indexArr.append(index)
     index = index+1
     
    #update arrays without 0 volume elements, to pass into openCMISS
    elems = temp_elems[indexArr,:]
    for i in range(len(elems)):
       elems[i] = [x+1 for x in elems[i]]
    element_array = range(1, len(elems)+1)
    node_array = range(1, len(node_loc)+1)

    return {'nodes': node_loc, 'elems': elems, 'element_array':element_array,'node_array': node_array}



def gen_rectangular_node(x_width, y_width, z_width, nnod_x, nnod_y, nnod_z):
      
    # Create linspaces for x y and z coordinates
    x = np.linspace(-x_width / 2.0, x_width / 2.0, nnod_x)  # linspace for x axis
    y = np.linspace(-y_width / 2.0, y_width / 2.0, nnod_y)  # linspace for y axis
    z = np.linspace(-z_width / 2.0, z_width / 2.0, nnod_z)  # linspace for z axis
    node_loc_temp = np.vstack(np.meshgrid(y, z, x)).reshape(3, -1).T  # generate nodes for rectangular mesh
    node_loc = np.zeros((len(node_loc_temp),3))
    for i in range(0,len(node_loc)):
        node_loc[i][0] = node_loc_temp[i][2]
        node_loc[i][1] = node_loc_temp[i][0]
        node_loc[i][2] = node_loc_temp[i][1]
    
    return node_loc

