#!/usr/bin/env python
import numpy as np

from . import pg_utilities


# Contains code to generate placental shapes for generic placental models (i.e. from literature measures without
# specific data from an individual

def equispaced_data_in_ellipsoid(n, volume, thickness, ellipticity):
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
    # Generates equally spaced data points on the surface of an ellipsoid with the following inputs
    # n=number of data points which we aim to generate
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
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

def generate_rectangular_mesh(x_min,x_max,y_min,y_max,z_min,z_max,nel_x,nel_y,nel_z,x_width,y_width,z_width):
    x = np.linspace(x_min,x_max, x_max/x_width+1)#linspace for x axis
    y = np.linspace(y_min, y_max, y_max/y_width+1)#linspace for y axis
    z = np.linspace(z_min, z_max, z_max/z_width+1)#linspace for z axis
    nodes = np.vstack(np.meshgrid(y,z,x)).reshape(3,-1)#generate rectangular mesh

    y=np.array(nodes[0,:])#y coordinates of nodes in rectangular mesh
    z=np.array(nodes[1,:])#z coordinates of nodes in rectangular mesh
    x=np.array(nodes[2,:])#x coordinates of nodes in rectangular mesh

    #Generating the element connectivity of each cube element, 8 nodes for each 3D cube element
    nodeOfelement=np.zeros((8,nel_x*nel_y*nel_z))#this stores element array of each mesh grid element
    E1=0

    for K1 in range (1,nel_z+1):
       for J1 in range (1,nel_y+1):
          for I1 in range(1,nel_x+1):
           
            nodeOfelement[0,E1] = I1+(nel_x+1)*(J1-1)+(nel_x+1)*(nel_y+1)*(K1-1)#node 1 of 3D cube
            nodeOfelement[1,E1] = nodeOfelement[0,E1]+1#node 2 of 3D cube
            nodeOfelement[2,E1] = nodeOfelement[0,E1]+nel_x+1#node 3 of 3D cube
            nodeOfelement[3,E1] = nodeOfelement[2,E1]+1#node 4 of 3D cube
            nodeOfelement[4,E1] = nodeOfelement[0,E1]+(nel_x+1)*(nel_y+1)#node 5 of 3D cube
            nodeOfelement[5,E1] = nodeOfelement[1,E1]+(nel_x+1)*(nel_y+1)#node 6 of 3D cube
            nodeOfelement[6,E1] = nodeOfelement[2,E1]+(nel_x+1)*(nel_y+1)#node 7 of 3D cube
            nodeOfelement[7,E1] = nodeOfelement[3,E1]+(nel_x+1)*(nel_y+1)#node 8 of 3D cube
            
            E1 = E1+1

    nodeOfelement=nodeOfelement.T
    return {'nodeOfelement': nodeOfelement, 'x_coor': x,'y_coor':y,'z_coor':z,'total_mesh_el':nel_x*nel_y*nel_z,'total_mesh_node':(nel_x+1)*(nel_y+1)*(nel_z+1)}

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

    # Create linspaces for x y and z coordinates
    x = np.linspace(-x_width / 2.0, x_width / 2.0, nnod_x)  # linspace for x axis
    y = np.linspace(-y_width / 2.0, y_width / 2.0, nnod_y)  # linspace for y axis
    z = np.linspace(-z_width / 2.0, z_width / 2.0, nnod_z)  # linspace for z axis
    node_loc_temp = np.vstack(np.meshgrid(y, z, x)).reshape(3, -1).T  # generate nodes for rectangular mesh

    node_loc = np.zeros((nnod_x*nnod_y*nnod_z,3))
    for i in range(0,len(node_loc)):
        node_loc[i][0] = node_loc_temp[i][2]
        node_loc[i][1] = node_loc_temp[i][0]
        node_loc[i][2] = node_loc_temp[i][1]

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
