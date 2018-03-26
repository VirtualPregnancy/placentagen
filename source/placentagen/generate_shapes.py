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


def umbilical_seed_geometry(volume, thickness, ellipticity, insertion_x, insertion_y, umb_artery_distance,
                            umb_artery_length):
    # Creating a basis for a branching geometry which assumes a simple umbilical cord structure, with two arteries,

    # Calulate axis dimensions of ellipsoid with given volume, thickness and ellipticity
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    # initialise node and element arrays
    node_loc = np.zeros((6, 4))
    elems = np.zeros((5, 3))
    elem_upstream = np.zeros((5, 3))
    elem_downstream = np.zeros((5, 3))

    # basic umbilical artery structure
    node_loc[0][0] = 0
    node_loc[0][1] = insertion_x
    node_loc[0][2] = insertion_y
    node_loc[0][3] = pg_utilities.z_from_xy(node_loc[0][1], node_loc[0][2], x_radius, y_radius,
                                            z_radius) + umb_artery_length + 3.0  # dummy branch 3mm long by default

    # node 2 is 3 mm up from node 1 in the z direction
    node_loc[1][0] = 1
    node_loc[1][1] = insertion_x
    node_loc[1][2] = insertion_y
    node_loc[1][3] = pg_utilities.z_from_xy(node_loc[0][1], node_loc[0][2], x_radius, y_radius,
                                            z_radius) + umb_artery_length

    # node 3 & 4 is the start of the 'umbilical artery'
    node_loc[2][0] = 2
    node_loc[2][1] = insertion_x
    node_loc[2][2] = insertion_y - umb_artery_distance / 2.0
    node_loc[2][3] = pg_utilities.z_from_xy(node_loc[0][1], node_loc[0][2], x_radius, y_radius,
                                            z_radius) + umb_artery_length
    node_loc[3][0] = 3
    node_loc[3][1] = insertion_x
    node_loc[3][2] = insertion_y + umb_artery_distance / 2.0
    node_loc[3][3] = pg_utilities.z_from_xy(node_loc[0][1], node_loc[0][2], x_radius, y_radius,
                                            z_radius) + umb_artery_length

    # node 5 and 6 'hit' the chorionic plate.
    node_loc[4][0] = 4
    node_loc[4][1] = insertion_x
    node_loc[4][2] = insertion_y - umb_artery_distance / 2.0
    node_loc[4][3] = pg_utilities.z_from_xy(node_loc[4][1], node_loc[4][2], x_radius, y_radius, z_radius)
    node_loc[5][0] = 5
    node_loc[5][1] = insertion_x
    node_loc[5][2] = insertion_y + umb_artery_distance / 2.0
    node_loc[5][3] = pg_utilities.z_from_xy(node_loc[5][1], node_loc[5][2], x_radius, y_radius, z_radius)
    # element locations
    elems[0, :] = [0, 0, 1]
    elem_upstream[0][0] = 0
    elem_downstream[0][0] = 2
    elem_downstream[0][1] = 1
    elem_downstream[0][2] = 2

    elems[1, :] = [1, 1, 2]
    elem_upstream[1][0] = 1
    elem_upstream[1][1] = 0
    elem_downstream[1][0] = 1
    elem_downstream[1][1] = 3

    elems[2, :] = [2, 1, 3]
    elem_upstream[2][0] = 1
    elem_upstream[2][1] = 0
    elem_downstream[2][0] = 1
    elem_downstream[2][1] = 4

    elems[3, :] = [3, 2, 4]
    elem_upstream[3][0] = 1
    elem_upstream[3][1] = 1
    elem_downstream[3][0] = 0

    elems[4, :] = [4, 3, 5]
    elem_upstream[4][0] = 1
    elem_upstream[4][1] = 2
    elem_downstream[4][0] = 0

    return {'umb_nodes': node_loc, 'umb_elems': elems, 'elem_up': elem_upstream, 'elem_down': elem_downstream}
