#!/usr/bin/env python
import numpy as np


def equispaced_data_in_ellipsoid(n, volume, thickness, ellipticity):
    # Generates equally sapced data points in an ellipsoid with the following inputs
    # n=number of data points which we aim to generate
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
    data_spacing = (volume / n) ** (1.0 / 3.0)
    radii = calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    # Aiming to generate seed points that fill a cuboid encompasing the placental volume then remove seed points that are
    # external to the ellipsoid

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
    count = 0
    for i in range(len(data_coords)):  # Loop through grid
        coord_check = (data_coords[i][0] / x_radius) ** 2 + (data_coords[i][1] / y_radius) ** 2 + (
                data_coords[i][2] / z_radius) ** 2  # check the point is in an ellipsoid
        if (coord_check) < 1.0:  # Has to be strictly in the ellipsoid
            Edata[count, :] = data_coords[i, :]  # add to data array
            count = count + 1
    Edata.resize(count, 3)  # resize data array to correct size

    print('Data points within ellipsoid allocated. Total = ' + str(len(Edata)))

    return Edata


def uniform_data_on_ellipsoid(n, volume, thickness, ellipticity, random_seed):
    # Generates equally spaced data points on the surface of an ellipsoid with the following inputs
    # n=number of data points which we aim to generate
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
    radii = calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    area_estimate = np.pi * x_radius * y_radius
    data_spacing = 0.85 * np.sqrt(area_estimate / n)

    chorion_data = np.zeros((n, 3))
    np.random.seed(random_seed)
    generated_seed = 0
    acceptable_attempts = n * 10  # try not to have too many failures
    attempts = 0

    while generated_seed < n and attempts < acceptable_attempts:
        # generate random x-y coordinates between -1 and 1
        random_number = np.random.uniform(-1.0, 1.0, 2)
        new_x = random_number[0] * x_radius
        new_y = random_number[1] * y_radius
        # check if new coordinate is on the ellipse
        if ((new_x / x_radius) ** 2 + (new_y / y_radius) ** 2) < 1:  # on the surface
            if (generated_seed) == 0:
                generated_seed = generated_seed + 1
                new_z = z_from_xy(new_x, new_y, x_radius, y_radius, z_radius)
                chorion_data[generated_seed - 1][:] = [new_x, new_y, new_z]
            else:
                reject = False
                for j in range(0, generated_seed):
                    distance = 0.0
                    for i in range(0, 3):
                        distance = distance + (chorion_data[j - 1][i] - new_x) ** 2
                    distance = np.sqrt(distance)
                    if distance <= data_spacing:
                        reject = True
                        break
                if reject == False:
                    generated_seed = generated_seed + 1
                    new_z = z_from_xy(new_x, new_y, x_radius, y_radius, z_radius)
                    chorion_data[generated_seed - 1][:] = [new_x, new_y, new_z]

        attempts = attempts + 1
    chorion_data.resize(generated_seed, 3)  # resize data array to correct size
    print('Data points on ellipsoid allocated. Total = ' + str(len(chorion_data)))

    return chorion_data


def umbilical_seed_geometry(volume, thickness, ellipticity, insertion_x, insertion_y, umb_artery_distance,
                            umb_artery_length):
    radii = calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    node_loc = np.zeros((6, 4))
    elems = np.zeros((5, 3))
    node_loc[0][0] = 1
    node_loc[0][1] = insertion_x
    node_loc[0][2] = insertion_y
    node_loc[0][3] = z_radius * np.sqrt(
        1.0 - (node_loc[0][1] / x_radius) ** 2 - (
                node_loc[0][2] / y_radius) ** 2) + umb_artery_length + 3.0  # dummy branch 3mm long by default

    # node 2 is 3 mm up from node 1 in the z direction
    node_loc[1][0] = 2
    node_loc[1][1] = insertion_x
    node_loc[1][2] = insertion_y
    node_loc[1][3] = z_radius * np.sqrt(
        1.0 - (node_loc[0][1] / x_radius ** 2) - (node_loc[0][2] / y_radius) ** 2) + umb_artery_length

    # node 3 & 4 is the start of the 'umbilical artery'
    node_loc[2][0] = 3
    node_loc[2][1] = insertion_x
    node_loc[2][2] = insertion_y - umb_artery_distance / 2.0
    node_loc[2][3] = z_radius * np.sqrt(
        1.0 - (node_loc[0][1] / x_radius) ** 2 - (node_loc[0][2] / y_radius) ** 2) + umb_artery_length
    node_loc[3][0] = 4
    node_loc[3][1] = insertion_x
    node_loc[3][2] = insertion_y + umb_artery_distance / 2.0
    node_loc[3][3] = z_radius * np.sqrt(
        1.0 - (node_loc[0][1] / x_radius) ** 2 - (node_loc[0][2] / y_radius) ** 2) + umb_artery_length

    # node 5 and 6 'hit' the chorionic plate.
    node_loc[4][0] = 5
    node_loc[4][1] = insertion_x
    node_loc[4][2] = insertion_y - umb_artery_distance / 2.0
    node_loc[4][3] = z_radius * np.sqrt(1.0 - (node_loc[4][1] / x_radius) ** 2 - (node_loc[4][2] / y_radius) ** 2)
    node_loc[5][0] = 6
    node_loc[5][1] = insertion_x
    node_loc[5][2] = insertion_y + umb_artery_distance / 2.0
    node_loc[5][3] = z_radius * np.sqrt(1.0 - (node_loc[5][1] / x_radius) ** 2 - (node_loc[5][2] / y_radius) ** 2)
    print(x_radius)
    elems[0, :] = [1, 1, 2]
    elems[1, :] = [2, 2, 3]
    elems[2, :] = [3, 2, 4]
    elems[3, :] = [4, 3, 5]
    elems[4, :] = [5, 4, 6]

    return {'umb_nodes': node_loc, 'umb_elems': elems}


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
