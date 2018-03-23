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


def umbilical_seed_geometry(volume, thickness, ellipticity, insertion_x, insertion_y, umb_artery_distance):
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
        1.0 - (node_loc[0][0] / x_radius) ** 2 - (node_loc[0][1] / y_radius) ** 2) + 23.0

    # node 2 is 3 mm up from node 1 in the z direction
    node_loc[1][0] = 2
    node_loc[1][1] = insertion_x
    node_loc[1][2] = insertion_y
    node_loc[1][3] = z_radius * np.sqrt(
        1.0 - (node_loc[0][0] / x_radius ** 2) - (node_loc[0][1] / y_radius) ** 2) + 20.0

    # node 3 & 4 is the start of the 'umbilical artery'
    node_loc[2][0] = 3
    node_loc[2][1] = insertion_x
    node_loc[2][2] = insertion_y - umb_artery_distance / 2.0
    node_loc[2][3] = z_radius * np.sqrt(
        1.0 - (node_loc[0][0] / x_radius) ** 2 - (node_loc[0][1] / y_radius) ** 2) + 20.0
    node_loc[3][0] = 4
    node_loc[3][1] = insertion_x
    node_loc[3][2] = insertion_y + umb_artery_distance / 2.0
    node_loc[3][3] = z_radius * np.sqrt(
        1.0 - (node_loc[0][0] / x_radius) ** 2 - (node_loc[0][1] / y_radius) ** 2) + 20.0

    # node 5 and 6 'hit' the chorionic plate.
    node_loc[4][0] = 5
    node_loc[4][1] = insertion_x
    node_loc[4][2] = insertion_y - umb_artery_distance / 2.0
    node_loc[4][3] = z_radius * np.sqrt(1.0 - (node_loc[4][0] / x_radius) ** 2 - (node_loc[4][1] / y_radius) ** 2)
    node_loc[5][0] = 6
    node_loc[5][1] = insertion_x
    node_loc[5][2] = insertion_y + umb_artery_distance / 2.0
    node_loc[5][3] = z_radius * np.sqrt(1.0 - (node_loc[5][0] / x_radius) ** 2 - (node_loc[5][1] / y_radius) ** 2)

    elems[0, :] = [1, 1, 2]
    elems[1, :] = [2, 2, 3]
    elems[2, :] = [3, 2, 4]
    elems[3, :] = [4, 3, 5]
    elems[4, :] = [5, 4, 6]

    return{'umb_nodes':node_loc, 'umb_elems':elems}


def calculate_ellipse_radii(volume, thickness, ellipticity):
    pi = np.pi
    z_radius = thickness / 2.0
    x_radius = np.sqrt(volume * 3.0 / (4.0 * pi * ellipticity * z_radius))
    y_radius = ellipticity * x_radius

    return {'x_radius': x_radius, 'y_radius': y_radius, 'z_radius': z_radius}
