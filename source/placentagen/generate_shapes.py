#!/usr/bin/env python
import numpy as np


def equispaced_data_in_ellipsoid(n, volume, thickness, ellipticity):
    # Generates equally sapced data points in an ellipsoid with the following inputs
    # n=number of data points which we aim to generate
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
    pi = np.pi
    data_spacing = (volume / n) ** (1. / 3)
    z_radius = thickness / 2.0
    x_radius = np.sqrt(volume * 3.0 / (4.0 * pi * ellipticity * z_radius))
    y_radius = ellipticity * x_radius

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
