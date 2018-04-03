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


def angle_two_vectors(vector1,vector2):
    vector1_u=vector1 / np.linalg.norm(vector1)
    vector2_u=vector2 / np.linalg.norm(vector2)

    if(np.equal(vector1_u,vector2_u)).all():
        print('true')
        angle = 0.0
    else:
        dotprod = np.dot(vector1_u,vector2_u)
        if np.isclose(1.0,dotprod):
            angle = 0
        else:
            angle = np.arccos(dotprod)


    return angle

