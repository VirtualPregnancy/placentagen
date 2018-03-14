#!/usr/bin/env python


def equispaced_data_in_ellipsoid(n, volume, thickness, ellipticity):
    # Generates equally sapced data points in an ellipsoid with the following inputs
    # n=number of data points which we aim to generate
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
    estimated_spacing = (volume / n) ** (1. / 3)
    print(volume)
    print(n)
    print(estimated_spacing)
