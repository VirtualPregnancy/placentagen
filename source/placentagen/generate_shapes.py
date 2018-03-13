#!/usr/bin/env python
import numpy as np
import math

def equispaced_data_in_ellipsoid():
    # Set problem parameters:
    volume = float(428)
    thickness = float(2.68)
    ellipse = float(1.68)
    n = input('no. of grid nodes (682000 for optimal mesh):')
    estimated_spacing=(volume/n)^(1/3)
    print(estimated_spacing)



