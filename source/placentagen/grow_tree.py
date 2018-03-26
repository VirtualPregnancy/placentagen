#!/usr/bin/env python
import numpy as np

def grow_chorionic_surface(volume, thickness, ellipticity, datapoints, initial_geom):
    # We can estimate the number of elements in the generated model based on the number of data (seed points) to
    #  pre-allocate data arrays.
    est_generation = int(np.ceil(np.log(len(datapoints))/np.log(2)))
    total_estimated = 0

    for i in range(0,est_generation+1):
        total_estimated = total_estimated + 2**i

    num_elems_new = len(initial_geom["umb_elems"]) + total_estimated
    num_nodes_new = len(initial_geom["umb_elems"]) + total_estimated
    # Pre-allocation of data arrays
    elem_directions = np.zeros((num_elems_new,3))
    elem_order = np.zeros((num_elems_new,3))
    node_loc = np.zeros((num_nodes_new , 4))
    elems = np.zeros((num_elems_new , 3))
    elem_upstream = np.zeros((num_elems_new, 3))
    elem_downstream = np.zeros((num_elems_new, 3))



