#!/usr/bin/env python
import numpy as np


def grow_chorionic_surface(volume, thickness, ellipticity, datapoints, initial_geom):
    # We can estimate the number of elements in the generated model based on the number of data (seed points) to
    #  pre-allocate data arrays.
    est_generation = int(np.ceil(np.log(len(datapoints)) / np.log(2)))
    total_estimated = 0

    for i in range(0, est_generation + 1):
        total_estimated = total_estimated + 2 ** i
    num_elems_old = len(initial_geom["umb_elems"])
    num_nodes_old = len(initial_geom["umb_nodes"])
    num_elems_new = num_elems_old + total_estimated
    num_nodes_new = num_nodes_old + total_estimated
    # Pre-allocation of data arrays
    elem_directions = np.zeros((num_elems_new, 3))
    elem_order = np.zeros((num_elems_new, 3))
    node_loc = np.zeros((num_nodes_new, 4))
    node_loc[0:num_nodes_old][:] = initial_geom["umb_nodes"]
    elems = np.zeros((num_elems_new, 3))
    elems[0:num_elems_old][:] = initial_geom["umb_elems"]
    elem_upstream = np.zeros((num_elems_new, 3))
    elem_upstream[0:num_elems_old][:] = initial_geom['elem_up']
    elem_downstream = np.zeros((num_elems_new, 3))
    elem_downstream[0:num_elems_old][:] = initial_geom['elem_down']

    # Calculate the generations for the initial branches
    # Calculate the direction of the initial branches
    for ne in range(0, num_elems_old):
        if elem_upstream[ne][0] == 0.0:  # This is the stem (inlet) vessel
            elem_order[ne][0] = 1
        else:
            ne0 = elem_upstream[ne][1]
            elem_order[ne][0] = elem_order[ne0][0] + 1
        node_in = elems[ne][1]
        node_out = elems[ne][2]
        elem_directions[ne][:] = calc_branch_direction(node_loc[node_out][1:4] - node_loc[node_in][1:4])

    parentlist = group_elem_parent_term(0, initial_geom['elem_down'])


def calc_branch_direction(vector):
    length = 0.0
    for i in range(0, len(vector)):
        length = length + vector[i] ** 2
    length = np.sqrt(length)

    branch_direction = vector / length

    return branch_direction


def group_elem_parent_term(ne_parent, elem_downstream):
    parentlist = np.zeros(10)
    ne_old = np.zeros(10)
    ntemp_list = np.zeros(10)
    ne_temp = np.zeros(10)

    NT_BNS = 1
    ne_old[0] = ne_parent
    ne_count = 0
    ntemp_list[ne_count] = 0
    while NT_BNS != 0:
        num_nodes = NT_BNS
        NT_BNS = 0
        for m in range(1, num_nodes + 1):
            ne0 = ne_old[m]
            for n in range(0, int(elem_downstream[ne0][0])):
                NT_BNS = NT_BNS + 1
                ne_temp[NT_BNS] = elem_downstream[ne0][n + 1]
        for n in range(1, NT_BNS + 1):
            ne_old[n] = ne_temp[n]
            ne_count = ne_count + 1
            ntemp_list[ne_count] = ne_temp[n]
    ntemp_list[0] = ne_count
    num_in_list = 0
    for noelem in range(1, ne_count + 1):
        ne = ntemp_list[noelem]
        if elem_downstream[ne][0] == 0:
            num_in_list = num_in_list + 1
            parentlist[num_in_list - 1] = ne
    parentlist.resize(num_in_list)

    print('Grouped by parent,' + str(ne_parent) + ' No in parent list,' + str(num_in_list))

    return parentlist
