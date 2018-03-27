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

    # local arrays
    nstem = np.zeros((num_elems_new, 2))
    ld = np.zeros(len(datapoints))
    ld_np = np.zeros(len(datapoints))
    ldtmp1 = np.zeros(len(datapoints))

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
    # initialise ne_old (list of old terminals) to list of terminals in current geometry
    parentlist = group_elem_parent_term(0, initial_geom['elem_down'])
    ne_old = parentlist
    nt_bns = len(ne_old)  # Curerent number of terminals
    n_elm = nt_bns
    for n in range(0, len(parentlist)):
        nstem[n][0] = ne_old[n]

    # Initialise LD array to map each seed point to a parent branch.
    # For a single parent, all seed points will initially be mapped to
    # it; for multiple parents data_to_mesh is called to calculate
    # the closest parent end-point to each seed point.
    ld = ld + ne_old[0]
    ld_np = ld
    ld = data_to_mesh(ld, datapoints, ne_old, node_loc, elems)

    # Assign temporary arrays
    for nd in range(0, len(datapoints)):
        if ld[nd] != 0:
            ne_min = ld[nd]
            nstem[ne_min][1] = nstem[ne_min][1] + 1
            ldtmp1[nd] = ld[nd]

    # Dealing with cases when there are no datapoints assigned with an element
    nt_bns = 0
    for n in range(0, len(parentlist)):
        if nstem[parentlist[n]][1] != 0:
            ne_old[nt_bns] = parentlist[n]
            nt_bns = nt_bns + 1
    n_elm = nt_bns
    n_elm_temp = n_elm

    # remove terminal elements with only one data point associated and corresponding data pomt
    # from the group
    for n in range(0, n_elm):
        ne_min = ne_old[n]
        if nstem[ne_min][1] < 2:
            for nd in range(0, len(datapoints)):
                ldtmp1[nd] = ld[nd]
                ld[nd] = 0
                ne_old[n] = 0
                n_elm_temp = n_elm_temp - 1
    numtb = 0
    for n in range(0, n_elm):
        if ne_old[n] == 0:
            i = 0
            while (n + i < n_elm) and (ne_old[n + i] == 0):
                i = i + 1
            for m in range(n, n_elm - 1):
                ne_old[m] = ne_old[m + i]
            num_tb = num_tb + 1
    n_elm = n_elm_temp

    #START OF BIFURCATING DISTRIBUTATIVE ALGORITHM

    print('         gen       #brn       total#       #term      #data')


def calc_branch_direction(vector):
    length = 0.0
    for i in range(0, len(vector)):
        length = length + vector[i] ** 2
    length = np.sqrt(length)

    branch_direction = vector / length

    return branch_direction


def dist_two_vectors(vector1, vector2):
    dist = np.sqrt((vector1[0] - vector2[0]) ** 2 + (vector1[1] - vector2[1]) ** 2 + (vector1[2] - vector2[2]) ** 2)

    return dist


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


def data_to_mesh(ld, datapoints, parentlist, node_loc, elems):
    # Assigns data(seed) points to the closest ending of branches in the current generation.
    for nd in range(0, len(datapoints)):
        if ld[nd] != 0:
            min_dist = 1e10
            for noelem in range(0, len(parentlist)):
                ne = parentlist[noelem]
                np = elems[ne][2]
                dist = dist_two_vectors(node_loc[np][1:4], datapoints[nd])
                if dist < min_dist:
                    ne_min = ne
                    min_dist = dist
            ld[nd] = ne_min

    return ld

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
