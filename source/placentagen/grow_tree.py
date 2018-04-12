#!/usr/bin/env python
import numpy as np

from . import pg_utilities


def grow_large_tree(angle_max, angle_min, fraction, min_length, point_limit,
                    volume, thickness, ellipticity, datapoints, initial_geom):
    # Calulate axis dimensions of ellipsoid with given volume, thickness and ellipticity
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    # We can estimate the number of elements in the generated model based on the number of data (seed points) to
    #  pre-allocate data arrays.
    est_generation = int(np.ceil(np.log(len(datapoints)) / np.log(2)))
    total_estimated = 0

    for i in range(0, est_generation + 1):
        total_estimated = total_estimated + 2 ** i
    # Define the total number of nodes and elements prior to growing, plus the new number expected
    num_elems_old = len(initial_geom["elems"])
    num_nodes_old = len(initial_geom["nodes"])
    num_elems_new = num_elems_old + total_estimated
    num_nodes_new = num_nodes_old + total_estimated
    original_data_length = len(datapoints)
    # Pre-allocation of data arrays
    # elem_directions = np.zeros((num_elems_new, 3))
    # elem_order = np.zeros((num_elems_new, 3))
    node_loc = np.zeros((num_nodes_new, 4))
    node_loc[0:num_nodes_old][:] = initial_geom["nodes"]
    elems = np.zeros((num_elems_new, 3), dtype=int)
    elems[0:num_elems_old][:] = initial_geom["elems"]
    elem_upstream = np.zeros((num_elems_new, 3), dtype=int)
    elem_upstream[0:num_elems_old][:] = initial_geom['elem_up']
    elem_downstream = np.zeros((num_elems_new, 3), dtype=int)
    elem_downstream[0:num_elems_old][:] = initial_geom['elem_down']

    # local arrays
    nstem = np.zeros((num_elems_new, 2), dtype=int)  # number of seeds per parent
    map_seed_to_elem = np.zeros(len(datapoints), dtype=int)  # seed to elem - initial array for groupings
    tb_list = np.zeros(2 * len(datapoints), dtype=int)  # number of terminal bronchioles

    # Set initial values for local and global nodes and elements
    ne = num_elems_old - 1  # current maximum element number
    nnod = num_nodes_old - 1  # current maximum node number
    numtb = 0  # count of terminals

    parentlist = group_elem_parent_term(0, initial_geom['elem_down'])  # master parent list

    # Initialise LD array to map each seed point to a parent branch.
    # For a single parent, all seed points will initially be mapped to
    # it; for multiple parents data_to_mesh is called to calculate
    # the closest parent end-point to each seed point.
    map_seed_to_elem = map_seed_to_elem + parentlist[0]
    map_seed_to_elem = data_to_mesh(map_seed_to_elem, datapoints, parentlist, node_loc, elems)

    for npar in range(0, len(parentlist)):
        print('Generating children for parent ' + str(npar) + '(elem #' + str(parentlist[npar]) + ') of a total of ' + str(len(parentlist)))
        current_parent = parentlist[npar]
        num_next_parents = 1
        data_current_parent = np.zeros((len(datapoints), 3))
        nd_for_parent = 0
        for nd in range(0, len(datapoints)):
            if map_seed_to_elem[nd] == current_parent:
                data_current_parent[nd_for_parent] = datapoints[nd][:]
                map_seed_to_elem[nd] = 0
                datapoints[nd][:] = 0
                nd_for_parent = nd_for_parent + 1

        datapoints = datapoints[np.nonzero(map_seed_to_elem)] #remove the datapoints you are currently analysing from master list
        map_seed_to_elem = map_seed_to_elem[np.nonzero(map_seed_to_elem)] #remove parent you are currently analysing from master list
        data_current_parent.resize(nd_for_parent, 3)
        local_parent_temp = np.zeros(num_elems_new, dtype=int)  # rezero local parent arrays
        local_parent = np.zeros(num_elems_new, dtype=int)
        local_parent[0] = current_parent
        map_seed_to_elem_new = np.zeros(len(data_current_parent), dtype=int)
        map_seed_to_elem_new = map_seed_to_elem_new + current_parent
        remaining_data = len(data_current_parent)
        original_data = len(data_current_parent)

        # START OF BIFURCATING DISTRIBUTATIVE ALGORITHM
        # could make this an optional output at a later date
        # print('         newgens       #brn       total#       #term      #data')

        ngen = 0  # for output, look to have each generation of elements recordded

        while num_next_parents != 0:
            # for ok in range(0,2):
            ngen = ngen + 1  # increment generation from parent for output
            num_parents = num_next_parents  # update the number of current parents
            num_next_parents = 0  # reset the number of local parents for next iteration
            noelem_gen = 0
            for m in range(0, num_parents):
                ne_parent = local_parent[m]
                com = mesh_com(ne_parent, map_seed_to_elem_new, data_current_parent)
                np_start = int(elems[ne_parent][2])
                np_prt_start = int(elems[ne_parent][1])
                # Split the seed points by the plane defined by the parent branch and the com of the
                # seedpoints attached to it
                split_data = data_splitby_plane(map_seed_to_elem_new, data_current_parent, com,
                                                node_loc[np_prt_start][1:4],
                                                node_loc[np_start][1:4], ne_parent, ne, point_limit)
                # Check that there ar enough seedpoints in both groups to proceed
                # Note long term one could allow one group to continue and the other not to
                if split_data['enough_points'][0] and split_data['enough_points'][1]:
                    for n in range(0, 2):
                        branch = True
                        # Calculate centre of mass to grow toward
                        com = mesh_com(ne + 1, map_seed_to_elem_new, data_current_parent)
                        # location of start of new element
                        start_node_loc = node_loc[np_start][1:4]
                        # length of new branch
                        length_new = np.linalg.norm(fraction * (com - start_node_loc))
                        # check that the new branch is long enough and if not need to adjust length and remove associated data point
                        if length_new < min_length:
                            length_new = min_length
                            branch = False
                        # calculate location of end node
                        end_node_loc = start_node_loc + length_new * (com - start_node_loc) / np.linalg.norm(
                            (com - start_node_loc))
                        # Checks that branch angles are appropriate
                        end_node_loc = mesh_check_angle(angle_min, angle_max, node_loc[elems[ne_parent][1]][1:4],
                                                        start_node_loc, end_node_loc, ne_parent, ne + 1)

                        # Create new elements and nodes
                        elems[ne + 1][0] = ne + 1  # creating new element
                        elems[ne + 1][1] = np_start  # starts at this node
                        elems[ne + 1][2] = nnod + 1  # ends at this node
                        elem_upstream[ne + 1][0] = 1  # The new element has one parent
                        elem_upstream[ne + 1][1] = ne_parent
                        elem_downstream[ne + 1][0] = 0  # the new element currently has no children

                        elem_downstream[ne_parent][0] = elem_downstream[ne_parent][
                                                            0] + 1  # the parent element gets this as a child
                        elem_downstream[ne_parent][elem_downstream[ne_parent][0]] = ne + 1

                        node_loc[nnod + 1][0] = nnod + 1
                        node_loc[nnod + 1][1:4] = end_node_loc
                        ne = ne + 1
                        nnod = nnod + 1

                        noelem_gen = noelem_gen + 1  # Only used for runtime output, number of new elements created in current generation

                        nstem[ne][0] = nstem[ne_parent][0]

                        if branch:
                            # We are making a new branch so this one becomes a parent for next time
                            local_parent_temp[num_next_parents] = ne
                            num_next_parents = num_next_parents + 1
                        else:
                            # Need to modify down the line not to branch in the next generation
                            local_parent_temp[num_next_parents] = ne
                            num_next_parents = num_next_parents + 1
                            # numtb = numtb + 1

                else:  # Not ss so no splitting
                    # Not enough seed points in the set during the split parent branch becomes a terminal
                    tb_list[numtb] = ne_parent
                    numtb = numtb + 1
                    min_dist = 1.0e10
                    count_data = 0
                    for nd in range(0, len(data_current_parent)):
                        if map_seed_to_elem_new[nd] != 0:
                            if map_seed_to_elem_new[nd] == ne + 1:
                                map_seed_to_elem_new[nd] = ne_parent  # give back to parent
                            if map_seed_to_elem_new[nd] == ne + 2:
                                map_seed_to_elem_new[nd] = ne_parent
                            dist = dist_two_vectors(data_current_parent[nd][:], node_loc[int(elems[ne_parent][2])][1:4])
                            if dist < min_dist:
                                if map_seed_to_elem_new[nd] == ne_parent:
                                    count_data = count_data + 1
                                    nd_min = nd
                                    min_dist = dist
                    if count_data != 0:  # If there were any data points associated
                        map_seed_to_elem_new[nd_min] = 0
                        remaining_data = remaining_data - 1

            # .......Copy the temporary list of branches to NE_OLD. These become the
            # .......parent elements for the next branching

            for n in range(0, num_next_parents):
                local_parent[n] = local_parent_temp[n]
                nstem[int(local_parent[n])][1] = 0  # initia count of data points

            if remaining_data < original_data:  # only need to reallocate data if we have lost some data points
                original_data = remaining_data
                # reallocate datapoints
                map_seed_to_elem_new = data_to_mesh(map_seed_to_elem_new, data_current_parent,
                                                    local_parent[0:num_next_parents], node_loc,
                                                    elems)
            # Could make this an optional output at a later date
            # print('   ' + str(ngen) + '   ' + str(noelem_gen) + '   ' + str(ne) +
            #      '   ' + str(numtb) + '   ' + str(remaining_data))

    elems.resize(ne + 1, 3, refcheck=False)
    elem_upstream.resize(ne + 1, 3, refcheck=False)
    elem_downstream.resize(ne + 1, 3, refcheck=False)
    node_loc.resize(nnod + 1, 4, refcheck=False)

    return {'nodes': node_loc, 'elems': elems, 'elem_up': elem_upstream, 'elem_down': elem_downstream}


def data_with_parent(current_parent, map_seed_to_elem, datapoints):
    new_datapoints = np.zeros((len(datapoints), 3))
    nd_for_parent = 0
    for nd in range(0, len(datapoints)):
        if map_seed_to_elem[nd] == current_parent:
            new_datapoints[nd_for_parent] = datapoints[nd][:]
            map_seed_to_elem[nd] = 0
            datapoints[nd][:] = 0
            nd_for_parent = nd_for_parent + 1
    datapoints = datapoints[np.nonzero(map_seed_to_elem)]
    map_seed_to_elem = map_seed_to_elem[np.nonzero(map_seed_to_elem)]
    new_datapoints.resize(nd_for_parent, 3)

    return new_datapoints


def grow_chorionic_surface(angle_max, angle_min, fraction, min_length, point_limit,
                           volume, thickness, ellipticity, datapoints, initial_geom, sorv):
    # Calulate axis dimensions of ellipsoid with given volume, thickness and ellipticity
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    # We can estimate the number of elements in the generated model based on the number of data (seed points) to
    #  pre-allocate data arrays.
    est_generation = int(np.ceil(np.log(len(datapoints)) / np.log(2)))
    total_estimated = 0

    for i in range(0, est_generation + 1):
        total_estimated = total_estimated + 2 ** i
    # Define the total number of nodes and elements prior to growing, plus the new number expected
    num_elems_old = len(initial_geom["elems"])
    num_nodes_old = len(initial_geom["nodes"])
    num_elems_new = num_elems_old + total_estimated
    num_nodes_new = num_nodes_old + total_estimated
    original_data = len(datapoints)
    # Pre-allocation of data arrays
    # elem_directions = np.zeros((num_elems_new, 3))
    # elem_order = np.zeros((num_elems_new, 3))
    node_loc = np.zeros((num_nodes_new, 4))
    node_loc[0:num_nodes_old][:] = initial_geom["nodes"]
    elems = np.zeros((num_elems_new, 3), dtype=int)
    elems[0:num_elems_old][:] = initial_geom["elems"]
    elem_upstream = np.zeros((num_elems_new, 3), dtype=int)
    elem_upstream[0:num_elems_old][:] = initial_geom['elem_up']
    elem_downstream = np.zeros((num_elems_new, 3), dtype=int)
    elem_downstream[0:num_elems_old][:] = initial_geom['elem_down']

    # local arrays
    nstem = np.zeros((num_elems_new, 2), dtype=int)
    local_parent_temp = np.zeros(num_elems_new, dtype=int)
    local_parent = np.zeros(num_elems_new, dtype=int)
    map_seed_to_elem = np.zeros(len(datapoints), dtype=int)
    tb_list = np.zeros(2 * len(datapoints), dtype=int)  # number of terminal bronchioles

    # initialise local_parent (list of old terminals) to list of terminals in current geometry
    parentlist = group_elem_parent_term(0, initial_geom['elem_down'])
    local_parent[0:len(parentlist)] = parentlist

    num_parents = len(parentlist)  # Curerent number of terminals
    num_next_parents = num_parents
    for n in range(0, len(parentlist)):
        nstem[n][0] = local_parent[n]

    # Initialise LD array to map each seed point to a parent branch.
    # For a single parent, all seed points will initially be mapped to
    # it; for multiple parents data_to_mesh is called to calculate
    # the closest parent end-point to each seed point.
    map_seed_to_elem = map_seed_to_elem + local_parent[0]
    map_seed_to_elem = data_to_mesh(map_seed_to_elem, datapoints, parentlist, node_loc, elems)
    # Assign temporary arrays
    for nd in range(0, len(datapoints)):
        if map_seed_to_elem[nd] != 0:
            ne_min = map_seed_to_elem[nd]  # element associated with this seed
            nstem[ne_min][1] = nstem[ne_min][1] + 1  # This element has x data pints

    # Dealing with cases when there are no datapoints assigned with an element
    num_parents = 0  # intialise number of parents
    for n in range(0, len(parentlist)):
        if nstem[parentlist[n]][1] != 0:  # if parent has a data point allocated then its  in local parent list
            local_parent[num_parents] = parentlist[n]
            num_parents = num_parents + 1
    num_next_parents = num_parents
    n_elm_temp = num_next_parents
    # remove terminal elements with only one data point associated and corresponding data pomt
    # from the group
    numtb = 0
    for n in range(0, num_next_parents):
        ne_min = local_parent[n]  # this is the parent in the local list (has at least one datapoint)
        if nstem[ne_min][1] < 2:  # If only one datapint
            for nd in range(0, len(datapoints)):
                if (map_seed_to_elem[nd] == ne_min):
                    map_seed_to_elem[nd] = 0  # remove datapoint from the list
                    local_parent[n] = 0  # Not a parent, wont branch
                    n_elm_temp = n_elm_temp - 1  #
                    tb_list[numtb] = ne_min
                    numtb = numtb + 1

    for n in range(0, num_next_parents):
        if local_parent[n] == 0:
            i = 0
            while (n + i < num_next_parents) and (local_parent[n + i] == 0):
                i = i + 1
            for m in range(n, num_next_parents - 1):
                local_parent[m] = local_parent[m + i]
    num_next_parents = n_elm_temp

    remaining_data = 0
    for nd in range(0, len(datapoints)):
        if map_seed_to_elem[nd] > 0:
            remaining_data = remaining_data + 1
    # START OF BIFURCATING DISTRIBUTATIVE ALGORITHM
    # Could potentially make this an optional output at a later date
    # print('         newgens       #brn       total#       #term      #data')

    # Set initial values for local and global nodes and elements
    ne = num_elems_old - 1  # current maximum element number
    nnod = num_nodes_old - 1  # current maximum node number

    ngen = 0  # for output, look to have each generation of elements recordded

    while num_next_parents != 0:
        # for ok in range(0,2):
        ngen = ngen + 1  # increment generation from parent for output
        num_parents = num_next_parents  # update the number of current parents
        num_next_parents = 0  # reset the number of local parents for next iteration
        noelem_gen = 0
        for m in range(0, num_parents):
            ne_parent = local_parent[m]
            com = mesh_com(ne_parent, map_seed_to_elem, datapoints)
            np_start = int(elems[ne_parent][2])
            np_prt_start = int(elems[ne_parent][1])
            # Split the seed points by the plane defined by the parent branch and the com of the
            # seedpoints attached to it
            if sorv is 'surface':
                split_data = data_splitby_xy(map_seed_to_elem, datapoints, com, node_loc[np_start][1:3], ne_parent, ne,
                                             point_limit)
            else:
                split_data = data_splitby_plane(map_seed_to_elem, datapoints, com, node_loc[np_prt_start][1:4],
                                                node_loc[np_start][1:4], ne_parent, ne, point_limit)
            # Check that there ar enough seedpoints in both groups to proceed
            # Note long term one could allow one group to continue and the other not to
            if split_data['enough_points'][0] and split_data['enough_points'][1]:
                for n in range(0, 2):
                    branch = True
                    # Calculate centre of mass to grow toward
                    com = mesh_com(ne + 1, map_seed_to_elem, datapoints)
                    # location of start of new element
                    start_node_loc = node_loc[np_start][1:4]
                    # length of new branch
                    length_new = np.linalg.norm(fraction * (com - start_node_loc))
                    # check that the new branch is long enough and if not need to adjust length and remove associated data point
                    if length_new < min_length:
                        length_new = min_length
                        branch = False
                    # calculate location of end node
                    end_node_loc = start_node_loc + length_new * (com - start_node_loc) / np.linalg.norm(
                        (com - start_node_loc))
                    # Checks that branch angles are appropriate
                    if sorv is 'surface':
                        # first node is 1st parent node
                        node1 = np.array(
                            [node_loc[elems[ne_parent][1]][0], node_loc[elems[ne_parent][1]][1], 0])
                        # second parent node (start of new branch
                        node2 = np.array([start_node_loc[0], start_node_loc[1], 0])
                        # end of new branch
                        node3 = np.array([end_node_loc[0], end_node_loc[1], 0])
                        end_node = mesh_check_angle(angle_min, angle_max, node1, node2, node3, ne_parent, ne + 1)
                        end_node_loc[0:2] = end_node[0:2]
                        end_node_loc[2] = pg_utilities.z_from_xy(end_node[0], end_node[1], x_radius, y_radius, z_radius)
                    elif sorv is 'volume':
                        end_node_loc = mesh_check_angle(angle_min, angle_max, node_loc[elems[ne_parent][1]][1:4],
                                                        start_node_loc, end_node_loc, ne_parent, ne + 1)

                    # Create new elements and nodes
                    elems[ne + 1][0] = ne + 1  # creating new element
                    elems[ne + 1][1] = np_start  # starts at this node
                    elems[ne + 1][2] = nnod + 1  # ends at this node
                    elem_upstream[ne + 1][0] = 1  # The new element has one parent
                    elem_upstream[ne + 1][1] = ne_parent
                    elem_downstream[ne + 1][0] = 0  # the new element currently has no children

                    elem_downstream[ne_parent][0] = elem_downstream[ne_parent][
                                                        0] + 1  # the parent element gets this as a child
                    elem_downstream[ne_parent][elem_downstream[ne_parent][0]] = ne + 1

                    node_loc[nnod + 1][0] = nnod + 1
                    node_loc[nnod + 1][1:4] = end_node_loc
                    ne = ne + 1
                    nnod = nnod + 1

                    noelem_gen = noelem_gen + 1  # Only used for runtime output, number of new elements created in current generation

                    nstem[ne][0] = nstem[ne_parent][0]

                    if branch:
                        # We are making a new branch so this one becomes a parent for next time
                        local_parent_temp[num_next_parents] = ne
                        num_next_parents = num_next_parents + 1
                    else:
                        # Note that this needs to be updated so that the branch and its datapoints get deleted from the list and dont grow nest time
                        # but for python speed reasons we should leave here for now
                        local_parent_temp[num_next_parents] = ne
                        num_next_parents = num_next_parents + 1
                        # numtb = numtb + 1

            else:  # Not ss so no splitting
                # Not enough seed points in the set during the split parent branch becomes a terminal
                tb_list[numtb] = ne_parent
                numtb = numtb + 1
                min_dist = 1.0e10
                count_data = 0
                for nd in range(0, len(datapoints)):
                    if map_seed_to_elem[nd] != 0:
                        if map_seed_to_elem[nd] == ne + 1:
                            map_seed_to_elem[nd] = ne_parent  # give back to parent
                        if map_seed_to_elem[nd] == ne + 2:
                            map_seed_to_elem[nd] = ne_parent
                        dist = dist_two_vectors(datapoints[nd][:], node_loc[int(elems[ne_parent][2])][1:4])
                        if dist < min_dist:
                            if map_seed_to_elem[nd] == ne_parent:
                                count_data = count_data + 1
                                nd_min = nd
                                min_dist = dist
                if count_data != 0:  # If there were any data points associated
                    map_seed_to_elem[nd_min] = 0
                    remaining_data = remaining_data - 1

        # .......Copy the temporary list of branches to NE_OLD. These become the
        # .......parent elements for the next branching

        for n in range(0, num_next_parents):
            local_parent[n] = local_parent_temp[n]
            nstem[int(local_parent[n])][1] = 0  # initialaw count of data points

        if remaining_data < original_data:  # only need to reallocate data if we have lost some data points
            original_data = remaining_data
            # reallocate datapoints
            map_seed_to_elem = data_to_mesh(map_seed_to_elem, datapoints, local_parent[0:num_next_parents], node_loc,
                                            elems)
        # Could potentially make this an optional output at a later date
        # print('   ' + str(ngen) + '   ' + str(noelem_gen) + '   ' + str(ne) +
        #      '   ' + str(numtb) + '   ' + str(remaining_data))

    elems.resize(ne + 1, 3, refcheck=False)
    elem_upstream.resize(ne + 1, 3, refcheck=False)
    elem_downstream.resize(ne + 1, 3, refcheck=False)
    node_loc.resize(nnod + 1, 4, refcheck=False)

    return {'nodes': node_loc, 'elems': elems, 'elem_up': elem_upstream, 'elem_down': elem_downstream}


def refine_1D(initial_geom, from_elem):
    # Estimate new number of nodes and elements
    num_elems_old = len(initial_geom['elems'])
    num_nodes_old = len(initial_geom['nodes'])
    num_elems_new = 2 * num_elems_old - from_elem + 1
    num_nodes_new = num_nodes_old + num_elems_old - from_elem + 1

    elems_old = initial_geom['elems']
    node_loc_old = initial_geom['nodes']
    elem_upstream_old = initial_geom['elem_up']
    elem_downstream_old = initial_geom['elem_down']

    # allocate temporary arrays
    elems = np.zeros((num_elems_new, 3), dtype=int)
    node_loc = np.zeros((num_nodes_new, 4))
    elem_upstream = np.zeros((num_elems_new, 3), dtype=int)
    elem_downstream = np.zeros((num_elems_new, 3), dtype=int)
    old_2_new = np.zeros((num_elems_old, 2), dtype=int)

    for nnod in range(0, num_nodes_new):
        node_loc[nnod][0] = nnod  # will be a numeric list of nodes

    # copy non-split elements and nodes into the arrays
    if from_elem != 0:  # nothing needs to be done if from elem is zero
        node_loc[0][:] = initial_geom['nodes'][0][:]
        elems[0:from_elem][:] = initial_geom['elems'][0:from_elem][:]
        elem_upstream[0:from_elem][:] = initial_geom['elem_up'][0:from_elem][:]
        elem_downstream[0:from_elem][:] = initial_geom['elem_down'][0:from_elem][:]
        for ne in range(0, from_elem):
            # grab second node from each non-split element
            node_loc[elems_old[ne][2]][:] = node_loc_old[elems_old[ne][2]][:]
            old_2_new[ne][0] = ne
            old_2_new[ne][1] = ne

    new_elems = from_elem - 1
    nnod = from_elem
    for ne in range(from_elem, num_elems_old):
        # split element in two at half way point
        old_2_new[ne][0] = new_elems + 1
        old_2_new[ne][1] = new_elems + 2
        in_node = elems[old_2_new[elem_upstream_old[ne][1]][1]][2]
        in_point = node_loc_old[elems_old[ne][1]][1:4]
        out_point = node_loc_old[elems_old[ne][2]][1:4]
        mid_point = np.zeros(3)
        for i in range(0, 3):
            mid_point[i] = (in_point[i] + out_point[i]) / 2.0
        # Create a new node at midpoint
        nnod = nnod + 1
        node_loc[nnod][1:4] = mid_point
        # Next node is outpoint
        nnod1 = nnod + 1
        node_loc[nnod1][1:4] = out_point
        # create first new element
        new_elems = new_elems + 1
        elems[new_elems][0] = new_elems
        elems[new_elems][1] = in_node
        elems[new_elems][2] = nnod

        new_elems = new_elems + 1
        elems[new_elems][0] = new_elems
        elems[new_elems][1] = nnod
        elems[new_elems][2] = nnod1

        nnod = nnod1

    node_loc.resize(nnod + 1, 4, refcheck=False)
    elems.resize(new_elems + 1, 3, refcheck=False)

    elem_connectivity = pg_utilities.element_connectivity_1D(node_loc, elems)

    return {'nodes': node_loc, 'elems': elems, 'elem_up': elem_connectivity['elem_up'],
            'elem_down': elem_connectivity['elem_down']}


def add_stem_villi(initial_geom, from_elem, sv_length):
    # Estimate new number of nodes and elements
    num_elems_old = len(initial_geom['elems'])
    num_nodes_old = len(initial_geom['nodes'])

    elems_old = initial_geom['elems']
    node_loc_old = initial_geom['nodes']
    elem_upstream_old = initial_geom['elem_up']
    elem_downstream_old = initial_geom['elem_down']

    # calculate number pf new elements to add
    num_2_add = 0
    for ne in range(from_elem, num_elems_old):
        if elem_downstream_old[ne][0] == 1:
            num_2_add = num_2_add + 1
    # Allocate arrays
    num_elems_new = num_elems_old + num_2_add
    num_nodes_new = num_nodes_old + num_2_add
    elems = np.zeros((num_elems_new, 3), dtype=int)
    node_loc = np.zeros((num_nodes_new, 4))

    node_loc[0:num_nodes_old][:] = initial_geom['nodes']
    elems[0:num_elems_old][:] = initial_geom['elems']

    # Create new elements
    nnod = num_nodes_old - 1
    ne0 = num_elems_old - 1
    for ne in range(from_elem, num_elems_old):
        if elem_downstream_old[ne][0] == 1:
            nnod = nnod + 1
            ne0 = ne0 + 1
            node_loc[nnod][0] = nnod
            elems[ne0][0] = ne0
            elems[ne0][1] = elems[ne][2]
            elems[ne0][2] = nnod

            node_loc[nnod][1:3] = node_loc[elems[ne][2]][1:3]
            node_loc[nnod][3] = node_loc[elems[ne][2]][3] - sv_length

    elem_connectivity = pg_utilities.element_connectivity_1D(node_loc, elems)

    return {'nodes': node_loc, 'elems': elems, 'elem_up': elem_connectivity['elem_up'],
            'elem_down': elem_connectivity['elem_down']}


def mesh_check_angle(angle_min, angle_max, node1, node2, node3, ne_parent, myno):
    normal_to_plane = np.zeros(3)
    vector_cross_n = np.zeros(3)
    vector1 = (node2 - node1)
    vector1_u = vector1 / np.linalg.norm(vector1)
    vector2 = (node3 - node2)
    vector2_u = vector2 / np.linalg.norm(vector2)
    # angle = pg_utilities.angle_two_vectors(vector1, vector2)
    # if angle == 0:
    #    print('in',vector1 / np.linalg.norm(vector1),vector2 / np.linalg.norm(vector2))

    if (np.isclose(vector1_u, vector2_u)).all() or (np.isclose(vector1_u, -1.0 * vector2_u)).all():
        #   perturb new node slightly to 'misalign vectors and allow for normal to plane to be calculated
        length = np.linalg.norm(vector2)
        node3[0] = node3[0] + np.random.uniform(-0.01 * length, 0.01 * length)
        node3[1] = node3[1] + np.random.uniform(-0.01 * length, 0.01 * length)
        node3[2] = node3[2] + np.random.uniform(-0.01 * length, 0.01 * length)
        vector2 = node3 - node2

    # want to rotate vector 2 wrt vector 1
    angle = pg_utilities.angle_two_vectors(vector1, vector2)
    if angle == 0:
        print('out', vector1 / np.linalg.norm(vector1), vector2 / np.linalg.norm(vector2))

    normal_to_plane[0] = (vector2[1] * vector1[2] - vector2[2] * vector1[1])
    normal_to_plane[1] = (vector2[0] * vector1[2] - vector2[2] * vector1[0])
    normal_to_plane[2] = (vector2[0] * vector1[1] - vector2[1] * vector1[0])

    normal_to_plane_u = normal_to_plane / np.linalg.norm(normal_to_plane)

    vector_cross_n[0] = (vector1[1] * normal_to_plane[2] - vector1[2] * normal_to_plane[1])
    vector_cross_n[1] = (vector1[0] * normal_to_plane[2] - vector1[2] * normal_to_plane[0])
    vector_cross_n[2] = (vector1[0] * normal_to_plane[1] - vector1[1] * normal_to_plane[0])

    dotprod = np.dot(vector2, vector_cross_n)

    if dotprod < 0:
        angle = -1 * angle

    if abs(angle) < angle_min:
        # Need to adjust node3 to get a angle equal to  angle_min
        angle0 = abs(angle)
        angle_rot = (angle0 - angle_min)  # amount of angle to add

        # need to rotate around axis normal to the plane that the original and the vector make
        # unit vector normal to plane:
        R = np.zeros((3, 3))
        R[0][0] = np.cos(angle_rot) + normal_to_plane_u[0] ** 2 * (1 - np.cos(angle_rot))
        R[0][1] = normal_to_plane_u[0] * normal_to_plane_u[1] * (1 - np.cos(angle_rot)) - normal_to_plane_u[2] * np.sin(
            angle_rot)
        R[0][2] = normal_to_plane_u[0] * normal_to_plane[2] * (1 - np.cos(angle_rot)) + normal_to_plane_u[1] * np.sin(
            angle_rot)
        R[1][0] = normal_to_plane_u[0] * normal_to_plane_u[1] * (1 - np.cos(angle_rot)) + normal_to_plane_u[2] * np.sin(
            angle_rot)
        R[1][1] = np.cos(angle_rot) + normal_to_plane_u[1] ** 2 * (1 - np.cos(angle_rot))
        R[1][2] = normal_to_plane_u[1] * normal_to_plane_u[2] * (1 - np.cos(angle_rot)) - normal_to_plane_u[0] * np.sin(
            angle_rot)
        R[2][0] = normal_to_plane_u[0] * normal_to_plane[2] * (1 - np.cos(angle_rot)) - normal_to_plane_u[1] * np.sin(
            angle_rot)
        R[2][1] = normal_to_plane_u[1] * normal_to_plane_u[2] * (1 - np.cos(angle_rot)) + normal_to_plane_u[0] * np.sin(
            angle_rot)
        R[2][2] = np.cos(angle_rot) + normal_to_plane_u[2] ** 2 * (1 - np.cos(angle_rot))
        nu_vec = np.zeros(3)
        nu_vec[0] = R[0][0] * vector2[0] + R[0][1] * vector2[1] + R[0][2] * vector2[2]
        nu_vec[1] = R[1][0] * vector2[0] + R[1][1] * vector2[1] + R[1][2] * vector2[2]
        nu_vec[2] = R[2][0] * vector2[0] + R[2][1] * vector2[1] + R[2][2] * vector2[2]

        node3 = node2 + nu_vec

    elif abs(angle) > angle_max:
        # Need to adjust node3 to get a angle equal to  angle_max
        angle0 = abs(angle)
        angle_rot = (angle0 - angle_max)  # amount of angle to add

        # need to rotate around axis normal to the plane that the original and the vector make
        # unit vector normal to plane:
        R = np.zeros((3, 3))
        R[0][0] = np.cos(angle_rot) + normal_to_plane_u[0] ** 2 * (1 - np.cos(angle_rot))
        R[0][1] = normal_to_plane_u[0] * normal_to_plane_u[1] * (1 - np.cos(angle_rot)) - normal_to_plane_u[2] * np.sin(
            angle_rot)
        R[0][2] = normal_to_plane_u[0] * normal_to_plane[2] * (1 - np.cos(angle_rot)) + normal_to_plane_u[1] * np.sin(
            angle_rot)
        R[1][0] = normal_to_plane_u[0] * normal_to_plane_u[1] * (1 - np.cos(angle_rot)) + normal_to_plane_u[2] * np.sin(
            angle_rot)
        R[1][1] = np.cos(angle_rot) + normal_to_plane_u[1] ** 2 * (1 - np.cos(angle_rot))
        R[1][2] = normal_to_plane_u[1] * normal_to_plane_u[2] * (1 - np.cos(angle_rot)) - normal_to_plane_u[0] * np.sin(
            angle_rot)
        R[2][0] = normal_to_plane_u[0] * normal_to_plane[2] * (1 - np.cos(angle_rot)) - normal_to_plane_u[1] * np.sin(
            angle_rot)
        R[2][1] = normal_to_plane_u[1] * normal_to_plane_u[2] * (1 - np.cos(angle_rot)) + normal_to_plane_u[0] * np.sin(
            angle_rot)
        R[2][2] = np.cos(angle_rot) + normal_to_plane_u[2] ** 2 * (1 - np.cos(angle_rot))

        nu_vec = np.zeros(3)
        nu_vec[0] = R[0][0] * vector2[0] + R[0][1] * vector2[1] + R[0][2] * vector2[2]
        nu_vec[1] = R[1][0] * vector2[0] + R[1][1] * vector2[1] + R[1][2] * vector2[2]
        nu_vec[2] = R[2][0] * vector2[0] + R[2][1] * vector2[1] + R[2][2] * vector2[2]

        node3 = node2 + nu_vec

    return node3


def data_splitby_xy(ld, datapoints, x0, x1, ne_parent, ne_current, point_limit):
    # Decides which side of a line a random point is on.
    # !Given a directed line from point p0(x0, y0) to p1(x1, y1), you can use the following condition to decide whether a point p2(x2, y2) is on the left of the line, on the right, or on the same line:
    # value = (x1 - x0)(y2 - y0) - (x2 - x0)(y1 - y0)
    # if value > 0, p2 is on the left side of the line.
    #
    npoints = 0
    dat1 = 0
    dat2 = 0
    ss = [True, True]
    for nd in range(0, len(datapoints)):
        nsp = ld[nd]

        if nsp == ne_parent:  # data point belongs to this element
            npoints = npoints + 1
            checkvalue = (x1[0] - x0[0]) * (datapoints[nd][1] - x0[1]) - (datapoints[nd][0] - x0[0]) * (x1[1] - x0[1])
            if checkvalue >= 0:
                dat1 = dat1 + 1
                ld[nd] = ne_current + 1
            else:
                dat2 = dat2 + 1
                ld[nd] = ne_current + 2
    if npoints < point_limit:
        ss[0] = False
        ss[1] = False
    else:
        ss[0] = True
        ss[1] = True
        if dat1 < point_limit:
            ss[0] = False
        if dat2 < point_limit:
            ss[0] = False

    return {'enough_points': ss}


def data_splitby_plane(ld, datapoints, x0, x1, x2, ne_parent, ne_current, point_limit):
    # Decides which side of a plane a random point is on.

    # x0 is location of com
    # x1 is start node of parent
    # x2 is end node of parent

    npoints = 0
    dat1 = 0
    dat2 = 0
    ss = [True, True]

    colinear = pg_utilities.check_colinear(x0, x1, x2)
    if colinear:
        print('WARNING: Colinear is true')
    plane = pg_utilities.plane_from_3_pts(x0, x1, x2, False)
    for nd in range(0, len(datapoints)):
        nsp = ld[nd]
        if nsp == ne_parent:  # data point belongs to this element
            npoints = npoints + 1
            checkvalue = 0.0
            for i in range(0, 3):
                checkvalue = checkvalue + plane[i] * datapoints[nd][i]
            checkvalue = -1.0 * checkvalue - plane[3]
            if checkvalue >= 0:
                dat1 = dat1 + 1
                ld[nd] = ne_current + 1
            else:
                dat2 = dat2 + 1
                ld[nd] = ne_current + 2
    if npoints < point_limit:
        ss[0] = False
        ss[1] = False
    else:
        ss[0] = True
        ss[1] = True
        if dat1 < point_limit:
            ss[0] = False
        if dat2 < point_limit:
            ss[0] = False

    return {'enough_points': ss}


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
    parentlist = np.zeros(1000, dtype=int)
    ne_old = np.zeros(1000, dtype=int)
    ntemp_list = np.zeros(1000, dtype=int)
    ne_temp = np.zeros(1000, dtype=int)

    NT_BNS = 1
    ne_old[0] = ne_parent
    ne_count = 0
    ntemp_list[ne_count] = 0
    while NT_BNS != 0:
        num_nodes = NT_BNS
        NT_BNS = 0
        for m in range(1, num_nodes + 1):
            ne0 = int(ne_old[m])
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
        ne = int(ntemp_list[noelem])
        if elem_downstream[ne][0] == 0:
            num_in_list = num_in_list + 1
            parentlist[num_in_list - 1] = ne
    parentlist.resize(num_in_list)
    ne_old.resize(num_in_list)
    ntemp_list.resize(num_in_list)
    ne_temp.resize(num_in_list)
    print('Grouped by parent,' + str(ne_parent) + ' No in parent list,' + str(num_in_list))

    return parentlist


def data_to_mesh(ld, datapoints, parentlist, node_loc, elems):
    # Assigns data(seed) points to the closest ending of branches in the current generation.
    for nd in range(0, len(datapoints)):
        if ld[nd] != 0:
            ne_min = 0
            min_dist = 1e10
            for noelem in range(0, len(parentlist)):
                ne = int(parentlist[noelem])
                np = int(elems[ne][2])
                dist = dist_two_vectors(node_loc[np][1:4], datapoints[nd])
                if dist < min_dist:
                    ne_min = ne
                    min_dist = dist
            ld[nd] = ne_min

    return ld


def umbilical_seed_geometry(volume, thickness, ellipticity, insertion_x, insertion_y, umb_artery_distance,
                            umb_artery_length, datapoints):
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
    elem_downstream[3][0] = 2
    elem_downstream[3][1] = 5
    elem_downstream[3][2] = 6

    elems[4, :] = [4, 3, 5]
    elem_upstream[4][0] = 1
    elem_upstream[4][1] = 2
    elem_downstream[4][0] = 2
    elem_downstream[4][1] = 7
    elem_downstream[4][2] = 8

    # split datapoints by x and y insertion points
    ld = np.zeros(len(datapoints))
    for nd in range(0, len(datapoints)):
        if (datapoints[nd][0] > insertion_x) and (datapoints[nd][1] > insertion_y):
            ld[nd] = 1
        elif (datapoints[nd][0] > insertion_x) and (datapoints[nd][1] < insertion_y):
            ld[nd] = 2
        elif (datapoints[nd][0] < insertion_x) and (datapoints[nd][1] < insertion_y):
            ld[nd] = 3
        else:
            ld[nd] = 4

    # Calculate centre of mass for each seedpoint set
    com1 = mesh_com(1, ld, datapoints)
    com2 = mesh_com(2, ld, datapoints)
    com3 = mesh_com(3, ld, datapoints)
    com4 = mesh_com(4, ld, datapoints)

    ##CREATE NEW NODES ON CHORION SURFACE HALF WAY TO GROUP COM

    node_loc = np.append(node_loc, np.zeros((4, 4)), axis=0)
    elems = np.append(elems, np.zeros((4, 3)), axis=0)
    elem_upstream = np.append(elem_upstream, np.zeros((4, 3)), axis=0)
    elem_downstream = np.append(elem_downstream, np.zeros((4, 3)), axis=0)
    node_loc[6][0] = 6
    node_loc[6][1:4] = com1 / 2.0
    node_loc[6][3] = pg_utilities.z_from_xy(node_loc[6][1], node_loc[6][2], x_radius, y_radius, z_radius)

    elems[7, :] = [7, 5, 6]
    elem_upstream[7][0] = 1
    elem_upstream[7][1] = 4
    elem_downstream[7][0] = 0

    node_loc[7][0] = 7
    node_loc[7][1:4] = com2 / 2.0
    node_loc[7][3] = pg_utilities.z_from_xy(node_loc[7][1], node_loc[7][2], x_radius, y_radius, z_radius)

    elems[5, :] = [5, 4, 7]
    elem_upstream[5][0] = 1
    elem_upstream[5][1] = 3
    elem_downstream[5][0] = 0

    node_loc[8][0] = 8
    node_loc[8][1:4] = com3 / 2.0
    node_loc[8][3] = pg_utilities.z_from_xy(node_loc[8][1], node_loc[8][2], x_radius, y_radius, z_radius)

    elems[6, :] = [6, 4, 8]
    elem_upstream[6][0] = 1
    elem_upstream[6][1] = 3
    elem_downstream[6][0] = 0

    node_loc[9][0] = 9
    node_loc[9][1:4] = com4 / 2.0
    node_loc[9][3] = pg_utilities.z_from_xy(node_loc[9][1], node_loc[9][2], x_radius, y_radius, z_radius)

    elems[8, :] = [8, 5, 9]
    elem_upstream[8][0] = 1
    elem_upstream[8][1] = 4
    elem_downstream[8][0] = 0

    return {'nodes': node_loc, 'elems': elems, 'elem_up': elem_upstream, 'elem_down': elem_downstream}


def mesh_com(ld_val, ld, datapoints):
    dat = 0

    com = np.zeros(3)
    for nd in range(0, len(datapoints)):
        nsp = ld[nd]
        if nsp == ld_val:
            dat = dat + 1
            for nj in range(0, 3):
                com[nj] = com[nj] + datapoints[nd][nj]
    if dat != 0:
        com = com / dat

    return com
