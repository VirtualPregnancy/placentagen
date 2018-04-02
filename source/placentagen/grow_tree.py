#!/usr/bin/env python
import numpy as np

from . import pg_utilities


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
    ne_temp = np.zeros(num_elems_new)
    ne_old = np.zeros(1000)
    ld = np.zeros(len(datapoints))
    ld_np = np.zeros(len(datapoints))
    ldtmp1 = np.zeros(len(datapoints))

    numtbsum = 0

    # Calculate the generations for the initial branches
    # Calculate the direction of the initial branches

    for ne in range(0, num_elems_old):
        if elem_upstream[ne][0] == 0.0:  # This is the stem (inlet) vessel
            elem_order[ne][0] = 1
        else:
            ne0 = int(elem_upstream[ne][1])
            elem_order[ne][0] = elem_order[ne0][0] + 1
        node_in = int(elems[ne][1])
        node_out = int(elems[ne][2])
        elem_directions[ne][:] = calc_branch_direction(node_loc[node_out][1:4] - node_loc[node_in][1:4])
    # initialise ne_old (list of old terminals) to list of terminals in current geometry
    parentlist = group_elem_parent_term(0, initial_geom['elem_down'])
    ne_old[0:len(parentlist)] = parentlist
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
            ne_min = int(ld[nd])
            nstem[ne_min][1] = nstem[ne_min][1] + 1
            ldtmp1[nd] = ld[nd]

    # Dealing with cases when there are no datapoints assigned with an element
    nt_bns = 0
    for n in range(0, len(parentlist)):
        if nstem[int(parentlist[n])][1] != 0:
            ne_old[nt_bns] = parentlist[n]
            nt_bns = nt_bns + 1
    n_elm = nt_bns
    n_elm_temp = n_elm

    # remove terminal elements with only one data point associated and corresponding data pomt
    # from the group
    for n in range(0, n_elm):
        ne_min = int(ne_old[n])
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
            numtb = numtb + 1
    n_elm = n_elm_temp

    ld_num = 0
    for nd in range(0, len(datapoints)):
        if ld[nd] > 0:
            ld_num = ld_num + 1
    # START OF BIFURCATING DISTRIBUTATIVE ALGORITHM

    print('         gen       #brn       total#       #term      #data')

    # Set initial values for local and global nodes and elements
    ne = num_elems_old - 1  # current maximum element number
    nnod = num_nodes_old - 1  # current maximum node number
    ne_start = ne

    ne_parent = int(parentlist[0])
    ngen = 3
    ne_min = ne

    while n_elm != 0:
    #for ok in range(0,2):
        ngen = ngen + 1  # increment generation from parent
        nt_bns = n_elm
        n_elm = 0
        numzero = 0
        noelem_gen = 0
        max_dist = 0.0

        for m in range(0, nt_bns):
            ne_parent = int(ne_old[m])
            com = mesh_com(ne_parent, ld, datapoints)
            ne_grnd_parent = elem_upstream[ne_parent][1]  # Assumes only one parent, true in diverging tree
            np_start = int(elems[ne_parent][2])
            np_prt_start = int(elems[ne_parent][1])
            np_grnd_start = int(elems[ne_grnd_parent][1])
            if elem_downstream[ne_grnd_parent][1] == ne_parent:
                if elem_downstream[ne_grnd_parent][2] == 0:
                    np4 = int(elems[elem_upstream[ne_grnd_parent][1]][1])
                else:
                    np4 = int(elems[elem_downstream[ne_grnd_parent][2]][2])
            else:
                np4 = int(elems[elem_downstream[ne_grnd_parent][1]][2])

            length_parent = dist_two_vectors(node_loc[int(elems[ne_parent][1])][1:4], node_loc[int(elems[ne_parent][2])][1:4])


            # Split the seed points by the plane defined by the parent branch and the com of the
            # seedpoints attached to it
            if sorv is 'surface':
                split_data = data_splitby_xy(ld, datapoints, com, node_loc[np_start][1:3], ne_parent, ne, point_limit)
            else:
                # data_splitby_plane
                print('hello')
            # Check that you are going to grow two branches
            if split_data['ss'][0] and split_data['ss'][1]:
                for n in range(0, 2):
                    com = mesh_com(ne + 1, ld, datapoints)

                    start_node_loc = node_loc[np_start][1:4]
                    length_new = np.linalg.norm(fraction * (com - start_node_loc))
                    if length_new <  min_length:
                        end_node_loc = start_node_loc + min_length * (com - start_node_loc)/np.linalg.norm((com - start_node_loc))
                    else:
                        end_node_loc = start_node_loc + fraction * (com - start_node_loc)
                    if sorv is 'surface':
                        end_node_loc[2] = pg_utilities.z_from_xy(end_node_loc[0], end_node_loc[1], x_radius, y_radius,
                                                                 z_radius)
                    # insert checks that the branches are valid here
                    branch = True
                    if sorv is 'surface':
                        node1=np.array([node_loc[elems[ne_parent][1]][1], node_loc[elems[ne_parent][1]][2], 0])
                        node2=np.array([start_node_loc[0], start_node_loc[1], 0])
                        node3=np.array([end_node_loc[0], end_node_loc[1], 0])
                        end_node=mesh_check_angle(angle_min,angle_max,node1,node2,node3,ne_parent,ne+1)
                        end_node_loc[0:2] = end_node[0:2]
                        end_node_loc[2] = pg_utilities.z_from_xy(end_node[0],end_node[1],x_radius,y_radius,z_radius)
                    # Create new elements and nodes
                    elems[ne + 1][0] = ne + 1
                    elems[ne + 1][1] = np_start
                    elems[ne + 1][2] = nnod + 1
                    elem_upstream[ne + 1][0] = 1
                    elem_upstream[ne + 1][1] = ne_parent
                    elem_downstream[ne + 1][0] = 0

                    elem_downstream[ne_parent][0] = elem_downstream[ne_parent][0] + 1
                    elem_downstream[ne_parent][elem_downstream[ne_parent][0]] = ne + 1

                    node_loc[nnod + 1][0] = nnod + 1
                    node_loc[nnod + 1][1:4] = end_node_loc
                    ne = ne + 1
                    nnod = nnod + 1

                    elem_order[ne][0] = elem_order[ne][0] + 1
                    noelem_gen = noelem_gen + 1 #Only used for runtime output

                    nstem[ne][0] = nstem[ne_parent][0]

                    if branch and split_data['ss'][n]:
                        ne_temp[n_elm] = ne
                        n_elm = n_elm + 1
                    else:
                        numtb = numtb + 1


            else:  # Not ss so no splitting
                # Not enough seed points in the set during the split so ne_parent becomes a terminal branch
                # Find teh closest seed point and remove
                numtb = numtb + 1
                min_dist = 1.0e10
                for nd in range(0, len(datapoints)):
                    if ld[nd] != 0:
                        if ld[nd] == ne + 1:
                            ld[nd] = ne_parent #give back to parent
                        if ld[nd] == ne + 2:
                            ld[nd] = ne_parent
                        dist = dist_two_vectors(datapoints[nd][:], node_loc[elems[ne_parent][2]][1:4])
                        if dist < min_dist:
                            if ld[nd] == ne_parent:
                                nd_min = nd
                                min_dist = dist
                ldtmp1[nd_min] = ne_parent
                ld[nd_min] = 0
                ld_num = ld_num - 1

        # .......Copy the temporary list of branches to NE_OLD. These become the
        # .......parent elements for the next branching

        for n in range(0, n_elm):
            ne_old[n] = ne_temp[n]
            nstem[ne_old[n]][1] = 0  # initialaw count of data points

        # reallocate datapoints
        ld = data_to_mesh(ld, datapoints, ne_old[0:n_elm], node_loc, elems)


        print('   ' +str(ngen) + '   ' + str(noelem_gen) + '   ' + str(ne) +
              '   ' + str(numtb) + '   ' + str(ld_num) + '   ' + str(
            numtb + ld_num + numtbsum))
        numtbsum = numtbsum + numtb

    elems.resize(ne + 1, 3, refcheck=False)
    elem_upstream.resize(ne + 1, 3, refcheck=False)
    elem_downstream.resize(ne + 1, 3, refcheck=False)
    node_loc.resize(nnod + 1, 4, refcheck=False)

    return {'nodes': node_loc, 'elems': elems, 'elem_up': elem_upstream, 'elem_down': elem_downstream}

def mesh_check_angle(angle_min,angle_max,node1,node2,node3,ne_parent,myno):
    normal_to_plane = np.zeros(3)
    vector_cross_n = np.zeros(3)
    vector1=(node2-node1)
    vector1_u=vector1/np.linalg.norm(node2-node1)
    vector2=(node3-node2)
    vector2_u=vector2/np.linalg.norm(node3-node2)

    if (np.equal(vector1_u,vector2_u)).all():
        node3 [0] = node3[0]*.99
        node3 [1] = node3[1]*1.01
        vector2=node3-node2


    #want to rotate vector 2 wrt vector 1
    angle=pg_utilities.angle_two_vectors(vector1,vector2)


    normal_to_plane [0] = (vector2[1]*vector1[2] - vector2[2]*vector1[1])
    normal_to_plane [1] = (vector2[0]*vector1[2] - vector2[2]*vector1[0])
    normal_to_plane [2] = (vector2[0]*vector1[1] - vector2[1]*vector1[0])

    normal_to_plane_u = normal_to_plane / np.linalg.norm(normal_to_plane)

    vector_cross_n [0] = (vector1[1]*normal_to_plane[2] - vector1[2]*normal_to_plane[1])
    vector_cross_n [1] = (vector1[0]*normal_to_plane[2] - vector1[2]*normal_to_plane[0])
    vector_cross_n [2] = (vector1[0]*normal_to_plane[1] - vector1[1]*normal_to_plane[0])

    dotprod=np.dot(vector2,vector_cross_n)

    if dotprod < 0:
        angle = -1* angle

    if abs(angle) < angle_min:
        # Need to adjust node3 to get a angle equal to  angle_min
        angle0=abs(angle)
        angle_rot=(angle0-angle_min) #amount of angle to add

        #need to rotate around axis normal to the plane that the original and the vector make
        #unit vector normal to plane:
        R=np.zeros((3,3))
        R[0][0] = np.cos(angle_rot) + normal_to_plane_u[0]**2*(1-np.cos(angle_rot))
        R[0][1] = normal_to_plane_u[0]*normal_to_plane_u[1]*(1-np.cos(angle_rot))-normal_to_plane_u[2]*np.sin(angle_rot)
        R[0][2] = normal_to_plane_u[0]*normal_to_plane[2]*(1-np.cos(angle_rot))+normal_to_plane_u[1]*np.sin(angle_rot)
        R[1][0] =normal_to_plane_u[0]*normal_to_plane_u[1]*(1-np.cos(angle_rot))+normal_to_plane_u[2]*np.sin(angle_rot)
        R[1][1] = np.cos(angle_rot) + normal_to_plane_u[1]**2*(1-np.cos(angle_rot))
        R[1][2] = normal_to_plane_u[1]*normal_to_plane_u[2]*(1-np.cos(angle_rot))-normal_to_plane_u[0]*np.sin(angle_rot)
        R[2][0] = normal_to_plane_u[0]*normal_to_plane[2]*(1-np.cos(angle_rot))-normal_to_plane_u[1]*np.sin(angle_rot)
        R[2][1] = normal_to_plane_u[1]*normal_to_plane_u[2]*(1-np.cos(angle_rot))+normal_to_plane_u[0]*np.sin(angle_rot)
        R[2][2] = np.cos(angle_rot) + normal_to_plane_u[2]**2*(1-np.cos(angle_rot))
        nu_vec=np.zeros(3)
        nu_vec [0] = R[0][0]*vector2[0]+R[0][1]*vector2[1]+R[0][2]*vector2[2]
        nu_vec[1] = R[1][0] * vector2[0] + R[1][1] * vector2[1] + R[1][2] * vector2[2]
        nu_vec[2] = R[2][0] * vector2[0] + R[2][1] * vector2[1] + R[2][2] * vector2[2]

        node3 = node2 + nu_vec

    elif abs(angle) > angle_max:
        # Need to adjust node3 to get a angle equal to  angle_max
        angle0=abs(angle)
        angle_rot=(angle0-angle_max) #amount of angle to add

        #need to rotate around axis normal to the plane that the original and the vector make
        #unit vector normal to plane:
        R=np.zeros((3,3))
        R[0][0] = np.cos(angle_rot) + normal_to_plane_u[0]**2*(1-np.cos(angle_rot))
        R[0][1] = normal_to_plane_u[0]*normal_to_plane_u[1]*(1-np.cos(angle_rot))-normal_to_plane_u[2]*np.sin(angle_rot)
        R[0][2] = normal_to_plane_u[0]*normal_to_plane[2]*(1-np.cos(angle_rot))+normal_to_plane_u[1]*np.sin(angle_rot)
        R[1][0] =normal_to_plane_u[0]*normal_to_plane_u[1]*(1-np.cos(angle_rot))+normal_to_plane_u[2]*np.sin(angle_rot)
        R[1][1] = np.cos(angle_rot) + normal_to_plane_u[1]**2*(1-np.cos(angle_rot))
        R[1][2] = normal_to_plane_u[1]*normal_to_plane_u[2]*(1-np.cos(angle_rot))-normal_to_plane_u[0]*np.sin(angle_rot)
        R[2][0] = normal_to_plane_u[0]*normal_to_plane[2]*(1-np.cos(angle_rot))-normal_to_plane_u[1]*np.sin(angle_rot)
        R[2][1] = normal_to_plane_u[1]*normal_to_plane_u[2]*(1-np.cos(angle_rot))+normal_to_plane_u[0]*np.sin(angle_rot)
        R[2][2] = np.cos(angle_rot) + normal_to_plane_u[2]**2*(1-np.cos(angle_rot))

        nu_vec=np.zeros(3)
        nu_vec [0] = R[0][0]*vector2[0]+R[0][1]*vector2[1]+R[0][2]*vector2[2]
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
    nd1_1st = 0
    nd2_1st = 0
    nde_change1 = 0
    nde_change2 = 0
    for nd in range(0, len(datapoints)):
        nsp = ld[nd]

        if nsp == ne_parent:  # data point belongs to this element
            npoints = npoints + 1
            checkvalue = (x1[0] - x0[0]) * (datapoints[nd][1] - x0[1]) - (datapoints[nd][0] - x0[0]) * (x1[1] - x0[1])
            if checkvalue >= 0:
                if dat1 == 0:
                    nd1_1st = nd
                dat1 = dat1 + 1
                ld[nd] = ne_current + 1
                nde_change1 = ne_current + 1
            else:
                if dat2 == 0:
                    nd2_1st = nd
                dat2 = dat2 + 1
                ld[nd] = ne_current + 2
                nde_change2 = ne_current + 2
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

    return {'nde_change1': nde_change1, 'nde_change2': nde_change2, 'ss': ss}


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

    return {'umb_nodes': node_loc, 'umb_elems': elems, 'elem_up': elem_upstream, 'elem_down': elem_downstream}


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
