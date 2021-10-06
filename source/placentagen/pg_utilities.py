#!/usr/bin/env python
import numpy as np


def rotation_matrix_3d(axis, angle_rot):
    axis = axis / np.linalg.norm(axis)  # normalise just in case
    R = np.zeros((3, 3))
    R[0][0] = np.cos(angle_rot) + axis[0] ** 2 * (1 - np.cos(angle_rot))
    R[0][1] = axis[0] * axis[1] * (1 - np.cos(angle_rot)) - axis[2] * np.sin(angle_rot)
    R[0][2] = axis[0] * axis[2] * (1 - np.cos(angle_rot)) + axis[1] * np.sin(angle_rot)
    R[1][0] = axis[0] * axis[1] * (1 - np.cos(angle_rot)) + axis[2] * np.sin(angle_rot)
    R[1][1] = np.cos(angle_rot) + axis[1] ** 2 * (1 - np.cos(angle_rot))
    R[1][2] = axis[1] * axis[2] * (1 - np.cos(angle_rot)) - axis[0] * np.sin(angle_rot)
    R[2][0] = axis[0] * axis[2] * (1 - np.cos(angle_rot)) - axis[1] * np.sin(angle_rot)
    R[2][1] = axis[1] * axis[2] * (1 - np.cos(angle_rot)) + axis[0] * np.sin(angle_rot)
    R[2][2] = np.cos(angle_rot) + axis[2] ** 2 * (1 - np.cos(angle_rot))

    return R


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
    coord_check = (x / x_radius) ** 2. + (y / y_radius) ** 2. + (z / z_radius) ** 2.
    if coord_check < 1.0:
        in_ellipsoid = True

    return in_ellipsoid


def check_on_ellipsoid(x, y, z, x_radius, y_radius, z_radius):
    zero_tol = 1.e-10
    on_ellipsoid = False  # default to false
    coord_check = (x / x_radius) ** 2. + (y / y_radius) ** 2. + (z / z_radius) ** 2.
    if abs(coord_check - 1.0) < zero_tol:
        on_ellipsoid = True

    return on_ellipsoid


def check_in_on_ellipsoid(x, y, z, x_radius, y_radius, z_radius):
    zero_tol = 1.e-10
    in_ellipsoid = False  # default to false
    coord_check = (x / x_radius) ** 2. + (y / y_radius) ** 2. + (z / z_radius) ** 2.
    if coord_check < 1.0:
        in_ellipsoid = True
    elif abs(coord_check - 1.0) < zero_tol:
        in_ellipsoid = True

    return in_ellipsoid


def angle_two_vectors(vector1, vector2):
    vector1_u = vector1 / np.linalg.norm(vector1)
    vector2_u = vector2 / np.linalg.norm(vector2)

    if (np.equal(vector1_u, vector2_u)).all():  # vectors are parallel
        angle = 0.0
    elif (np.equal(vector1_u, -1.0 * vector2_u)).all():  # vectors are anti-parrallel.
        angle = np.pi
    else:
        dotprod = np.dot(vector1_u, vector2_u)
        if np.isclose(1.0, dotprod):
            # can't do arccos of 1
            angle = np.sqrt(2. * np.abs(1. - dotprod))  # small angle approximation to cos near theta =1
        else:
            angle = np.arccos(dotprod)

    return angle


def element_connectivity_1D(node_loc, elems):
    # Initialise connectivity arrays
    num_elems = len(elems)
    elem_upstream = np.zeros((num_elems, 3), dtype=int)
    elem_downstream = np.zeros((num_elems, 3), dtype=int)
    num_nodes = len(node_loc)
    elems_at_node = np.zeros((num_nodes, 4), dtype=int)
    # determine elements that are associated with each node
    for ne in range(0, num_elems):
        for nn in range(1, 3):
            nnod = elems[ne][nn]
            elems_at_node[nnod][0] = elems_at_node[nnod][0] + 1
            elems_at_node[nnod][elems_at_node[nnod][0]] = ne
    # assign connectivity
    for ne in range(0, num_elems):
        nnod2 = elems[ne][2]  # second node in elem
        for noelem in range(1, elems_at_node[nnod2][0] + 1):
            ne2 = elems_at_node[nnod2][noelem]
            if ne2 != ne:
                elem_upstream[ne2][0] = elem_upstream[ne2][0] + 1
                elem_upstream[ne2][elem_upstream[ne2][0]] = ne
                elem_downstream[ne][0] = elem_downstream[ne][0] + 1
                elem_downstream[ne][elem_downstream[ne][0]] = ne2

    return {'elem_up': elem_upstream, 'elem_down': elem_downstream}

def group_elem_parent(ne_parent, elem_downstream):

    ne_old = np.zeros(5000, dtype=int)
    ntemp_list = np.zeros(5000, dtype=int)
    ne_temp = np.zeros(5000, dtype=int)

    NT_BNS = 1 #Initialising to arbritary non-zero value, number of branches downstream of current
    ne_old[1] = ne_parent #first element
    ne_count = 0 #How many elements
    ntemp_list[ne_count] = ne_parent
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

    ntemp_list.resize(ne_count+1,refcheck=False)

    return ntemp_list


def plane_from_3_pts(x0, x1, x2, normalise):
    #    PLANE_FROM_3_PTS finds the equation of a plane in three
    #    dimensions and a vector normal to the plane from three
    #    non-colinear points.
    #    normalise = False for raw normal and plane equation
    #    normalise = True for unit normal and plane equation
    #    The coefficients represent aX + bY + cZ + d = 0
    #    NORML(1)=a,NORML(2)=b,NORML(3)=c,NORML(4)=d

    norml = np.zeros(4)
    diff1 = x1 - x0
    diff2 = x1 - x2


    norml[0] = diff1[1] * diff2[2] - diff1[2] * diff2[1]
    norml[1] = diff1[2] * diff2[0] - diff1[0] * diff2[2]
    norml[2] = diff1[0] * diff2[1] - diff1[1] * diff2[0]

    if normalise:
        norml[0:3] = norml[0:3] / np.linalg.norm(norml[0:3])

    norml[3] = 0.0
    for nj in range(0, 3):
        norml[3] = norml[3] - norml[nj] * x0[nj]

    return norml


def check_colinear(x0, x1, x2):
    colinear = False
    vector1 = (x1 - x0) / np.linalg.norm(x1 - x0)
    vector2 = (x1 - x2) / np.linalg.norm(x1 - x2)
    array_test1 = np.isclose(vector1, vector2)
    array_test2 = np.isclose(vector1, -1. * vector2)
    if all(array_test1) is True:
        colinear = True
    elif all(array_test2) is True:
        colinear = True

    return colinear


def samp_gr_for_node_loc(rectangular_mesh):
    # calculate parameter of generated rectangular mesh
    # Input:
    # - rectangular_mesh: generated mesh
    # Outputs:
    # - startx:starting point in x axis
    # - starty:starting point in y axis 
    # - startz:starting point in z axis
    # - xside:lenght of each el in x axis
    # - yside:lenght of each el in y axis
    # - zside:lenght of each el in z axis
    # - nelem_x:number of element in x axis
    # - nelem_y:number of element in y axis
    # - nelem_z:number of element in z axis

    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']
    startx = np.min(nodes[:, 0])
    xside = nodes[elems[0][8]][0] - nodes[elems[0][1]][0]
    endx = np.max(nodes[:, 0])
    nelem_x = (endx - startx) / xside
    starty = np.min(nodes[:, 1])
    yside = nodes[elems[0][8]][1] - nodes[elems[0][1]][1]
    endy = np.max(nodes[:, 1])
    nelem_y = (endy - starty) / yside
    startz = np.min(nodes[:, 2])
    zside = nodes[elems[0][8]][2] - nodes[elems[0][1]][2]
    endz = np.max(nodes[:, 2])
    nelem_z = (endz - startz) / zside

    return startx, starty, startz, xside, yside, zside, nelem_x, nelem_y, nelem_z




def sort_elements(v1, v2):
    ######
    # Function: takes a list of node pairs (v1, v2) and creates a list of nodes and elements
    # Inputs: v1 - N x 3 array of start node coordinates
    #         v2 - N x 3 array of end node coordinates
    # Outputs: nodes - an M x 3 array giving cartesian coordinates (x,y,z) for the node locations in the tree
    #          elems - an N x 3 array, the first colum in the element number, the second two columns are the index of the start and end node
    #          elements that start and end at the same node are given a value of -1
    ######
    Nelem = len(v1)
    print("number of elements",Nelem)
    elems = np.zeros([Nelem, 3])
    nodes = np.zeros([Nelem * 2, 3])  # max number of nodes possible

    iN = 0  # node index

    # go through first node list
    for iE in range(0, Nelem):

        v = v1[iE, :]
        index = is_member(v, nodes[0:iN][:])  # see if the node is in the nodes list

        if index == -1:  # if not, create a new node
            nodes[iN, :] = v
            index = iN
            iN = iN + 1

        # else just use index of existing node
        elems[iE, 1] = int(index)
        elems[iE, 0] = int(iE)  # first column of elements is just the element number

    # go through second node list
    for iE in range(0, Nelem):

        v = v2[iE, :]
        index = is_member(v, nodes[0:iN, :])

        if index == -1:
            nodes[iN, :] = v
            index = iN
            iN = iN + 1

        elems[iE, 2] = int(index)

        if elems[iE][1] == elems[iE][2]:
            elems[iE, 0:2] = int(-1)

    nodes = nodes[0:iN:1][:]  # truncate based on how many nodes were actually assigned

    return (elems, nodes)



def locate_node(startx, starty, startz, xside, yside, zside, nelem_x, nelem_y, nelem_z, coord_node):
    # calculate where a give point/node is located in a given rectangular mesh
    # Inputs:
    # - startx:starting point in x axis
    # - starty:starting point in y axis 
    # - startz:starting point in z axis
    # - xside:lenght of each el in x axis
    # - yside:lenght of each el in y axis
    # - zside:lenght of each el in z axis
    # - nelem_x:number of element in x axis
    # - nelem_y:number of element in y axis
    # - nelem_z:number of element in z axis
    # - coord_node: the node that needs to be located
    # Outputs:
    # - nelem:number of element in rectangular mesh where the node/point is located  
    xelem_num = np.floor((coord_node[0] - startx) / xside)
    if xelem_num == int(nelem_x):
        xelem_num = xelem_num - 1
    yelem_num = np.floor((coord_node[1] - starty) / yside)
    if yelem_num == int(nelem_y):
        yelem_num = yelem_num - 1
    zelem_num = np.floor((coord_node[2] - startz) / zside)
    if zelem_num == int(nelem_z):
        zelem_num = zelem_num - 1
    nelem = int(xelem_num + yelem_num * nelem_x + zelem_num * (
                nelem_x * nelem_y))  # this is the element where the point/node located
    return nelem

def renumber_geom(nodes,elems):
    #renumbers nodes and elements in case of skipped nodes in external editing
    nodes_list = nodes['nodes'][:,0].astype(int)
    elems_temp = elems['elems']
    total_nodes = nodes['total_nodes']
    total_elems = elems['total_elems']
    nod_array = np.copy(nodes['nodes'])
    el_array = np.copy(elems_temp)
    total_el = 0
    for ne in range(0,total_elems):
        new_np1 = np.where(nodes_list == elems_temp[ne,1])[0][0]
        new_np2 = np.where(nodes_list == elems_temp[ne,2])[0][0]

        if new_np1 != new_np2:#removing accidental node to node connections
            #np.where(el_array[:,1]==new_np1 and el_array[:,2]==new_np2))
            el_array[total_el,0] = total_el
            if new_np1 <   new_np2:
                el_array[total_el,1] = new_np1
                el_array[total_el, 2] = new_np2
            else:
                el_array[total_el,1] = new_np2
                el_array[total_el, 2] = new_np1
            total_el = total_el + 1
            unq, count = np.unique(el_array[0:total_el,1:3], axis=0, return_counts=True)
            repeated_groups = unq[count > 1]
            if len(repeated_groups) > 0:
                print('removing duplicate')
                el_array[total_el, :] = 0
                total_el = total_el - 1
        else:
            print('removing isolated')


    for nnod in range(0,total_nodes):
        nod_array[nnod,0] = nnod

    return {'total_elems': total_el, 'elems': el_array[0:total_el,:], 'total_nodes':total_nodes, 'nodes': nod_array}

def fix_branch_direction(first_node,elems_at_node,elems,seen_elements,branch_id,branches,old_parent_list,inlet_branch):
    new_parent_list = np.zeros(2,dtype = int)
    continuing = False
    elem = elems_at_node[first_node][1]
    connected_elems_no = elems_at_node[first_node][0]  # number of elements connected to this one
    branch_starts_at = first_node
    loop_parent = len(elems)+1
    while connected_elems_no ==2 or inlet_branch: #continuing branch
        inlet_branch = False
        if first_node in np.asarray(old_parent_list) and first_node != branch_starts_at:
            connected_elems_no=1
            loop_parent = first_node
        else:
            for i in range(0, connected_elems_no):
                elem = elems_at_node[first_node][i + 1]  # elements start at column index 1
                if not seen_elements[elem]:
                    branch_id[elem] = branches
                    if elems[elem][1] != first_node:
                        # swap nodes
                        elems[elem][2] = elems[elem][1]
                        elems[elem][1] = first_node
                    seen_elements[elem] = True
                    first_node = elems[elem][2]
                    connected_elems_no = elems_at_node[first_node][0]  # number of elements connected to this one
                    if connected_elems_no == 3:
                        new_parents = 0
                        for i in range(0, connected_elems_no):
                            elem = elems_at_node[first_node][i + 1]  # elements start at column index 1
                            if not seen_elements[elem]:
                                new_parent_list[new_parents] = elem
                                new_parents = new_parents + 1
                                branch_id[elem] = branches + new_parents
                                if elems[elem][1] != first_node:
                                    # swap nodes
                                    elems[elem][2] = elems[elem][1]
                                    elems[elem][1] = first_node
                            seen_elements[elem] = True
                            continuing = True
                    if connected_elems_no == 1:
                        continuing = False
                        break

    return new_parent_list,continuing,loop_parent,elem

def fix_elem_direction(inlet_node,nodes,elems):
    # populate the elems_at_node array listing the elements connected to each node
    num_nodes = len(nodes)
    num_elems = len(elems)
    elems_at_node = np.zeros((num_nodes, 10), dtype=int)
    for i in range(0, num_elems):
        elems_at_node[elems[i][1]][0] = elems_at_node[elems[i][1]][0] + 1
        j = elems_at_node[elems[i][1]][0]
        elems_at_node[elems[i][1]][j] = elems[i][0]
        elems_at_node[elems[i][2]][0] = elems_at_node[elems[i][2]][0] + 1
        j = elems_at_node[elems[i][2]][0]
        elems_at_node[elems[i][2]][j] = elems[i][0]

    for i in range(0,num_nodes):
        if np.all(nodes[i,1:4]== inlet_node):
            first_node = i
            print("FOUND FIRST NODE",i)
    seen_elements = np.zeros((num_elems), dtype=bool)
    branch_id = np.zeros((num_elems),dtype = int)
    #first_node = inlet_node
    branches = 1
    continuing = True
    old_parent_list = first_node
    loop_list = np.zeros(1,dtype=int)
    loop_list[0] = (num_elems+1)
    branch_start = []
    branch_end = []

    branch_start = np.append(branch_start, elems_at_node[first_node,1])

    [new_parent_list,continuing,loop_parent,branch_end_elem] = fix_branch_direction(first_node, elems_at_node, elems, seen_elements,branch_id,branches,old_parent_list,True)
    branch_end = np.append(branch_end, branch_end_elem)
    while len(new_parent_list)>0:
        if len(new_parent_list) > 0:
            new_parent_list2 = []
            for parent in range(0,len(new_parent_list)):
                second_node = elems[new_parent_list[parent],2]
                if second_node not in loop_list:
                    branches = branches + 1
                    branch_start  = np.append(branch_start,new_parent_list[parent])
                    branch_id[new_parent_list[parent]] = branches
                    [branch_list, continuing,loop_parent,branch_end_elem] = fix_branch_direction(second_node, elems_at_node, elems, seen_elements,
                                                           branch_id, branches,elems[new_parent_list,2],False)
                    branch_end = np.append(branch_end, branch_end_elem)
                else:
                    print('the parent loop',second_node,loop_list)

                if loop_parent< num_elems:
                    loop_list = np.append(loop_list, [loop_parent], axis=0)
                if continuing:
                    new_parent_list2 = np.append(new_parent_list2,branch_list,axis = 0)

            if(len(new_parent_list2)>0):
                new_parent_list = new_parent_list2.astype(int)
            else:
                new_parent_list = []
            print('new generation',new_parent_list)

    return elems,branch_id,branch_start,branch_end

def remove_rows(main_array, arrays):
    ######
    # Function: Remove rows from both mainArray and Arrays at which main array has values less than zero
    # Inputs: mainArray - an N x M array of values
    #         Arrays - a list of arrays each with length N for their first axis
    # Outputs: for each row of mainArray for which the first element is below zero; this row is removed from mainArray and from each array
    ######

    i = 0

    while i < len(main_array):
        if main_array[i, 0] < 0:  # then get rid of row from all arrays

            for j in range(0, len(arrays)):
                array = arrays[j]
                array = np.delete(array, (i), axis=0)
                arrays[j] = array
            main_array = np.delete(main_array, (i), axis=0)

        else:
            i = i + 1

    return main_array, arrays


######
# Function: Swaps 2 rows in an array
# Inputs: array - a N x M array
#         row1 & row2 - the indices of the two rows to be swapped
# Outputs: array, with row1 and row2 swapped
######

def row_swap_2d(array, row1, row2):
    placeholder = np.copy(array[row1, :])
    array[row1, :] = array[row2, :]
    array[row2, :] = placeholder
    return array


######
# Function: Swaps 2 rows in an array
# Inputs: array - a N x 1 array
#         row1 & row2 - the indices of the two rows to be swapped
# Outputs: array, with row1 and row2 swapped
######

def row_swap_1d(array, row1, row2):
    placeholder = np.copy(array[row1])
    array[row1] = array[row2]
    array[row2] = placeholder
    return array


######
# Function: Finds first occurrence of a specified row of values in an array or returns -1 if the given row is not present
#           Similar to Matlab isMember function
# Inputs: matrix - an N x M array
#         v - a 1 x M array
# Outputs: index at which v first occurs in matrix, or else -1
######

def is_member(v, matrix):
    L = (np.shape(matrix))
    L = L[0]

    for i in range(0, L):
        if np.array_equal(v, matrix[i, :]):
            index = i
            return index
    return -1


######
# Function: Creates a 3D plot of branching tree
# Inputs: nodes - an M x 3 array giving cartesian coordinates (x,y,z) for the node locations in the tree
#         elems - an N x 3 array, the first colum in the element number, the second two columns are the index of the start and end node
#         colour - an N x 1 array where value determines colour of corresponding element
#         Nc - the maximum number of elements connected at a single node
# Outputs: 3D plot of tree, with radius proportional to radii and colour depending on the input array
######

def plot_vasculature_3d(nodes, elems, colour, radii):
    # initialize arrays
    Ne = len(elems)
    elems = elems[:, 1:3]
    x = np.zeros([Ne, 2])
    y = np.zeros([Ne, 2])
    z = np.zeros([Ne, 2])

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # scale colour and radii
    colour = (colour - min(colour)) / max(colour) * 255.
    radii = radii / max(radii) * 3.

    for i in range(0, Ne):
        # get start and end node
        nN1 = int(elems[i, 0])
        nN2 = int(elems[i, 1])

        # get coordinates of nodes
        x[i, 0] = nodes[nN1, 0]
        y[i, 0] = nodes[nN1, 1]
        z[i, 0] = nodes[nN1, 2]
        x[i, 1] = nodes[nN2, 0]
        y[i, 1] = nodes[nN2, 1]
        z[i, 1] = nodes[nN2, 2]

        colour_value = np.asarray(cm.jet(int(colour[i])))
        ax.plot(np.squeeze(x[i, :]), np.squeeze(y[i, :]), np.squeeze(z[i, :]), c=colour_value[0:3], linewidth=radii[i])

    plt.show()

    return 0


######
# Function: Finds the maximum number of elements that join at one node
# Inputs: elems - an N x 3 array containing element number in the first column and node indices in the second two columns
# Outputs: Nc - the maximum number of elements that join at one node
######

def find_maximum_joins(elems):
    elems = np.concatenate([np.squeeze(elems[:, 1]), np.squeeze(elems[:, 2])])
    elems = elems.astype(int)
    result = np.bincount(elems)
    Nc = (max(result)) + 1
    plt.plot(result)
    plt.show()
    # Warning if detect an unusual value
    if Nc > 12:
        print('Warning, large number of elements at one node: ' + str(Nc))
        Nc = 12

    return Nc



def find_strahler_ratio(Orders, Factor):

    x = Orders
    yData = np.log(Factor)
    #plt.plot(x, yData, 'k--', linewidth=1.5, label='Data')

    # fit line to data
    xFit = np.unique(Orders)
    yFit = np.poly1d(np.polyfit(x, yData, 1))(np.unique(x))
    #plt.plot(np.unique(x), yFit, label='linear fit')

    # Scaling Coefficient is gradient of the fit
    grad = (yFit[len(yFit) - 1] - yFit[0]) / (xFit[len(xFit) - 1] - xFit[0])
    grad=np.abs(grad)
    grad=np.exp(grad)

    # R^2 value
    yMean = [np.mean(yData) for y in yData]
    r2 = 1 - (sum((yFit - yData) * (yFit - yData)) / sum((yMean - yData) * (yMean - yData)))

    heading = ('Strahler Ratio = ' + str(grad))
    #plt.title(heading)
    #plt.legend()
    #plt.show()
    return grad, r2


######
# Function: Remove rows from both mainArray and Arrays at which main array has values less than zero
# Inputs: mainArray - an N x M array of values
#         Arrays - a list of arrays each with length N for their first axis
# Outputs: for each row of mainArray for which the first element is below zero; this row is removed from mainArray and from each array
######

def remove_rows(main_array, arrays):
    i = 0

    while i < len(main_array):
        if main_array[i, 0] < 0:  # then get rid of row from all arrays

            for j in range(0, len(arrays)):
                array = arrays[j]
                array = np.delete(array, (i), axis=0)
                arrays[j] = array
            main_array = np.delete(main_array, (i), axis=0)

        else:
            i = i + 1

    return main_array, arrays

######
# Function: Swaps 2 rows in an array
# Inputs: array - a N x M array
#         row1 & row2 - the indices of the two rows to be swapped
# Outputs: array, with row1 and row2 swapped
######

def row_swap_2d(array, row1, row2):
    placeholder = np.copy(array[row1, :])
    array[row1, :] = array[row2, :]
    array[row2, :] = placeholder
    return array

######
# Function: Swaps 2 rows in an array
# Inputs: array - a N x 1 array
#         row1 & row2 - the indices of the two rows to be swapped
# Outputs: array, with row1 and row2 swapped
######

def row_swap_1d(array, row1, row2):
    placeholder = np.copy(array[row1])
    array[row1] = array[row2]
    array[row2] = placeholder
    return array