#!/usr/bin/env python
import numpy as np
from . import pg_utilities
import sys
from numpy import matlib

"""
.. module:: analyse_tree
  :synopsis: One sentence synopis (brief) could appear in module index.

:synopsis:A longer synopsis that could appear on the home page for that module in documentation.

"""


def calc_terminal_branch(node_loc, elems):
    """ What this function does

    Inputs:
       - input name: A description of the input

    Returns:
       - output name: a description of the output
            - can do sub points if there are complicated arrays


    A way you might want to use me is:

    >>> include a simple example of how to call the function

    Tell us what the function will do
    """
    # This function generates a list of terminal nodes associated with a branching geometry
    # inputs are node locations and elements
    num_elems = len(elems)
    num_nodes = len(node_loc)
    elem_cnct = pg_utilities.element_connectivity_1D(node_loc, elems)

    terminal_branches = np.zeros(num_elems, dtype=int)
    terminal_nodes = np.zeros(num_nodes, dtype=int)

    num_term = 0
    for ne in range(0, num_elems):
        if elem_cnct['elem_down'][ne][0] == 0:  # no downstream element
            terminal_branches[num_term] = ne
            terminal_nodes[num_term] = elems[ne][2]  # node that leaves the terminal element
            num_term = num_term + 1

    terminal_branches = np.resize(terminal_branches, num_term)
    terminal_nodes = np.resize(terminal_nodes, num_term)

    print('Total number of terminals assessed, num_terminals =  ' + str(num_term))

    return {'terminal_elems': terminal_branches, 'terminal_nodes': terminal_nodes, 'total_terminals': num_term}


def evaluate_orders(node_loc, elems):
    # calculates generations, Horsfield orders, Strahler orders for a given tree
    # Works for diverging trees only
    # Inputs are:
    # node_loc = array with location of nodes
    # elems = array with location of elements
    num_elems = len(elems)
    # Calculate connectivity of elements
    elem_connect = pg_utilities.element_connectivity_1D(node_loc, elems)
    elem_upstream = elem_connect['elem_up']
    elem_downstream = elem_connect['elem_down']
    # Initialise order definition arrays
    strahler = np.zeros(len(elems), dtype=int)
    horsfield = np.zeros(len(elems), dtype=int)
    generation = np.zeros(len(elems), dtype=int)

    # Calculate generation of each element
    maxgen = 1  # Maximum possible generation
    for ne in range(0, num_elems):
        ne0 = elem_upstream[ne][1]
        if ne0 != 0:
            # Calculate parent generation
            n_generation = generation[ne0]
            if elem_downstream[ne0][0] == 1:
                # Continuation of previous element
                generation[ne] = n_generation
            elif elem_downstream[ne0][0] >= 2:
                # Bifurcation (or morefurcation)
                generation[ne] = n_generation + 1
        else:
            generation[ne] = 1  # Inlet
        maxgen = np.maximum(maxgen, generation[ne])

    # Now need to loop backwards to do ordering systems
    
    for ne in range(num_elems - 1, -1, -1):
        
        n_horsfield = np.maximum(horsfield[ne], 1)
        n_children = elem_downstream[ne][0]
        if n_children == 1:
            if generation[elem_downstream[ne][1]] == 0:
                n_children = 0
        temp_strahler = 0
        strahler_add = 1
        if n_children >= 2:  # Bifurcation downstream
            temp_strahler = strahler[elem_downstream[ne][1]]  # first daughter
            for noelem in range(1, n_children + 1):
                ne2 = elem_downstream[ne][noelem]
                temp_horsfield = horsfield[ne2]
                if temp_horsfield > n_horsfield:
                    n_horsfield = temp_horsfield
                if strahler[ne2] < temp_strahler:
                    strahler_add = 0
                elif strahler[ne2] > temp_strahler:
                    strahler_add = 0
                    temp_strahler = strahler[ne2]  # strahler of highest daughter
            n_horsfield = n_horsfield + 1
        elif n_children == 1:
            ne2 = elem_downstream[ne][1]  # element no of daughter
            n_horsfield = horsfield[ne2]
            strahler_add = strahler[ne2]
        horsfield[ne] = n_horsfield
        strahler[ne] = temp_strahler + strahler_add
        
    return {'strahler': strahler, 'horsfield': horsfield, 'generation': generation}


def define_radius_by_order(node_loc, elems, system, inlet_elem, inlet_radius, radius_ratio):
    # This function defines radii in a branching tree by 'order' of the vessel
    # Inputs are:
    # node_loc: The nodes in the branching tree
    # elems: The elements in the branching tree
    # system: 'strahler','horsfield' or 'generation' to define vessel order
    # inlet_elem: element number that you want to define as having inlet_radius
    # inlet_radius: the radius of your inlet vessel
    # radius ratio: Strahler or Horsfield type ratio, defines the slope of log(order) vs log(radius)
    num_elems = len(elems)
    radius = np.zeros(num_elems)  # initialise radius array
    # Evaluate orders in the system
    orders = evaluate_orders(node_loc, elems)
    elem_order = orders[system]
    ne = inlet_elem
    n_max_ord = elem_order[ne]
    radius[ne] = inlet_radius

    for ne in range(0, num_elems):
        
        radius[ne] = 10. ** (np.log10(radius_ratio) * (elem_order[ne] - n_max_ord) + np.log10(inlet_radius))
      
    
    return radius


def tree_statistics(node_loc, elems, radius, orders):
    # Caclulates tree statistics for a given tree and prints to terminal
    # Inputs are:
    # node_loc: The nodes in the branching tree
    # elems: The elements in the branching tree
    # radius: per element radius
    # orders: per element order

    num_elems = len(elems)
    diameters = 2.0 * radius
    connectivity = pg_utilities.element_connectivity_1D(node_loc, elems)
    elem_upstream = connectivity['elem_up']
    elem_downstream = connectivity['elem_down']
    num_schemes = 3
    generation = orders['generation']
    horsfield = orders['horsfield']
    strahler = orders['strahler']

    index = np.zeros(3, dtype=int)

    # local arrays
    # length array
    lengths = np.zeros(num_elems)
    # ratios: index i is order, index j is ratios of branching (j=1), length (j=2), and diameter (j=3)

    nbranches = np.zeros((num_schemes + 1, num_elems))

    # j = 0 length, j = 1 diameter, j = 4 L/D
    branches = np.zeros((5, num_elems))

    for ne in range(0, num_elems):
        np1 = elems[ne][1]
        np2 = elems[ne][2]
        point1 = node_loc[np1][1:4]
        point2 = node_loc[np2][1:4]
        lengths[ne] = np.linalg.norm(point1 - point2)

    ntotal = 0
    num_dpp = 0
    num_llp = 0
    N = 1
    for ne in range(0, num_elems):
        num_upstream = elem_upstream[ne][0]
        ne0 = elem_upstream[ne][1]
        index[0] = generation[ne]
        index[1] = horsfield[ne]
        index[2] = strahler[ne]
        add = False
        if (num_upstream == 0):  # nothing upstream
            add = True
        elif generation[ne0] != index[0]:  # something upstream and not the
            add = True
        if add:
            N = N + 1
            for i in range(0, num_schemes):
                nbranches[i][N - 1] = index[i]
            if num_upstream != 0:
                nbranches[num_schemes][N - 1] = strahler[ne0]  # strahler order of parent
            else:
                nbranches[num_schemes][N - 1] = 0

        # Add length of all segments along branch, and calculate mean diameter
        n_segments = 1

        mean_diameter = diameter[ne]
        branches[0][N - 1] = lengths[ne]
        ne_next = ne
        while elem_downstream[ne_next][0] == 1:
            ne_next = elem_downstream[ne_next][1]
            branches[0][N - 1] = branches[0][N - 1] + lengths[ne_next]
            mean_diameter = mean_diameter + diameters[ne_next]
            n_segments = n_segments + 1
            print(n_segments)


def terminals_in_sampling_grid_fast(rectangular_mesh, terminal_list, node_loc):
    # This function counts the number of terminals in a sampling grid element, will only work with
    # rectangular mesh created as in generate_shapes.gen_rectangular_mesh
    # inputs are:
    # Rectangular mesh - the sampling grid
    # terminal_list - a list of terminal,
    # node_loc - location of nodes
    num_terminals = terminal_list['total_terminals']
    terminals_in_grid = np.zeros(len(rectangular_mesh['elems']), dtype=int)
    terminal_elems = np.zeros(num_terminals, dtype=int)
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

    for nt in range(0, num_terminals):
        coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
        xelem_num = np.floor((coord_terminal[0] - startx) / xside)
        yelem_num = np.floor((coord_terminal[1] - starty) / yside)
        zelem_num = np.floor((coord_terminal[2] - startz) / zside)
        nelem = int(xelem_num + (yelem_num) * nelem_x + (zelem_num) * (nelem_x * nelem_y))
        terminals_in_grid[nelem] = terminals_in_grid[nelem] + 1
        terminal_elems[nt] = nelem  # record what element the terminal is in
    return {'terminals_in_grid': terminals_in_grid, 'terminal_elems': terminal_elems}


def terminals_in_sampling_grid(rectangular_mesh, placenta_list, terminal_list, node_loc):
    # This function counts the number of terminals in a sampling grid element
    # inputs are:
    # Rectangular mesh - the sampling grid
    # terminal_list - a list of terminal,
    # node_loc - location of nodes
    num_sample_elems = len(placenta_list)
    num_terminals = terminal_list['total_terminals']
    terminals_in_grid = np.zeros(len(rectangular_mesh['elems']), dtype=int)
    terminal_mapped = np.zeros(num_terminals, dtype=int)
    terminal_elems = np.zeros(num_terminals, dtype=int)

    for ne_i in range(0, num_sample_elems):
        # First node has min x,y,z and last node has max x,y,z
        ne = placenta_list[ne_i]
        if placenta_list[ne_i] > 0:  # There is some placenta in this element (assuming none in el 0)
            first_node = rectangular_mesh['elems'][ne][1]
            last_node = rectangular_mesh['elems'][ne][8]
            min_coords = rectangular_mesh['nodes'][first_node][0:3]
            max_coords = rectangular_mesh['nodes'][last_node][0:3]
            for nt in range(0, num_terminals):
                if terminal_mapped[nt] == 0:
                    in_element = False
                    coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
                    if coord_terminal[0] >= min_coords[0]:
                        if coord_terminal[0] < max_coords[0]:
                            if coord_terminal[1] >= min_coords[1]:
                                if coord_terminal[1] < max_coords[1]:
                                    if coord_terminal[2] >= min_coords[2]:
                                        if coord_terminal[2] < max_coords[2]:
                                            in_element = True
                    if in_element:
                        terminals_in_grid[ne] = terminals_in_grid[ne] + 1
                        terminal_mapped[nt] = 1
                        terminal_elems[nt] = ne
   
    return {'terminals_in_grid': terminals_in_grid, 'terminal_elems': terminal_elems}


def ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, num_test_points):
    # This subroutine calculates the placental volume associated with each element in a samplling grid
    # inputs are:
    # rectangular_mesh = the sampling grid nodes and elements
    # volume = placental volume
    # thickness = placental thickness
    # ellipiticity = placental ellipticity
    # num_test_points = resolution of integration quadrature
    total_elems = rectangular_mesh['total_elems']
    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']
    
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    # Initialise the array that defines the volume of placenta in each grid element
    pl_vol_in_grid = np.zeros(total_elems)
    non_empty_loc = np.zeros(total_elems, dtype=int)
    non_empty_count = 0

    for ne in range(0, len(elems)):  # looping through elements
        count_in_range = 0
        nod_in_range = np.zeros(8, dtype=int)
        # define range of x, y , and z in the element
        startx = nodes[elems[ne][1]][0]
        endx = nodes[elems[ne][8]][0]
        starty = nodes[elems[ne][1]][1]
        endy = nodes[elems[ne][8]][1]
        startz = nodes[elems[ne][1]][2]
        endz = nodes[elems[ne][8]][2]
        for nod in range(1, 9):
            check_in_range = pg_utilities.check_in_ellipsoid(nodes[elems[ne][nod]][0], nodes[elems[ne][nod]][1],
                                                             nodes[elems[ne][nod]][2], x_radius, y_radius, z_radius)
            check_on_range = pg_utilities.check_on_ellipsoid(nodes[elems[ne][nod]][0], nodes[elems[ne][nod]][1],
                                                             nodes[elems[ne][nod]][2], x_radius, y_radius, z_radius)
            if check_in_range or check_on_range:
                count_in_range = count_in_range + 1
                nod_in_range[nod - 1] = 1
        if count_in_range == 8:  # if all 8 nodes are inside the ellipsoid
            non_empty_loc[non_empty_count] = ne
            non_empty_count = non_empty_count + 1
            pl_vol_in_grid[ne] = (endx - startx) * (endy - starty) * (
                    endz - startz)  # the placental vol in that samp_grid_el is same as vol of samp_grid_el
        elif count_in_range == 0:  # if all 8 nodes are outside the ellpsiod
            # since this samp_grid_el is completely outside, the placental vol is zero
            pl_vol_in_grid[ne] = 0
        else:  # if some nodes in and some nodes out, the samp_grid_el is at the edge of ellipsoid
            # Use trapezoidal quadrature to caculate the volume under the surface of the ellipsoid in each element
            non_empty_loc[non_empty_count] = ne
            non_empty_count = non_empty_count + 1
            # need to map to positive quadrant
            repeat = False
            if (startz < 0 and endz <= 0):
                # need to project to positive z axis
                startz = abs(nodes[elems[ne][8]][2])
                endz = abs(nodes[elems[ne][1]][2])
            elif (startz < 0 and endz > 0):
                # Need to split into components above and below the axis and sum the two
                startz = 0
                endz = abs(nodes[elems[ne][1]][2])
                startz_2 = 0
                endz_2 = nodes[elems[ne][8]][2]
                repeat = True
            xVector = np.linspace(startx, endx, num_test_points)
            yVector = np.linspace(starty, endy, num_test_points)
            xv, yv = np.meshgrid(xVector, yVector)
            zv = z_radius ** 2 * (1 - (xv / x_radius) ** 2 - (yv / y_radius) ** 2)
            for i in range(num_test_points):
                for j in range(num_test_points):
                    if zv[i, j] <= startz ** 2:
                        zv[i, j] = startz ** 2
                    zv[i, j] = np.sqrt(zv[i, j])
                    if zv[i, j] > endz:
                        zv[i, j] = endz
                    elif zv[i, j] < startz:
                        zv[i, j] = startz
            intermediate = np.zeros(num_test_points)
            for i in range(0, num_test_points):
                intermediate[i] = np.trapz(zv[:, i], xVector)
            Value1 = np.trapz(intermediate, yVector)
            pl_vol_in_grid[ne] = (Value1 - startz * (endx - startx) * (endy - starty))
            if repeat:
                xVector = np.linspace(startx, endx, num_test_points)
                yVector = np.linspace(starty, endy, num_test_points)
                xv, yv = np.meshgrid(xVector, yVector)
                zv = z_radius ** 2 * (1 - (xv / x_radius) ** 2 - (yv / y_radius) ** 2)
                for i in range(num_test_points):
                    for j in range(num_test_points):
                        if zv[i, j] <= startz_2 ** 2:
                            zv[i, j] = startz_2 ** 2
                        zv[i, j] = np.sqrt(zv[i, j])
                        if zv[i, j] > endz_2:
                            zv[i, j] = endz_2
                        elif zv[i, j] < startz_2:
                            zv[i, j] = startz_2
                intermediate = np.zeros(num_test_points)
                for i in range(0, num_test_points):
                    intermediate[i] = np.trapz(zv[:, i], xVector)
                Value1 = np.trapz(intermediate, yVector)
                pl_vol_in_grid[ne] = pl_vol_in_grid[ne] + (Value1 - startz_2 * (endx - startx) * (
                        endy - starty))

    print('Number of Non-empty cells: ' + str(non_empty_count))
    print('Total number of cells: ' + str(total_elems))
    non_empty_loc = np.resize(non_empty_loc, non_empty_count)

    return {'pl_vol_in_grid': pl_vol_in_grid, 'non_empty_rects': non_empty_loc}


def cal_br_vol_samp_grid(rectangular_mesh, branch_nodes, branch_elems,branch_radius, volume, thickness, ellipticity, start_elem):
    '''
    This subroutine is to:
    1. calculate total volume of branches in each samp_grid_el (to use when calculate vol_frac/porosity)
    2. calculate total daimeter variable of branches in each samp_grid_el(to use when calculate wt_diam)
    3. count number of branches in each samp_grid_el
    4. calculate the volume of individual branch and total vol of all branches in the whole tree
         
    '''

    #Define the resolution of cylinder for analysis
    num_points_xy = 8
    num_points_z = 8
    # Define information about sampling grid required to place data points in correct locations
    total_sample_elems = rectangular_mesh['total_elems']
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
    
    #Define the placental ellipsoid
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)  # calculate radii of ellipsoid
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    
    unit_cyl_points = np.zeros((num_points_xy*num_points_xy*num_points_z,3))
    #Define a cylinder of points of radius 1 and length 1
    x = np.linspace(-1,1,num_points_xy)
    y = np.linspace(-1,1,num_points_xy)
    num_accepted = 0
    for k in range(0,num_points_z+1):
        for i in range(0,num_points_xy):
            for j in range(0,num_points_xy):
                if(x[i]**2 + y[j]**2)<=1:
                    new_z = 1 / np.double(num_points_z) * k
                    unit_cyl_points[num_accepted][0] = x[i]
                    unit_cyl_points[num_accepted][1] = y[j]
                    unit_cyl_points[num_accepted][2] = new_z
                    num_accepted = num_accepted+1
    unit_cyl_points.resize(num_accepted, 3, refcheck=False)
    cyl_points = np.copy(unit_cyl_points)
    cylindervector = np.array([0.0,0.0,1.0])
    
    ###Define and initialise arrays to be populated

    #The volume of each branch
    vol_each_br = np.zeros(len(branch_elems))
    # Array for total volume of sampling grid in each element
    total_vol_samp_gr = np.zeros(total_sample_elems)
    # Array for diameter variable of sampling grid in each element (this variable is to be used for weighted diameter calculation)
    total_diameter_samp_gr = np.zeros(total_sample_elems)
    #initialise counters
    branch_count = 0
    volume_outside_ellipsoid = 0.0
    volume_inside_ellipsoid = 0.0

    for ne in range(start_elem,len(branch_elems)):#len(branch_elems)):  # looping for all branchs in tree
        
        node1 = branch_nodes[branch_elems[ne][1]][1:4]  # coor of start node of a branch element
        node2 = branch_nodes[branch_elems[ne][2]][1:4]  # coor of end node of a branch element
        node1in = pg_utilities.check_in_on_ellipsoid(node1[0], node1[1], node1[2], x_radius, y_radius, z_radius)
        node2in = pg_utilities.check_in_on_ellipsoid(node2[0], node2[1], node2[2], x_radius, y_radius, z_radius)
        
        if not node1in and not node2in:
            print('Warning, element ' + str(ne) + 'is not in ellipsoid, if this is not expected check your geometry')
            print('Skipping this element from analysis')
            continue
        elif not node1in or not node2in:
            print('Warning, element ' + str(ne) + 'has one node not in the ellipsoid.')
            print('The first node ' + str(node1) + ' is ' + srt(node1in) + ' (True means inside).')
            print('The second node ' + str(node2) + ' is ' + srt(node2in) + ' (True means inside).')
            print('Skipping this element from analysis')
            continue

        branch_vector = node2-node1
        r = branch_radius[ne]
        length = np.linalg.norm(branch_vector)
        vol_each_br[ne] = np.pi*length*r**2.0
        vol_per_point = vol_each_br[ne]/(np.double(num_accepted))
        
        cyl_points[:,0:2]=unit_cyl_points[:,0:2]*r
        cyl_points[:,2] = unit_cyl_points[:,2]*length

        desiredvector = branch_vector / np.linalg.norm(branch_vector)

        rotation_axis = np.cross(desiredvector,cylindervector)

        if np.linalg.norm(rotation_axis) == 0:#aligned
            if node2[2]-node1[2] < 0:
                cyl_points[:, 2] = -1.0*cyl_points[:, 2]
        else:
            angle = pg_utilities.angle_two_vectors(cylindervector,desiredvector)
            rotation_mat = pg_utilities.rotation_matrix_3d(rotation_axis,angle)
            cyl_points = np.array(np.matrix(cyl_points)*np.matrix(rotation_mat))

        cyl_points[:, 0] = cyl_points[:, 0] + node1[0]
        cyl_points[:, 1] = cyl_points[:, 1] + node1[1]
        cyl_points[:, 2] = cyl_points[:, 2] + node1[2]

        # Array for vol distribution of inidvidual branch (not total)
        vol_distribution_each_br=np.zeros(total_sample_elems, dtype=float)
               
        for nt in range(0, num_accepted):
        
            coord_point = cyl_points[nt][0:3]
            inside=pg_utilities.check_in_on_ellipsoid(coord_point[0], coord_point[1], coord_point[2], x_radius, y_radius, z_radius)
            if inside:
                xelem_num = np.floor((coord_point[0] - startx) / xside)
                yelem_num = np.floor((coord_point[1] - starty) / yside)
                zelem_num = np.floor((coord_point[2] - startz) / zside)
                nelem = int(xelem_num + (yelem_num) * nelem_x + (zelem_num) * (nelem_x * nelem_y))
                total_vol_samp_gr[nelem] = total_vol_samp_gr[nelem] + vol_per_point
                vol_distribution_each_br[nelem]=vol_distribution_each_br[nelem]+vol_per_point
                volume_inside_ellipsoid = volume_inside_ellipsoid + vol_per_point
            else:
                #Data points lie outside the ellipsoid - this is OK in some cases, so the code shouldn't exit. However,
                #users should be able to check how much is outside of ellipsoid if they believe their branching geometry
                #is set up NOT to go outside the ellipsoid at all.
                volume_outside_ellipsoid = volume_outside_ellipsoid + vol_per_point
        
        total_diameter_samp_gr = total_diameter_samp_gr + vol_distribution_each_br*2*r#this variable is calculated as summation of diameter * vol of branch in grid (to be used for weight_diam)
    
    percent_outside = volume_outside_ellipsoid/np.sum(total_vol_samp_gr)*100.0
    
    print('Analysis complete ' + str(percent_outside) + '% of analysed points lie outside the ellipsoid.')
    print('Total branch volume analysed ' + str(volume_outside_ellipsoid + np.sum(total_vol_samp_gr)) + ' (' + str(np.sum(vol_each_br)) +')')
    
    return {'br_vol_in_grid': total_vol_samp_gr,'br_diameter_in_grid':total_diameter_samp_gr}



def terminal_villous_volume(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute):
    """ This function calculates the average volume of a terminal villous based on structural
    characteristics measured in the literature.

    Inputs:
       - num_int_gens: Number of generations of intermediate villous per terminal 'stem' villois
       - num_convolutes: Number of terminal convolutes per intermediate villous
       - len_int: Length of a typical intermediate villous
       - rad_int: Radius of a typical intermediate villous
       - len_convolute: Length of a typical terminal convolute
       - rad_convolute: Radius of a typical terminal convolute

    Returns:
       - term_vill_volume: Typical volume of a terminal villous


    A way you might want to use me is:

    >>> num_int_gens = 3
    >>> num_convolutes = 10
    >>> len_int = 1.5 #mm
    >>> rad_int = 0.03 #mm
    >>> len_convolute = 3.0 #mm
    >>> rad_convolute = 0.025 #mm
    >>> terminal_villous_volume(num_int_gens,num_convolutes,len_int,rad_int,len_convulute,rad_convolute)

    This will take the normal average data from Leiser et al (1990, IBBN:3805554680) and calculate
    average volume of terminal villi to be ~1.77 mm^3
    """

    #Each terminal stem villous branches to two immature intermediate villi
    #and then to three generations of mature intermediate villi each with ~10 terminal conduits
    num_ints = 1
    term_vill_volume = 0.0
    for i in range(0,4):
        num_ints = num_ints*2.0
        vol_ints = num_ints*np.pi*len_int*rad_int**2.0
        if i > 0:
            vol_convolutes = num_ints*num_convolutes*np.pi*len_convolute*rad_convolute**2.0
        else:
            vol_convolutes = 0.0
        term_vill_volume = term_vill_volume + vol_ints + vol_convolutes
    
    return term_vill_volume

def tissue_vol_in_samp_gr(term_vill_volume,br_vol_in_grid,terminals_in_grid):

    tissue_vol = br_vol_in_grid +term_vill_volume*terminals_in_grid
 
    return tissue_vol


def terminal_villous_diameter(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute):
  
    """
    The concept to calculate terminal villous diameter follows the same as terminal_villous_volume calculation
    We multiply vol of each branch with diameter of each branch and summation of them to be able to calculate the weighted_diameter in the next subroutine
    """
    num_ints = 1
    term_vill_diameter = 0.0
    for i in range(0,4):
        num_ints = num_ints*2.0
        diameter_ints = num_ints*(np.pi*len_int*rad_int**2.0)*2*rad_int
        if i > 0:
            diameter_convolutes = num_ints*num_convolutes*(np.pi*len_convolute*rad_convolute**2.0)*2*rad_convolute
        else:
            diameter_convolutes = 0.0
        term_vill_diameter = term_vill_diameter + diameter_ints + diameter_convolutes
      
    return term_vill_diameter


def weighted_diameter_in_samp_gr(term_vill_diameter,br_diameter_in_grid,terminals_in_grid,tissue_vol):
    """
    Weighted_diameter is calculated as:

    weighted_diameter each sampling grid = (d1*v1+d2*v2+d3*v3+...+dn*vn)/(v1+v2+v2+...+vn)
    The numerater is calculated by summation of br_diameter_in_grid (output from cal_br_vol_samp_grid)  + term_vill_diameter (output from terminal_villous_diameter) *terminals_in_grid
    The denominator comes from tissue_vol (output of tissue_vol_in_samp_gr)
    """

    tissue_diameter = br_diameter_in_grid +term_vill_diameter*terminals_in_grid
    np.seterr(divide='ignore', invalid='ignore')
    weighted_diameter = np.nan_to_num(tissue_diameter/tissue_vol)

    return weighted_diameter


