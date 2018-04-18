#!/usr/bin/env python
import numpy as np
from . import pg_utilities
import time


def calc_terminal_branch(node_loc, elems):
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


    print('length of data', len(placenta_list))

    for ne_i in range(0, num_sample_elems):
        # First node has min x,y,z and last node has max x,y,z
        ne = placenta_list[ne_i]
        start = time.time()
        if placenta_list[ne_i] > 0: #There is some placenta in this element (assuming none in el 0)
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
                        print(nt, in_element)
                        terminals_in_grid[ne] = terminals_in_grid[ne] + 1
                        terminal_mapped[nt] = 1
        end = time.time()
        print('Time taken for this element: ' + str(ne))
        print(end-start, np.sum(terminal_mapped))
    return terminals_in_grid


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

    #Initialise the array that defines the volume of placenta in each grid element
    pl_vol_in_grid = np.zeros(total_elems)
    non_empty_loc = np.zeros(total_elems, dtype = int)
    non_empty_count = 0

    for ne in range(0, len(elems)):  # looping through elements
        count_in_range = 0
        nod_in_range = np.zeros(8, dtype=int)
        #define range of x, y , and z in the element
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
            pl_vol_in_grid[ne] = (endx - startx) * (endy - starty) * (endz - startz)  # the placental vol in that samp_grid_el is same as vol of samp_grid_el
        elif count_in_range == 0:  # if all 8 nodes are outside the ellpsiod
            # since this samp_grid_el is completely outside, the placental vol is zero
            pl_vol_in_grid[ne] = 0
        else:  # if some nodes in and some nodes out, the samp_grid_el is at the edge of ellipsoid
        #Use trapezoidal quadrature to caculate the volume under the surface of the ellipsoid in each element
            non_empty_loc[non_empty_count] = ne
            non_empty_count = non_empty_count + 1
            # need to map to positive quadrant
            repeat = False
            if (startz < 0 and endz <=0):
                #need to project to positive z axis
                startz = abs(nodes[elems[ne][8]][2])
                endz = abs(nodes[elems[ne][1]][2])
            elif(startz < 0 and endz > 0):
                #Need to split into components above and below the axis and sum the two
                startz = 0
                endz = abs(nodes[elems[ne][1]][2])
                startz_2=0
                endz_2 =nodes[elems[ne][8]][2]
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
                pl_vol_in_grid[ne] =pl_vol_in_grid[ne] + (Value1 - startz_2 * (endx - startx) * (
                                endy - starty))

    print('Number of Non-empty cells: ' + str(non_empty_count))
    print('Total number of cells: ' + str(total_elems))
    non_empty_loc = np.resize(non_empty_loc, non_empty_count)

    return {'pl_vol_in_grid': pl_vol_in_grid, 'non_empty_rects': non_empty_loc}
