#!/usr/bin/env python
import numpy as np
from . import pg_utilities

def calc_terminal_branch(node_loc,elems):
    #This function generates a list of terminal nodes associated with a branching geometry
    #inputs are node locations and elements
    num_elems = len(elems)
    num_nodes = len(node_loc)
    elem_cnct = pg_utilities.element_connectivity_1D(node_loc, elems)

    terminal_branches = np.zeros(num_elems, dtype = int)
    terminal_nodes = np.zeros(num_nodes, dtype = int)

    num_term = 0
    for ne in range(0,num_elems):
        if elem_cnct['elem_down'][ne][0] == 0: #no downstream element
            terminal_branches[num_term] = ne
            terminal_nodes[num_term] = elems[ne][2] #node that leaves the terminal element
            num_term = num_term + 1

    terminal_branches = np.resize(terminal_branches,num_term)
    terminal_nodes = np.resize(terminal_nodes,num_term)

    print('Total number of terminals assessed, num_terminals =  ' + str(num_term))

    return {'terminal_elems': terminal_branches, 'terminal_nodes': terminal_nodes, 'total_terminals': num_term}


def terminals_in_sampling_grid(rectangular_mesh,terminal_list,node_loc):
    #This function counts the number of terminals in a sampling grid element
    #inputs are:
    #Rectangular mesh - the sampling grid
    #terminal_list - a list of terminals
    #node_loc - location of nodes
    num_sample_elems = rectangular_mesh['total_elems']
    num_terminals = terminal_list['total_terminals']
    terminals_in_grid = np.zeros(num_sample_elems,dtype = int)

    for ne in range(0,num_sample_elems):
        #First node has min x,y,z and last node has max x,y,z
        first_node = rectangular_mesh['elems'][ne][1]
        last_node = rectangular_mesh['elems'][ne][8]
        min_coords = rectangular_mesh['nodes'][first_node][0:3]
        max_coords = rectangular_mesh['nodes'][last_node][0:3]
        for nt in range(0,num_terminals):
            in_element = [False, False, False]
            coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
            if coord_terminal[0] >= min_coords[0] and coord_terminal[0] < max_coords[0]:
                in_element[0] = True
                if coord_terminal[1] >= min_coords[1] and coord_terminal[1] < max_coords[1]:
                    in_element[1] = True
                    if coord_terminal[2] >= min_coords[2] and coord_terminal[2] < max_coords[2]:
                        in_element[2] = True
            if(np.all(in_element)):
                terminals_in_grid[ne] = terminals_in_grid[ne]+1

    return terminals_in_grid

def ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity,num_test_points):
    # This subroutine calculates the placental volume associated with each element in a samplling grid
    # inputs are:
    # rectangular_mesh = the sampling grid nodes and elements
    # volume = placental volume
    # thickness = placental thickness
    # ellipiticity = placental ellipticity
    # x_spacing,y_spacing,z_spaing (dont think these are needed)
    total_elems = rectangular_mesh['total_elems']
    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']
    x = nodes[:, 0]
    y = nodes[:, 1]
    z = nodes[:, 2]

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    pl_vol_in_grid = np.zeros(
        total_elems)  # 1st col stores the number of samp_grid_el, 2nd col stores 0/1/2 (depend on type of samp_grid_el), 3rd col stores the pl_vol in that sam_grid_el


    xVector = np.linspace(0.0, 1.0, num_test_points)
    yVector = np.linspace(0.0, 1.0, num_test_points)
    zVector = np.linspace(0.0, 1.0, num_test_points)

    calculating_nodes = np.vstack(np.meshgrid(xVector, yVector, zVector)).reshape(3, -1).T

    print(calculating_nodes)

    for ne in range(0, len(elems)):  # looping through elements
        count_in_range = 0
        nod_in_range = np.zeros(8, dtype=int)
        for nod in range(1, 9):
            check_in_range = pg_utilities.check_in_ellipsoid(nodes[elems[ne][nod]][0], nodes[elems[ne][nod]][1],
                                                             nodes[elems[ne][nod]][2], x_radius, y_radius, z_radius)
            check_on_range = pg_utilities.check_on_ellipsoid(nodes[elems[ne][nod]][0], nodes[elems[ne][nod]][1],
                                                             nodes[elems[ne][nod]][2], x_radius, y_radius, z_radius)
            if check_in_range or check_on_range:
                count_in_range = count_in_range + 1
                nod_in_range[nod - 1] = 1
        if count_in_range == 8:  # if all 8 nodes are inside the ellipsoid
            pl_vol_in_grid[ne] = 1.0  # the placental vol in that samp_grid_el is same as vol of samp_grid_el
        if count_in_range == 0:  # if all 8 nodes are outside the ellpsiod
            pl_vol_in_grid[ne] = 0  # since this samp_grid_el is completely outside, the placental vol is zero (there will be no tree here and no need to worried about the vol)
        else:  # if some nodes in and some nodes out, the samp_grid_el is at the edge of ellipsoid
            startx = nodes[elems[ne][1]][0]
            endx = nodes[elems[ne][8]][0]
            starty = nodes[elems[ne][1]][1]
            endy = nodes[elems[ne][8]][1]
            startz = nodes[elems[ne][1]][2]
            endz = nodes[elems[ne][8]][2]




    print(np.sum(pl_vol_in_grid))
    print(pl_vol_in_grid)

    return pl_vol_in_grid

