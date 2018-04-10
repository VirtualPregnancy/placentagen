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

