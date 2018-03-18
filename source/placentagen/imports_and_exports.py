#!/usr/bin/env python
import os


def import_exnode_tree(filename):
    print(filename)
    print(os.getcwd())
    # count nodes for check of correct number for the user, plus use in future arrays
    count_node = 0
    # Initialise array of node numbers and values
    node_array = []
    # open file
    with open(filename) as f:
        # loop through lines of file
        while 1:
            line = f.readline()
            if not line:
                break  # exit if done with all lines
            # identifying whether there is a node defined here
            line_type = str.split(line)[0]
            if (line_type == 'Node:'):  # line dedfines new node
                count_node = count_node + 1  # count the node
                count_atribute = 0  # intitalise attributes of the node (coordinates, radius)
                node_array.append([0, 0, 0, 0, 0, 0, 0])  # initialise a list of attributes for each node
                print(int(str.split(line)[1]))
                node_array[count_node - 1][count_atribute] = int(str.split(line)[1])
            else:
                line_num = is_float(line_type)  # checking if the line is a number
                if (line_num):  # it is a number
                    if not "index" in line:
                        count_atribute = count_atribute + 1
                        node_array[count_node - 1][count_atribute] = float(str.split(line)[0])

    total_nodes = count_node
    return {'total_nodes': total_nodes, 'node_array': node_array}



def import_exelem_tree(filename):
    print(filename)
    print(os.getcwd())

    # count element for check of correct number for the user, plus use in future arrays
    count_el = 0
    # Initialise array of el numbers and values
    el_array = []
    # open file
    with open(filename) as f:
        # loop through lines of file
        while 1:
            line = f.readline()
            if not line:
                break  # exit if done with all lines
            # identifying whether there is an element defined here
            line_type = str.split(line)[0]
           
            if (line_type == 'Element:'):  # line dedfines new el
                count_el = count_el + 1  # count the el
                count_atribute = 0  # intitalise attributes of the el (1st el, 2nd el)
                el_array.append([0, 0, 0, 0, 0, 0, 0])  # initialise a list of attributes for each el
                #print(int(str.split(line)[1]))
                el_array[count_el - 1][count_atribute] = int(str.split(line)[1])
                
            else:
                line_num = is_float(line_type)  # checking if the line is a number
                if (line_num):  # it is a number
                    if "#Values" not in line and "l.Lagrange" not in line and "0.1000000000000000E+01" not in line:
                         count_atribute = count_atribute + 1
                         el_array[count_el - 1][count_atribute]   = float(str.split(line)[0]) # first node of element
                         el_array[count_el - 1][count_atribute+1] = float(str.split(line)[1]) # 2nd node of element
                      

    total_el = count_el
    return {'total_el': total_el, 'el_array': el_array}

def is_float(str):
    try:
        num = float(str)
    except ValueError:
        return False
    return True
