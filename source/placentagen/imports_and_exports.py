#!/usr/bin/env python
import numpy as np

def export_ex_coords(data, groupname, filename, type):
    # Exports coordinates to exnode or exdata format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    # type = exnode or exdata
    data_length = len(
        data[0])  # if this is 3 then number nodes or data automatically if 4 then node numbers are given as
    # first entry
    data_num = len(data)
    filename = filename + '.' + type
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" #Fields=1\n")
    f.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
    f.write(" x.  Value index=1, #Derivatives=0\n")
    f.write(" y.  Value index=1, #Derivatives=0\n")
    f.write(" z.  Value index=1, #Derivatives=0\n")

    for x in range(0, data_num):
        if data_length is 4:
            f.write("Node:  "        "%s\n" % int(data[x][0] + 1))
            f.write("          %s\n" % data[x][1])
            f.write("          %s\n" % data[x][2])
            f.write("          %s\n" % data[x][3])
        else:
            f.write("Node:  "        "%s\n" % (x + 1))
            f.write("          %s\n" % data[x][0])
            f.write("          %s\n" % data[x][1])
            f.write("          %s\n" % data[x][2])
    f.close()


def export_exelem_1d(data, groupname, filename):
    # Exports element locations to exelem format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    data_num = len(data)
    filename = filename + '.exelem'
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" Shape.  Dimension=1\n")
    f.write(" #Scale factor sets= 1\n")
    f.write("   l.Lagrange, #Scale factors= 2\n")
    f.write(" #Nodes=           2\n")
    f.write(" #Fields=1\n")
    f.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
    f.write("   x.  l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 2\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   1\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   2\n")
    f.write("   y.  l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 2\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   1\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   2\n")
    f.write("   z.  l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 2\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   1\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   2\n")
    for x in range(0, data_num):
        f.write(" Element:            %s 0 0\n" % int(data[x][0] + 1))
        f.write("   Nodes:\n")
        f.write("                %s            %s\n" % (int(data[x][1] + 1), int(data[x][2] + 1)))
        f.write("   Scale factors:\n")
        f.write("       0.1000000000000000E+01   0.1000000000000000E+01\n")
    f.close()


def export_exelem_3d_linear(data, groupname, filename):
    # Exports element locations to exelem format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    data_num = len(data)
    filename = filename + '.exelem'
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" Shape. Dimension=3 line*line*line\n")
    f.write(" #Scale factor sets= 0\n")
    f.write(" #Nodes=           8\n")
    f.write(" #Fields=1\n")
    f.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
    f.write("   x.  l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 8\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      3.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      4.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      5.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      6.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      7.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      8.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("   y.  l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 8\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      3.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      4.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      5.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      6.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      7.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      8.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("   z.  l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 8\n")
    f.write("      1.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      2.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      3.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      4.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      5.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      6.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      7.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      8.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    for x in range(0, data_num):
        f.write(" Element:            %s 0 0\n" % int(data[x][0] + 1))
        f.write("   Nodes:")
        f.write(
            "                %s            %s            %s            %s            %s            %s            %s            %s\n" % (
                int(data[x][1] + 1), int(data[x][2] + 1), int(data[x][3] + 1), int(data[x][4] + 1), int(data[x][5] + 1),
                int(data[x][6] + 1), int(data[x][7] + 1), int(data[x][8] + 1)))

    f.close()

def export_exfield_3d_linear(data, groupname, fieldname, filename):
    # Exports element locations to exelem format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    data_num = len(data)
    filename = filename + '.exelem'
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" Shape. Dimension=3 line*line*line\n")
    f.write(" #Scale factor sets= 0\n")
    f.write(" #Nodes=           0\n")
    f.write(" #Fields=1\n")
    f.write(" 1) %s, field, rectangular cartesian, #Components=1\n" % fieldname)
    f.write("   %s.  l.Lagrange*l.Lagrange*l.Lagrange, no modify, grid based.\n" % fieldname)
    f.write("   #xi1=1 \n")
    f.write("   #xi2=1 \n")
    f.write("   #xi3=1 \n")
    for x in range(0, data_num):
        f.write(" Element:            %s 0 0\n" % int(x + 1))
        f.write("   Values:\n")
        f.write(
            "           %s       %s       %s       %s       %s       %s       %s       %s\n" % (
                data[x], data[x], data[x], data[x], data[x], data[x], data[x], data[x]))

    f.close()

def export_exfield_1d_linear(data, groupname, fieldname, filename):
    # Exports element locations to exelem format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    data_num = len(data)
    filename = filename + '.exelem'
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" Shape.  Dimension=1\n")
    f.write(" #Scale factor sets= 0\n")
    f.write(" #Nodes=           0\n")
    f.write(" #Fields=1\n")
    f.write(" 1) %s, field, rectangular cartesian, #Components=1\n" % fieldname)
    f.write("   %s.  l.Lagrange, no modify, grid based.\n" % fieldname)
    f.write("   #xi1=1 \n")
    for x in range(0, data_num):
        f.write(" Element:            %s 0 0\n" % int(x + 1))
        f.write("   Values:\n")
        f.write(
            "           %s       %s\n" % (
                data[x], data[x]))
    f.close()


def import_exnode_tree(filename):
    # count nodes for check of correct number for the user, plus use in future arrays
    count_node = 0
    # Initialise array of node numbers and values
    node_array = np.empty((0,7))
    # open file
    with open(filename) as f:
        # loop through lines of file
        while 1:
            line = f.readline()
            if not line:
                break  # exit if done with all lines
            # identifying whether there is a node defined here
            line_type = str.split(line)[0]
            if (line_type == 'Node:'):  # line defines new node
                count_node = count_node + 1  # count the node
                count_atribute = 0  # intitalise attributes of the node (coordinates, radius)
                node_array=np.append(node_array,np.zeros((1,7)),axis = 0)  # initialise a list of attributes for each node
                node_array[count_node - 1][count_atribute] = int(str.split(line)[1])-1
            else:
                line_num = is_float(line_type)  # checking if the line is a number
                if (line_num):  # it is a number
                    if not "index" in line:
                        count_atribute = count_atribute + 1
                        node_array[count_node - 1][count_atribute] = float(str.split(line)[0])
    total_nodes = count_node
    return {'total_nodes': total_nodes, 'nodes': node_array}


def import_exelem_tree(filename):
    # count element for check of correct number for the user, plus use in future arrays
    count_el = 0
    # Initialise array of el numbers and values
    el_array = np.empty((0,3),dtype = int)
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
                el_array = np.append(el_array, np.zeros((1, 3),dtype = int), axis=0)
                el_array[count_el - 1][count_atribute] = int(str.split(line)[1])-1
            else:
                line_num = is_float(line_type)  # checking if the line is a number
                if (line_num):  # it is a number
                    if "#Values" not in line and "l.Lagrange" not in line and "0.1000000000000000E+01" not in line:
                        count_atribute = count_atribute + 1
                        el_array[count_el - 1][count_atribute] = float(str.split(line)[0])-1  # first node of element
                        el_array[count_el - 1][count_atribute + 1] = float(str.split(line)[1])-1  # 2nd node of element

    total_el = count_el
    return {'total_elems': total_el, 'elems': el_array}


def is_float(str):
    try:
        num = float(str)
    except ValueError:
        return False
    return True
