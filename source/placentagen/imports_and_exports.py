#!/usr/bin/env python


import numpy as np
from . import pg_utilities
import warnings
import skimage
from skimage import io


def export_ex_coords(data, groupname, filename, type):
    # Exports coordinates to exnode or exdata format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    # type = exnode or exdata
    print('filename', filename)
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
        if data_length == 4:
            f.write("Node:  "        "%s\n" % int(data[x][0] + 1))
            f.write("          %s\n" % (data[x][1]))
            f.write("          %s\n" % (data[x][2]))
            f.write("          %s\n" % (data[x][3]))
        else:
            f.write("Node:  "        "%s\n" % (x + 1))
            f.write("          %s\n" % data[x][0])
            f.write("          %s\n" % data[x][1])
            f.write("          %s\n" % data[x][2])
    f.close()


def export_ip_coords(data,name,filename):
    # Write header
    type = "ipnode"
    data_num = len(data)
    filename = filename + '.' + type
    f = open(filename, 'w')
    f.write(" CMISS Version 2.1  ipnode File Version 2\n")
    f.write(" Heading: %s\n\n" % name)
    f.write(" The number of nodes is [1]: %s \n" % int(data_num))
    f.write(" Number of coordinates [3]: 3\n")
    f.write(" Do you want prompting for different versions of nj=1 [N]? N\n")
    f.write(" Do you want prompting for different versions of nj=2 [N]? N\n")
    f.write(" Do you want prompting for different versions of nj=3 [N]? N\n")
    f.write(" The number of derivatives for coordinate 1 is [0]: 0\n")
    f.write(" The number of derivatives for coordinate 2 is [0]: 0\n")
    f.write(" The number of derivatives for coordinate 3 is [0]: 0\n")

    # Write element values
    for x in range(0, data_num):
        f.write(" Node number [    1]:     %s\n" % int(x + 1))
        f.write(" The Xj(1) coordinate is [ 0.00000E+00]:  %s\n" % data[x][0])
        f.write(" The Xj(2) coordinate is [ 0.00000E+00]:  %s\n" % data[x][1])
        f.write(" The Xj(3) coordinate is [ 0.00000E+00]:  %s\n\n" % data[x][2])

    f.close()

    return 0

def export_ex_field(data, groupname, fieldname, filename, type):
    # Exports field to exnode or exdata format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    # type = exnode or exdata
    # first entry
    data_num = len(data)
    filename = filename + '.' + type
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" #Fields=1\n")
    f.write(" 1) %s, coordinate, rectangular cartesian, #Components=1\n" % fieldname)
    f.write(" %s.  Value index=1, #Derivatives=0\n" % fieldname)

    for x in range(0, data_num):
        f.write("Node:  "        "%s\n" % (x + 1))
        f.write("          %s\n" % data[x])
    f.close()


def export_nodal_rad_field(data, groupname, fieldname, filename, type, nodes, elems):
    # Exports coordinates to exnode or exdata format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    # type = exnode or exdata
    # first entry
    data_num = len(data)
    filename = filename + '.' + type
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" #Fields=1\n")
    f.write(" 1) %s, coordinate, rectangular cartesian, #Components=1\n" % fieldname)
    f.write(" %s.  Value index=1, #Derivatives=0\n" % fieldname)
    num_per_node = np.zeros(len(nodes))
    node_rad = np.zeros(len(nodes))
    for x in range(0, data_num):
        np1 = elems[x][1]
        np2 = elems[x][2]
        num_per_node[np1] = num_per_node[np1] + 1.
        num_per_node[np2] = num_per_node[np2] + 1.
        node_rad[np1] = node_rad[np1] + data[x]
        node_rad[np2] = node_rad[np2] + data[x]

    for y in range(0, len(nodes)):
        node_rad[y] = node_rad[y] / num_per_node[y]
        f.write("Node:  "        "%s\n" % (y + 1))
        f.write("          %s\n" % (node_rad[y]))
    f.close()

def export_ipfiel(data,filename):
    data_num = len(data)
    # write the node radius to an ipfiel file
    filename = filename + '.ipfiel'
    f = open(filename, 'w')

    f.write("CMISS Version 2.1  ipelem File Version 2\n")
    f.write("Heading:\n")
    f.write("\n")
    f.write("The number of nodes is [     {0}]:      {1}\n".format(data_num, data_num))
    f.write("Do you want prompting for different versions of field variable 1 [N]? Y\n")
    f.write("The number of derivatives for field variable 1 is [0]: 0\n")

    for n in range(0, data_num):
        f.write("\n")
        f.write("Node number [     {0}]:    {1}\n".format(n + 1, n + 1))
        f.write("The number of versions for field variable 1 is [1]:  1\n")
        # use scientific notation for radius values
        f.write("The field variable 1 value is [ {0:.5e}]:  {1:.5e}\n".format(data[n], data[n]))

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



def export_ipelem_1d(data, name, filename):
    # Write header
    type = "ipelem"
    data_num = len(data)
    filename = filename + '.' + type
    f = open(filename, 'w')
    f.write(" CMISS Version 2.1  ipelem File Version 2\n")
    f.write(" Heading: %s\n\n" % name)
    f.write(" The number of elements is [1]: %s \n\n" % int(data_num))

    # Write element values
    for x in range(0, data_num):
        f.write(" Element number [    1]:     %s\n" % int(x + 1))
        f.write(" The number of geometric Xj-coordinates is [3]: 3\n")
        f.write(" The basis function type for geometric variable 1 is [1]:  1\n")
        f.write(" The basis function type for geometric variable 2 is [1]:  1\n")
        f.write(" The basis function type for geometric variable 3 is [1]:  1\n")
        f.write(" Enter the 2 global numbers for basis 1: %s %s\n\n" % (int(data[x][1] + 1), int(data[x][2] + 1)))

    f.close()

    return 0


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


def export_exelem_3d_linear_list(data, list, groupname, filename):
    # Exports element locations to exelem format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    data_num = len(list)
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
        y = list[x]
        f.write(" Element:            %s 0 0\n" % int(data[x][0] + 1))
        f.write("   Nodes:")
        f.write(
            "                %s            %s            %s            %s            %s            %s            %s            %s\n" % (
                int(data[y][1] + 1), int(data[y][2] + 1), int(data[y][3] + 1), int(data[y][4] + 1),
                int(data[y][5] + 1),
                int(data[y][6] + 1), int(data[y][7] + 1), int(data[y][8] + 1)))

    f.close()


def export_exfield_3d_linear(data, groupname, fieldname, filename):
    # Exports element fields to exelem format
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


def export_exfield_3d_linear_list(data, list, groupname, fieldname, filename):
    # Exports element fields to exelem format when data is defined at a specified list of nodes
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    data_num = len(list)
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
        exp_data = data[list[x]]
        f.write(" Element:            %s 0 0\n" % int(x + 1))
        f.write("   Values:\n")
        f.write(
            "           %s       %s       %s       %s       %s       %s       %s       %s\n" % (
                exp_data, exp_data, exp_data, exp_data, exp_data, exp_data, exp_data, exp_data))

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


######
# Function: takes data from the csv and converts it to arrays
# Inputs: data_file - generated from the panadas read_csv function, containing results from imageJ image analysis
#         Arrays - a group of arrays each with length N for their first axis
# Outputs: nodes - an M x 3 array giving cartesian coordinates (x,y,z) for the node locations in the tree
#         elems - an N x 3 array, the first colum in the element number, the second two columns are the index of the start and end node
#         radii, length, euclidean_length - there are all an N x 1 array containing a property for each element
######

def import_imagej_skel_csv(data_file, keep_skeleton, what_skel):
    # If keep_skeleton = 0 we keep all the elements in the skeleton
    # Otherwise we can select one out from the dataset to analyse
    if what_skel == "less":  # implies we want <= value
        data_file = data_file[data_file.SkeletonID <= keep_skeleton]
    elif what_skel == "single":
        data_file = data_file[data_file.SkeletonID == keep_skeleton]
    elif what_skel == "all":
        print("reading all skeletons: could take a while")
    else:
        print("Not a valid option for reading skeltons")
        return

    # get skeleton properties as arrays
    euclid_length = data_file.Euclideandistance.values
    length = data_file.Branchlength.values
    radii = data_file.averageintensityinner3rd.values
    branch_id = data_file.SkeletonID.values

    print("sorting data")
    # get elem and node data
    data_file = data_file.drop(['SkeletonID', 'Branchlength', 'averageintensityinner3rd', 'Euclideandistance'], axis=1)
    data_file = data_file.values
    (elems, nodes) = pg_utilities.sort_elements(data_file[:, 0:3], data_file[:, 3:6])

    print("elements sorted")
    # get rid of dud elements
    (elems, [length, euclid_length, radii, branch_id]) = pg_utilities.remove_rows(elems, [length, euclid_length,
                                                                                          radii, branch_id])

    return {'nodes': nodes, 'elems': elems, 'radii': radii, 'length': length, 'euclidean length': euclid_length,
            'branch_id': branch_id}


def import_stemxy(stem_file):
    # reading in the stem vessel to map the spiral artery location
    stem_xy = open(stem_file, 'r')
    stem_coor = stem_xy.readlines()  # readlines
    startLines = range(0, len(stem_coor))

    for i in range(len(stem_coor)):
        stem_coor[i] = stem_coor[i].split()
    stem_xyList = []
    stem_elemList = []
    for i in startLines:
        node = []
        node.append(float(stem_coor[i][0]))  # x coor of stem villi
        node.append((float(stem_coor[i][1])))  # y coor of stem villi
        stem_xyList.append(node)
        elem = int(stem_coor[i][2]) - 1
        stem_elemList.append(elem)
    stem_xy.close()

    return {'stem_xy': stem_xyList, 'elem': stem_elemList}


def import_exnode_tree(filename):
    # count nodes for check of correct number for the user, plus use in future arrays
    count_node = 0
    # Initialise array of node numbers and values
    node_array = np.empty((0, 7))
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
                node_array = np.append(node_array, np.zeros((1, 7)),
                                       axis=0)  # initialise a list of attributes for each node
                node_array[count_node - 1][count_atribute] = int(str.split(line)[1]) - 1
            else:
                line_num = is_float(line_type)  # checking if the line is a number
                if (line_num):  # it is a number
                    if not "index" in line:
                        count_atribute = count_atribute + 1
                        node_array[count_node - 1][count_atribute] = float(str.split(line)[0])

    if (count_atribute < 7):
        node_array = np.delete(node_array, np.s_[count_atribute + 1:7], axis=1)
    total_nodes = count_node
    return {'total_nodes': total_nodes, 'nodes': node_array}


def import_exelem_tree(filename):
    # count element for check of correct number for the user, plus use in future arrays
    count_el = 0
    # Initialise array of el numbers and values
    el_array = np.empty((0, 3), dtype=int)
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
                el_array = np.append(el_array, np.zeros((1, 3), dtype=int), axis=0)
                el_array[count_el - 1][count_atribute] = int(str.split(line)[1]) - 1
            else:
                line_num = is_float(line_type)  # checking if the line is a number
                if (line_num):  # it is a number
                    if "#Values" not in line and "l.Lagrange" not in line and "1.000000000000000e+00" not in line and "0.1000000000000000E+01" not in line and '1.000000e+00' not in line:
                        count_atribute = count_atribute + 1
                        el_array[count_el - 1][count_atribute] = float(str.split(line)[0]) - 1  # first node of element
                        el_array[count_el - 1][count_atribute + 1] = float(
                            str.split(line)[1]) - 1  # 2nd node of element

    total_el = count_el
    return {'total_elems': total_el, 'elems': el_array}
    
def import_exelem_field(filename):
    # count element for check of correct number for the user, plus use in future arrays
    count_el = 0
    # Initialise array of el numbers and values
    el_array = np.empty((0))
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
                el_array = np.append(el_array, np.zeros((1)), axis=0)
                #el_array[count_el - 1][count_atribute] = int(str.split(line)[1]) - 1
            else:
                line_num = is_float(line_type)  # checking if the line is a number
                if (line_num):  # it is a number
                    if "#Values" not in line and "l.Lagrange" not in line:
                        #count_atribute = count_atribute + 1
                        el_array[count_el - 1] = float(str.split(line)[0])   # first node of element
                        #el_array[count_el - 1][count_atribute + 1] = float(
                        #    str.split(line)[1]) - 1  # 2nd node of element

    total_el = count_el
    return el_array


def is_float(str):
    try:
        num = float(str)
    except ValueError:
        return False
    return True


def export_exelem_3d_quadratic(data, groupname, filename):
    # Exports element locations to exelem format
    # data = array of data
    # groupname = what you want your data to be called in cmgui
    # filename = file name without extension
    data_num = len(data)
    filename = filename + '.exelem'
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write(" Shape.  Dimension=3\n")
    f.write(" #Scale factor sets= 1\n")
    f.write(" q.Lagrange*q.Lagrange*q.Lagrange, #Scale factors=27\n")
    f.write(" #Nodes=           27\n")
    f.write(" #Fields=1\n")
    f.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
    f.write("   x.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 27\n")
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
    f.write("      9.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      10.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      11.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      12.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      13.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      14.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      15.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      16.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      17.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      18.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      19.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      20.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      21.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      22.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      23.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      24.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      25.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      26.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      27.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("   y.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 27\n")
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
    f.write("      9.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      10.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      11.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      12.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      13.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      14.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      15.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      16.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      17.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      18.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      19.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      20.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      21.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      22.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      23.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      24.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      25.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      26.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      27.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("   z.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.\n")
    f.write("     #Nodes= 27\n")
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
    f.write("      9.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      10.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      11.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      12.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      13.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      14.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      15.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      16.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      17.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      18.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      19.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      20.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      21.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      22.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      23.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      24.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      25.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      26.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    f.write("      27.  #Values=1\n")
    f.write("       Value indices:     1\n")
    f.write("       Scale factor indices:   0\n")
    for x in range(0, data_num):
        f.write(" Element:            %s 0 0\n" % int(data[x][0] + 1))
        f.write("   Nodes:")
        f.write(
            "                %s            %s            %s            %s            %s            %s            %s            %s                %s            %s            %s            %s            %s            %s            %s            %s                 %s            %s            %s            %s            %s            %s            %s            %s                %s            %s            %s            \n" % (
                int(data[x][1] + 1), int(data[x][2] + 1), int(data[x][3] + 1), int(data[x][4] + 1), int(data[x][5] + 1),
                int(data[x][6] + 1), int(data[x][7] + 1), int(data[x][8] + 1), int(data[x][9] + 1),
                int(data[x][10] + 1), int(data[x][11] + 1), int(data[x][12] + 1), int(data[x][13] + 1),
                int(data[x][14] + 1), int(data[x][15] + 1), int(data[x][16] + 1), int(data[x][17] + 1),
                int(data[x][18] + 1), int(data[x][19] + 1), int(data[x][20] + 1), int(data[x][21] + 1),
                int(data[x][22] + 1), int(data[x][23] + 1), int(data[x][24] + 1), int(data[x][25] + 1),
                int(data[x][26] + 1), int(data[x][27] + 1)))
        f.write("Scale factors:\n")
        f.write(
            "1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00   1.0000000000000000E+00\n")

    f.close()


def export_exfield_3d_quadratic(data, groupname, fieldname, filename):
    # Exports element fields to exelem format
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
            "      %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s    %s\n" % (
                data[x], data[x], data[x], data[x], data[x], data[x], data[x], data[x], data[x], data[x],
                data[x], data[x], data[x], data[x], data[x], data[x], data[x], data[x], data[x], data[x],
                data[x], data[x], data[x], data[x], data[x], data[x], data[x]))

    f.close()

    ######
    # Function: Loads in a stack of images, located in path, and with naming convention name (goes slice at a time to avoid memory errors)
    #     Inputs: numImages - integer, number of images in the stack
    #             name - string for name of images. Note images must be numbered from 0
    #     Outputs: Image, a 3D BOOLEAN array containing image
    ######


def load_image_bool(name, numImages):
    # read in first image + get dimensions to initialize array
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        print(name.format(0))
        im = io.imread(name.format(0))
        gray_image = skimage.color.rgb2gray(im)
        skimage.img_as_bool(gray_image)
        Image = np.zeros([im.shape[0], im.shape[1], numImages], dtype=bool)
        Image[:, :, 0] = gray_image

        # load all slices
    for i in range(0, numImages):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            im = io.imread(name.format(i))
            gray_image = skimage.color.rgb2gray(im)
            skimage.img_as_bool(gray_image)
        Image[:, :, i] = gray_image

    print('Image ' + name + ' loaded. Shape: ' + str(Image.shape))
    return Image
