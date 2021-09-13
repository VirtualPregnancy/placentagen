#!/usr/bin/env python

"""This is the skeleton_to_tree.py module. This module contains code required to generate morphologically accurate tree
structures from skeletons derived from imaging.
"""

import numpy as np

from . import pg_utilities
from . import analyse_tree
from . import imports_and_exports


def assess_graph_structure(pixel_graph, coordinates,dimensions,groupname,outputfile):

    mydegrees = np.zeros(len(coordinates), dtype=int)
    ne = 0
    #closest = np.ones((len(inlets), 2)) * 1000.
    most_allowed_fications = 6
    count_isolated = 0
    count_terminal = 0
    count_bifurcations = 0
    count_morefications = 0
    count_continuing = 0
    node_kount = 0
    elem_kount = 0
    node_map = np.zeros(len(coordinates),dtype=int)
    graph_map = np.zeros(len(coordinates),dtype=int)
    nodal_degrees = np.zeros(len(coordinates))
    nodes = np.zeros((len(coordinates), 4))
    elems = np.zeros((len(coordinates), 3), dtype=int)
    potential_elems = np.zeros((len(coordinates),most_allowed_fications+1),dtype=int)

    for i in range(1, len(pixel_graph.indptr) - 1):
        num_attached = (pixel_graph.indptr[i + 1] - pixel_graph.indptr[i])  # looking at how many attachments it has
        if (num_attached == 0):
            count_isolated = count_isolated + 1
            mydegrees[i] = 0
        elif (num_attached == 1):
            count_terminal = count_terminal + 1
            mydegrees[i] = 1
        elif num_attached == 2:
            count_continuing = count_continuing + 1
            mydegrees[i] = 2
        elif num_attached == 3:
            count_bifurcations = count_bifurcations + 1
            mydegrees[i] = 3
        else:
            count_morefications = count_morefications + 1
            mydegrees[i] = num_attached
        count_new_vessel = 0
        if num_attached>=1 and num_attached<=most_allowed_fications:#possible element
            potential_elems[i,0] = num_attached
            potential_elems[i,1:num_attached+1]=pixel_graph.indices[pixel_graph.indptr[i]:pixel_graph.indptr[i+1]]

    for i in range(1, len(pixel_graph.indptr) - 1):
        for j in range(pixel_graph.indptr[i], pixel_graph.indptr[i + 1]):
            inew = pixel_graph.indices[j] #grab nodes
            #inew is  index of connected branch (old) indexed from one
            np_old = np.where(node_map == i)  # node number
            if(mydegrees[inew]>0):
               if inew not in node_map:
                   currentdegrees = mydegrees[inew]  # how many branches this node is connected to
                   count_new_vessel = count_new_vessel + 1  # create a new vessel segment
                   # Create new node
                   node_map[node_kount] = inew  # Mapping old to new node number
                   graph_map[inew]=node_kount
                   nodes[node_kount, 0] = node_kount  # new node number
                   if dimensions == 3:
                       nodes[node_kount, 1] = coordinates[inew, 0]# coordinates indexed to 'old' node number
                       nodes[node_kount, 2] = coordinates[inew, 1]
                       nodes[node_kount, 3] = coordinates[inew, 2]
                   else: #two dimensionional skeleton images
                       nodes[node_kount, 1] = coordinates[inew, 0]  # coordinates indexed to 'old' node number
                       nodes[node_kount, 2] = coordinates[inew, 1]
                       nodes[node_kount, 3] = 0.0  # dummy z-coord
                   nodal_degrees[node_kount] = currentdegrees
                   node_kount = node_kount + 1  # create a new node



    print('Num isolated points, removed ' + str(count_isolated))
    print('Num terminal nodes (inlets or outlets) ' + str(count_terminal))
    print('Num bifurcations ' + str(count_bifurcations))
    print('Num morefiurcations ' + str(count_morefications))
    print('Mostification ' + str(max(mydegrees)-1))

    for i in range(0,node_kount):
        inew = int(node_map[i])
        if i==0 or int(nodes[i,0]) not in elems[:,1] :
            if nodal_degrees[i] == 1.0:
                #In this case we always more froml
                elems[elem_kount, 0] = elem_kount
                elems[elem_kount, 1] = int(nodes[i, 0])
                elems[elem_kount, 2] = graph_map[potential_elems[inew,1]]
                elem_kount = elem_kount + 1
            else:
                for j in range(0,int(nodal_degrees[i]-1)):
                    connect_me = False
                    for pot in range(1,potential_elems[inew,0]+1):
                        connection = graph_map[potential_elems[inew,pot]]
                        if connection == 0:
                            connect_me = False
                        elif (connection <= int(nodes[i,0])):#moving forward through nodes
                            connect_me = False
                        else:
                            connect_me = True
                            potential_elems[inew, pot] = 0
                            break
                    if connect_me:
                        elems[elem_kount, 0] = elem_kount
                        elems[elem_kount, 1] = int(nodes[i,0])
                        elems[elem_kount, 2] = connection
                        elem_kount = elem_kount+1


    imports_and_exports.export_ex_coords(nodes[:,:][0:node_kount],groupname,outputfile,'exnode')
    imports_and_exports.export_ex_field(nodal_degrees[0:node_kount], groupname, 'degrees', outputfile + '_degrees', 'exnode')
    imports_and_exports.export_exelem_1d(elems[:, :][0:elem_kount], 'arteries', outputfile)