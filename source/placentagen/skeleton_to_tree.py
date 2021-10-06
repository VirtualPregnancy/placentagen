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
    kount_duplicates = 0
    for i in range(0,node_kount):
        inew = int(node_map[i]) #Original identifier for this node
        if i==0 or int(nodes[i,0]) not in elems[:,1]: #only adding a branch if we haven't seen this node before
            if nodal_degrees[i] == 1.0: #This node only has one branch coming from it
                #In this case we always include a branch from the current to the next node
                elems[elem_kount, 0] = elem_kount
                elems[elem_kount, 1] = int(nodes[i, 0])
                elems[elem_kount, 2] = graph_map[potential_elems[inew,1]]
                elem_kount = elem_kount + 1
            else:
                for j in range(0,int(nodal_degrees[i]-1)): #looping through all defined connections to this branch
                    connect_me = False
                    for pot in range(1,potential_elems[inew,0]+1):
                        connection = graph_map[potential_elems[inew,pot]] #node to connect to
                        if connection == 0: #zeroed if already dealt with
                            connect_me = False
                            #elif (connection == int(nodes[i,0])):#don't connect a node to itself
                            #potential_elems[inew, pot] = 0 #zeroing this potential connection as we have created an element for it
                            #connect_me = False
                        else:
                            connect_me = True
                            potential_elems[inew, pot] = 0 #zeroing this potential connection as we have created an element for it
                            break #break this logical, we'll add a branch
                    if connect_me: #There is a branch to connect
                        new_branch = [min(int(nodes[i,0]),connection), max(int(nodes[i,0]),connection)]
                        duplicate = False
                        for k in range(0,elem_kount):
                            if all(new_branch == elems[k,1:3]): #excluding duplicates
                                kount_duplicates = kount_duplicates + 1
                                print('Duplicate',new_branch,kount_duplicates)
                                duplicate = True
                        if not duplicate:
                            elems[elem_kount, 0] = elem_kount
                            elems[elem_kount, 1] = new_branch[0]#Order so that we have smallest node number first
                            elems[elem_kount, 2] = new_branch[1]
                            elem_kount = elem_kount+1

    ('Total number of elems',elem_kount)
    print('Total number of nodes', node_kount)
    print('Possible duplicates', kount_duplicates)
    imports_and_exports.export_ex_coords(nodes[:,:][0:node_kount],groupname,outputfile,'exnode')
    imports_and_exports.export_ex_field(nodal_degrees[0:node_kount], groupname, 'degrees', outputfile + '_degrees', 'exnode')
    imports_and_exports.export_exelem_1d(elems[:, :][0:elem_kount], 'arteries', outputfile)

def find_distances_using_normal(coord1, coord2, VolumeImage):

    # get centre line vector
    centre = np.double((coord1 - coord2))/ np.linalg.norm(np.double(coord1 - coord2))

    numSamples = 20
    distances = np.zeros(numSamples)
    normal = np.zeros(3)

    for i in range(0,numSamples):

        # Randomly assign normal vector, using the dot product rule (centre.normal==0) and avoiding div0 errors
        if centre[2]!= 0:
            normal[0] = np.random.rand() - 0.5
            normal[1] = np.random.rand() - 0.5
            normal[2] = -(centre[0] * normal[0] + centre[1] * normal[1])/ centre[2]

        elif centre[0] != 0:
            normal[2] = np.random.rand() - 0.5
            normal[1] = np.random.rand() - 0.5
            normal[0] = -(centre[2] * normal[2] + centre[1] * normal[1]) / centre[0]

        else: # centre[1]!= 0:
            normal[0] = np.random.rand() - 0.5
            normal[2] = np.random.rand() - 0.5
            normal[1] = -(centre[0] * normal[0] + centre[2] * normal[2]) / centre[1]

        normal = normal / np.linalg.norm(normal)

        # Find distances
        step = 0
        counter = 0
        currentValue = 1
        while (currentValue == 1) & (counter < 1000): # check if in vessel (plus arbitrary check)

             step = step + 0.1 # step update by 1/5 of a voxel (could increase in order to speed up)
             counter = counter + 1
             currentPosition = np.double(coord1) + step*normal # take step in direction of normal vector
             currentPosition = np.round(currentPosition)
             if(int(currentPosition[0])>np.shape(VolumeImage)[0]-1):
                 currentValue = 0
             elif(int(currentPosition[1])>np.shape(VolumeImage)[1]-1):
                 currentValue = 0
             elif(int(currentPosition[2])>np.shape(VolumeImage)[2]-1):
                 currentValue = 0
             else:
                 currentValue = VolumeImage[int(currentPosition[0]), int(currentPosition[1]), int(currentPosition[2])]
        distances[i] = step - 0.1

    return distances