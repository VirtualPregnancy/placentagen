#!/usr/bin/env python

"""This is the skeleton_to_tree.py module. This module contains code required to generate morphologically accurate tree
structures from skeletons derived from imaging.
"""

import numpy as np

from . import pg_utilities
from . import analyse_tree
from . import imports_and_exports


def create_graph_structure(pixel_graph, coordinates,dimensions,groupname,outputfile):
    mydegrees = np.zeros(len(coordinates), dtype=int)
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
    potential_elems = np.empty((0, 2), dtype=int)#np.zeros((len(coordinates)*2,2),dtype=int)
    print(potential_elems.shape)
    count_new_vessel = 0
    for i in range(1, len(pixel_graph.indptr) - 1): #Looping through all points

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

        if num_attached>=1 and num_attached<=most_allowed_fications:#possible element
            print(pixel_graph.indices[pixel_graph.indptr[i]:pixel_graph.indptr[i+1]],num_attached)
            attached_elems = pixel_graph.indices[pixel_graph.indptr[i]:pixel_graph.indptr[i+1]]
            for k in range(0,num_attached):
                print('ink',k, num_attached,attached_elems[k], attached_elems[k],count_new_vessel)
                if attached_elems[k] != i: #avoiding self connections (simple loops)
                    potential_elems = np.append(potential_elems, np.zeros((1, 2), dtype=int), axis=0)
                    potential_elems[count_new_vessel,0] = int(min(i,attached_elems[k]))
                    potential_elems[count_new_vessel,1]= int(max(i,attached_elems[k]))
                    count_new_vessel = count_new_vessel + 1
        #    print('Excluding',pixel_graph.indptr[i],pixel_graph.indptr[i+1])

    print(potential_elems.shape)

    potential_elems = np.unique(potential_elems, axis=0)

    print(potential_elems.shape)

    print("Creating nodes")
    for i in range(1, len(pixel_graph.indptr) - 1):
        for j in range(pixel_graph.indptr[i], pixel_graph.indptr[i + 1]):
            inew = pixel_graph.indices[j] #grab nodes
            #inew is  index of connected branch (old) indexed from one
            np_old = np.where(node_map == i)  # node number
            if(mydegrees[inew]>0):
               if inew not in node_map:
                   currentdegrees = mydegrees[inew]  # how many branches this node is connected to
                   #count_new_vessel = count_new_vessel + 1  # create a new vessel segment
                   # Create new node
                   node_map[node_kount] = inew  # Mapping old to new node number
                   graph_map[inew]=node_kount #maps old to new node number
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
    print("Number of Unique elements = ", potential_elems.shape[0])
    print('Total number of nodes', node_kount)
    elems =  np.zeros((potential_elems.shape[0], 3), dtype=int)
    for i in range(0,potential_elems.shape[0]):
        elems[i,0] = i
        elems[i,1] = graph_map[potential_elems[i,0]]
        elems[i,2] = graph_map[potential_elems[i,1]]

    imports_and_exports.export_ex_coords(nodes[:,:][0:node_kount],groupname,outputfile,'exnode')
    imports_and_exports.export_ex_field(nodal_degrees[0:node_kount], groupname, 'degrees', outputfile + '_degrees', 'exnode')
    imports_and_exports.export_exelem_1d(elems, 'arteries', outputfile)

    return elems, nodes[:,:][0:node_kount],nodal_degrees[0:node_kount]

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

def find_radius_euclidean(DistanceImage,elems,nodes):

    #Distance image contains a Euclidean distance map created from simple itk, so indexing is
    # z, y, x
    
    total_nodes = len(nodes)
    total_elems = len(elems)
    #In terms of node coordinates, y is up down
    radius = np.zeros(total_elems)
    for ne in range(0,total_elems):
        nod1 = elems[ne,1]
        nod2 = elems[ne,2]
        radius1 = DistanceImage[int(nodes[nod1,3]),int(nodes[nod1,2]),int(nodes[nod1,1])]
        radius2 = DistanceImage[int(nodes[nod2,3]),int(nodes[nod2,2]),int(nodes[nod2,1])]
        radius[ne]= (radius1+radius2)/2.
        #if radius[ne]>300.:
        #    radius[ne] = 0.

    print(np.min(radius), np.max(radius), np.mean(radius))

    return radius


#def remove_short_branches(nodes, elems,threshold)
#    num_elem = len(elems)
#    for ne in range(0,num_elem):
#         node1 = nodes[elems[ne,1],1:4]
#         node2 = nodes[elems[ne,2],1:4]
#         length = np.sqrt((node1[0]-node2[0])**2. + (node1[1]-node2[1])**2. (node1[2]-node2[2])**2.)

         
    
    #array = np.delete(array, (i), axis=0)

def find_radius_normal_projection(SkeletonImage, VolumeImage, elems, nodes):
    #, euclid_radii):
    ######
    # Function: Find radius by normal projection
    # Inputs: euclid_radii - radii according to shortest euclidean distance, Ne x 1
    #         SkeletonImage - a logical matrix with skeleton image, Nx x Ny x Nz
    #         VolumeImage - a logical matrix with volume image, Nx x Ny x Nz
    #         elems - and Ne x 3 araay with elems, must be in Strahler order
    #         nodes - Nn x 3 array with node coordinates (it is assumed all nodes are far from the edges of the image) - [Y, X, Z]
    # Outputs: normal_radii - radii according to normal projections, Ne x 1
    #######
    
    # switch nodes to agree with images
    placeHolder=np.copy(nodes[:,0])
    nodes[:, 0]=nodes[:,1]
    nodes[:, 1]=placeHolder

    NumElems = len(elems)
    normal_radii=np.zeros((NumElems))
    normal_radii_std=np.zeros((NumElems)) # within element variation in radius, this isn't needed but has been included just in case

    # Starts and terminal elements and works back up the tree through all the image voxels
    totalErrors=0

    for ne in range(NumElems-1, -1, -1):
        #print(ne)
        coord = np.squeeze(nodes[int(elems[ne, 2]),:]) # start of element
        endCoord = np.squeeze(nodes[int(elems[ne, 1]),:]) # end of element

        count = 0
        elementVoxelSet=np.zeros((1,3))
        elementVoxelSet[count,:]=coord
        errorForElement = 0

        while (~np.prod(coord == endCoord)) & (count < 1000): # ie not another junction voxel + arbitrary absolute check

            # find next coord
            x_start=(int(coord[0]) - 2)
            y_start=int(coord[1] - 2)
            z_start=int(coord[2] - 2)
            large_neighbourhood = np.copy(SkeletonImage[x_start:x_start+5, y_start:y_start+5, z_start:z_start+5]) # 5 x 5 x 5 region around currect coord
            # find next coord
            [nextCoord, error] = find_connected_voxels(large_neighbourhood, coord, endCoord)

            # update image
            SkeletonImage[int(coord[0]), int(coord[1]), int(coord[2])] = 0 # get rid of coordinates once they have been used

            if error:
                errorForElement=1 # will get an error if there is an error at any point in branch

            # update loop
            coord = nextCoord # move up element
            elementVoxelSet=np.append(elementVoxelSet, np.reshape(coord,[1,3]), axis=0)

        SkeletonImage[int(endCoord[0]), int(endCoord[1]), int(endCoord[2])] = 1 # keep the junction voxel as may encounter again

        if count >= 1000:
            print('stuck in loop error') # such as by jumping from one branch to another
            errorForElement = 1

        if errorForElement:
            totalErrors=totalErrors+1
            #shape=np.shape(elementVoxelSet)
            #for k in range (0, shape[0]):
            #    coord2=np.squeeze(elementVoxelSet[k,:])
            #    SkeletonImage[int(coord2[0]),int(coord2[1]),int(coord2[2])]=1 #idea was to return unsuccessfully tracked branches but didnt really work


        # only need to keep inner third of the element for radius calculations
        branchSize = np.shape(elementVoxelSet)
        branchSize=branchSize[0]
        elementVoxelSet = elementVoxelSet[int(np.ceil(branchSize / 3)):int(np.ceil(branchSize - branchSize / 3)), :]
        numVoxels = np.shape(elementVoxelSet)
        numVoxels=numVoxels[0]
        gap = 2 # determines how far we look ahead to get centre line direction (also determined how discretized angles are)

        # Estimate radius by normal projection
        if (numVoxels > gap)&(errorForElement == 0): # can go on to calculate radii

            distanceSet = np.zeros((numVoxels - gap, 1))

            for i in range(0,numVoxels - gap):

                coord1 = np.squeeze(elementVoxelSet[i,:])
                coord2 = np.squeeze(elementVoxelSet[i + gap,:])
                distanceSet[i] = np.mean(find_distances_using_normal(coord1, coord2, VolumeImage))

            # Find mean
            normal_radii[ne] = np.mean(distanceSet)
            normal_radii_std[ne] = np.std(distanceSet)
        else:
            normal_radii[ne] =-1
            normal_radii_std[ne] = -1

    print('Number of Elements that could not successfully be tracked: '+str(totalErrors))

    ## Compare radii to euclidean radii
    #normal_radii[normal_radii < 0] = euclid_radii[normal_radii<0]
    #euclid_radii[euclid_radii == 0] = normal_radii[euclid_radii == 0]
    #euclid_radii[euclid_radii == 0]=np.min(euclid_radii[euclid_radii>0]) # so no chance of div 0
    #difference = abs(normal_radii - euclid_radii)/ euclid_radii
    #cutoff = 0.33333 # distances that are larger than cutoff are not used, 1 means that the distance is the same magnitude as the euclidean distance
    #normal_radii[difference > cutoff] = euclid_radii[difference > cutoff]

    return (normal_radii)
