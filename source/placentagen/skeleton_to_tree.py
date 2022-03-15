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
    potential_elems = np.empty((0, 2), dtype=int)
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
            attached_elems = pixel_graph.indices[pixel_graph.indptr[i]:pixel_graph.indptr[i+1]]
            for k in range(0,num_attached):
                if attached_elems[k] != i: #avoiding self connections (simple loops)
                    potential_elems = np.append(potential_elems, np.zeros((1, 2), dtype=int), axis=0)
                    potential_elems[count_new_vessel,0] = int(min(i,attached_elems[k]))
                    potential_elems[count_new_vessel,1]= int(max(i,attached_elems[k]))
                    count_new_vessel = count_new_vessel + 1


    potential_elems = np.unique(potential_elems, axis=0)

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
    
def cut_loops(elems,nodes,branch_id,branch_start,branch_end,cycles,radii):
     delete_list = []
     for i in range(0,len(branch_start)):
         if cycles[i]:
             tmp_elems = elems[branch_id == i+1]
             tmp_radii = radii[branch_id == i+1]
             where=np.argmin(tmp_radii)
             delete_list = np.append(delete_list,tmp_elems[where,0])
             
     if delete_list != []:        
         delete_list = delete_list.astype(int)
         elems = np.delete(elems,delete_list,axis=0)
         radii = np.delete(radii,delete_list,axis=0)
     
     for ne in range(0,len(elems)): #renumber elems
        elems[ne,0] = ne
             
     return elems, radii

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
    
def find_inlet_auto(elems,nodes,radii,length_threshold):
    # populate the elems_at_node array listing the elements connected to each node
    num_nodes = len(nodes)
    num_elems = len(elems)
    elems_at_node = np.zeros((num_nodes, 10), dtype=int)
    possible_inlets = []
    for i in range(0, num_elems):
        elems_at_node[elems[i,1],0] = elems_at_node[elems[i,1],0] + 1
        j = elems_at_node[elems[i,1],0]
        elems_at_node[elems[i,1],j] = elems[i,0] #element number 1 at this node
        elems_at_node[elems[i,2],0] = elems_at_node[elems[i,2],0] + 1
        j = elems_at_node[elems[i,2],0]
        elems_at_node[elems[i,2],j] = elems[i,0]
        
    biggest_branch = 0 
    rad_max = 0.0
    for i in range(0,num_nodes):
        if i ==3218:
           print('branch branch', elems_at_node[i,0],elems_at_node[i,1],elems_at_node[i,2],nodes[i,0])

        if elems_at_node[i,0]==1:
            branch_length = 0.
            branch_radius = 0.
            num_elems = 0.

            ne = elems_at_node[i,1]
            if i == elems[ne,1]:
                first_node = elems[ne,2]
            else:
                first_node = elems[ne,1]
            connected_elems_no = elems_at_node[first_node,0]
            elem_list = []
            while connected_elems_no == 2:
               for k in range(0, connected_elems_no): #Goes through elems
                   elem = elems_at_node[first_node,k + 1]
                   if elem != ne and elem not in elem_list:
                      elem_list = np.append(elem_list,elem)
                      node_in = elems[elem,1]
                      node_out = elems[elem,2]
                      branch_length = branch_length + np.sqrt((nodes[node_in,1]-nodes[node_out,1])**2. + (nodes[node_in,2]-nodes[node_out,2])**2.+(nodes[node_in,3]-nodes[node_out,3])**2.)
                      branch_radius = branch_radius + radii[elem]
                      num_elems = num_elems+1
                      #moving through the elements
                      if first_node == elems[elem,1]:
                          first_node = elems[elem,2]
                      else:
                          first_node = elems[elem,1]
                      ne = elem
                      connected_elems_no = elems_at_node[first_node,0]
            if num_elems > 0:
                branch_radius = branch_radius/num_elems #Mean branch radius
            else:
                branch_radius = 0.
                
                
            if branch_length >= length_threshold:#exclude short branches
                if(branch_radius>=rad_max):
                    rad_max = branch_radius
                    biggest_branch = i

            
    print('Determined inlet, node: ', biggest_branch)
    return nodes[biggest_branch,:]

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


    return radius


def find_radius_normal_projection(SkeletonImage, VolumeImage, elems, nodes):
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
    
def fix_branch_direction(first_node,elems_at_node,elems,seen_elements,branch_id,branches,old_parent_list,inlet_branch):
    #This routine should correct branch direction correctly in a generated tree
    maxification = np.max(elems_at_node[:,0])
    new_parents = 0
    new_parent_list = np.zeros(int(maxification),dtype = int)
    continuing = False
    elem = elems_at_node[first_node][1]
    connected_elems_no = elems_at_node[first_node][0]  # number of elements connected to this one
    branch_starts_at = first_node
    loop_parent = len(elems)+1
    cycle = False
    
    while connected_elems_no ==2 or inlet_branch: #continuing branch
        inlet_branch = False #No longer an inlet branch, if it was
        #Checking whether first node in an element is a parent, and its not where the branch starts (should be detecting loops)
        if first_node in np.asarray(old_parent_list) and first_node != branch_starts_at:
            connected_elems_no=1
            loop_parent = first_node
        else:
            if first_node == branch_starts_at and seen_elements[elems_at_node[first_node][1:connected_elems_no+1]].all():
               connected_elems_no = 1                
            for i in range(0, connected_elems_no): #Goes through elems
                elem = elems_at_node[first_node][i + 1]  # elements start at column index 1
                if not seen_elements[elem]:
                    branch_id[elem] = branches #This is a new element and belongs in this branch
                    if elems[elem][1] != first_node: #Swap nodes if going 'wrong way'
                        # swap nodes
                        elems[elem][2] = elems[elem][1]
                        elems[elem][1] = first_node
                    seen_elements[elem] = True
                    first_node = elems[elem][2] #New first node is the second node of the element
                    connected_elems_no = elems_at_node[first_node][0]  # number of elements connected to this one
                    #Now we check if we've reached the end of a branch (either a branch point, a termination, or we have come back in a loop around)
                    if connected_elems_no >= 3: #Bifurcation or morefication? point
                        # we are going to create just the first element in each new branch which will give a seed of parents to give back to our main code
                        new_parents = 0
                        for i in range(0, connected_elems_no):
                            elem = elems_at_node[first_node][i + 1]  # elements start at column index 1
                            if not seen_elements[elem]:
                                new_parent_list[new_parents] = elem
                                new_parents = new_parents + 1
                                branch_id[elem] = branches + new_parents
                                if elems[elem][1] != first_node:
                                    # swap nodes
                                    elems[elem][2] = elems[elem][1]
                                    elems[elem][1] = first_node
                            seen_elements[elem] = True
                            continuing = True #We are going to grow some more branches from here
                    elif connected_elems_no == 1: #Terminal
                        continuing = False #This branch does not continue beyond this point
                        break
                    elif connected_elems_no == 2: 
                         #A check to see if a branch has looped back on itself or another existing branch. If it has we want to backtrack through the branch and break it in the middle. However, for now lets see what happens if we treat simply as a terminal.
                        if seen_elements[elems_at_node[first_node][1:connected_elems_no+1]].all():
                            continuing = False
                            cycle = True
                            connected_elems_no = 1 #trick into leaving main loop here
    if new_parents > 0:
       new_parent_list = new_parent_list[0:new_parents]    
    return new_parent_list,cycle,continuing,loop_parent,elem
    
def fix_elem_direction(inlet_node,elems,nodes):

    # populate the elems_at_node array listing the elements connected to each node
    num_nodes = len(nodes)
    num_elems = len(elems)
    elems_at_node = np.zeros((num_nodes, 10), dtype=int)
    branch_id = np.zeros((num_elems),dtype = int)
    branch_start = []
    branch_end = []
    cycle_bool = []
    for i in range(0, num_elems):
        elems_at_node[elems[i][1]][0] = elems_at_node[elems[i][1]][0] + 1
        j = elems_at_node[elems[i][1]][0]
        elems_at_node[elems[i][1]][j] = elems[i][0]
        elems_at_node[elems[i][2]][0] = elems_at_node[elems[i][2]][0] + 1
        j = elems_at_node[elems[i][2]][0]
        elems_at_node[elems[i][2]][j] = elems[i][0]
        
    
    for i in range(0,num_nodes):
        if np.all(nodes[i,1:4]== inlet_node):
            first_node = i
            print("FOUND FIRST NODE",i)
            break

    ############
    seen_elements = np.zeros((num_elems), dtype=bool)
    branches = 1
    continuing = True
    old_parent_list = first_node
    loop_list = np.zeros(1,dtype=int)
    loop_list[0] = (num_elems+1)
  
    branch_start = np.append(branch_start, elems_at_node[first_node,1]) #First element in branch

    [new_parent_list,cycle,continuing,loop_parent,branch_end_elem] = fix_branch_direction(first_node, elems_at_node, elems, seen_elements,branch_id,branches,old_parent_list,True)
    branch_end = np.append(branch_end, branch_end_elem)
    cycle_bool = np.append(cycle_bool, cycle)
    
    while len(new_parent_list)>0:
        if len(new_parent_list) > 0:
            new_parent_list2 = []
            for parent in range(0,len(new_parent_list)):
                second_node = elems[new_parent_list[parent],2]
                if second_node not in loop_list:
                    branches = branches + 1
                    branch_start  = np.append(branch_start,new_parent_list[parent])
                    branch_id[new_parent_list[parent]] = branches
                    
                    [branch_list,cycle, continuing,loop_parent,branch_end_elem] = fix_branch_direction(second_node, elems_at_node, elems, seen_elements,
                                                           branch_id, branches,elems[new_parent_list,2],False)
                    branch_end = np.append(branch_end, branch_end_elem)
                    cycle_bool = np.append(cycle_bool, cycle)

                if loop_parent< num_elems:
                    loop_list = np.append(loop_list, [loop_parent], axis=0) #Removes simple loops
                if continuing:
                    new_parent_list2 = np.append(new_parent_list2,branch_list,axis = 0)
                #print(continuing,new_parent_list2)
            if(len(new_parent_list2)>0):
                new_parent_list = new_parent_list2.astype(int)
            else:
                new_parent_list = []


    print("Number of branches identified: ", len(branch_start))
    
    return elems,branch_id,branch_start,branch_end,cycle_bool,seen_elements
  
def remove_disconnected(elems, euclid_radii, branch_id, seen):
    #remove disconnected elements (not part of fix elem direction as there may be two inlets)
    elems = elems[seen == True]
    euclid_radii = euclid_radii[seen == True]
    branch_id = branch_id[seen == True]
    for ne in range(0,len(elems)): #renumber elems
        elems[ne,0] = ne
        
    return elems, euclid_radii,branch_id
    
def remove_small_radius(elems,radii,threshold):
    #Creating a list of elements with radii less than the threshold
    delete_list = []
    for ne in range(0,len(elems)):
        if radii[ne]<=threshold:
            delete_list = np.append(delete_list,int(ne))
    
    
    if delete_list != []:
        delete_list = delete_list.astype(int)
        elems = np.delete(elems,delete_list,axis=0)
        radii = np.delete(radii,delete_list,axis=0)

    #Renumber the elements now some have been deleted      
    for ne in range(0,len(elems)): #renumber elems
        elems[ne,0] = ne
    
    
    print('We have deleted based on radius, removing:',len(delete_list), " Elements")
    return elems, radii
