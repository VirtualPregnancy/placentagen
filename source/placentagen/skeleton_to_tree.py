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
                   #simple itk, so indexing is# z, y, x
                       nodes[node_kount, 1] = coordinates[inew, 2]# coordinates indexed to 'old' node number
                       nodes[node_kount, 2] = coordinates[inew, 1]
                       nodes[node_kount, 3] = coordinates[inew, 0]
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
     
def delete_unused_nodes(nodes,elems):
        
    full_nodes = np.zeros(2*len(elems))
    for ne in range(0,len(elems)):
        full_nodes[2*ne] = elems[ne,1]
        full_nodes[2*ne+1] = elems[ne,2]
    used_nodes = np.unique(full_nodes)
        
    delete_list = []
    node_kount = 0
    node_map = np.zeros(len(nodes))
    for nn in range(0,len(nodes)):
       if nodes[nn,0] not in used_nodes:
          delete_list = np.append(delete_list,nn)
       else:
          node_map[nn]=node_kount#maps old node, to current count 
          node_kount = node_kount+1
    if delete_list != []:        
        delete_list = delete_list.astype(int)
        nodes = np.delete(nodes,delete_list,axis=0)

    for nn in range(0,len(nodes)):
        nodes[nn,0] = nn
    for ne in range(0,len(elems)):
        for nind in range(1,3):
            nn = elems[ne,nind] #Old node number
            elems[ne,nind] = int(node_map[nn])


    return nodes, elems
    
def find_distances_using_normal(coord1, coord2, VolumeImage):

    # get centre line vector
    centre = np.double((coord1 - coord2))/ np.linalg.norm(np.double(coord1 - coord2))
    numSamples = 10
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
    elems_at_node = np.zeros((num_nodes, 20), dtype=int)
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
        radius1 = DistanceImage[int(nodes[nod1,1]),int(nodes[nod1,2]),int(nodes[nod1,3])]
        radius2 = DistanceImage[int(nodes[nod2,1]),int(nodes[nod2,2]),int(nodes[nod2,3])]
        radius[ne]= (radius1+radius2)/2.


    return radius


def find_radius_normal_projection(SegmentationImage, elems, nodes, euclid_radii):
    ######
    # Function: Find radius by normal projection
    # Inputs: euclid_radii - radii according to shortest euclidean distance, Ne x 1
    #         SegmentationImage - a logical matrix with volume image, Nx x Ny x Nz
    #         elems - and Ne x 3 araay with elems, must be in Strahler order
    #         nodes - Nn x 3 array with node coordinates (it is assumed all nodes are far from the edges of the image) - [Y, X, Z]
    # Outputs: normal_radii - radii according to normal projections, Ne x 1
    #######
    
    # switch nodes to agree with images
    #placeHolder=np.copy(nodes[:,0])
    #nodes[:, 0]=nodes[:,1]
    #nodes[:, 1]=placeHolder

    num_elems = len(elems)
    normal_radii=np.zeros((num_elems))
    normal_radii_std=np.zeros((num_elems)) # within element variation in radius, this isn't needed but has been included just in case
    
    for ne in range(0,len(elems)):
        nnod1 = elems[ne,1]
        nnod2 = elems[ne,2]
        coord1 = np.array([nodes[nnod1,3],nodes[nnod1,2],nodes[nnod1,1]]) #Segmentation from simpleITK orders coordinates backwards
        coord2 = np.array([nodes[nnod2,3],nodes[nnod2,2],nodes[nnod2,1]])
        distance_set = find_distances_using_normal(coord1, coord2, SegmentationImage)
        normal_radii[ne] = np.mean(distance_set)

    ## Compare radii to euclidean radii
    normal_radii[normal_radii < 0] = euclid_radii[normal_radii<0]
    euclid_radii[euclid_radii == 0] = normal_radii[euclid_radii == 0]
    euclid_radii[euclid_radii == 0]=np.min(euclid_radii[euclid_radii>0]) # so no chance of div 0
    difference = abs(normal_radii - euclid_radii)/ euclid_radii
    cutoff = 0.33333 # distances that are larger than cutoff are not used, 1 means that the distance is the same magnitude as the euclidean distance
    normal_radii[difference > cutoff] = euclid_radii[difference > cutoff]

    return normal_radii
    
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
    branch_end_elem = elem

    if connected_elems_no>2: #Occasionally we have a single element branch, this is one of them
        #create just the first element in each new branch which will give a seed of parents to give back to our main code
        branch_end_elem = elem
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
    else:
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
                            #create just the first element in each new branch which will give a seed of parents to give back to our main code
                            branch_end_elem = elem
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
                            branch_end_elem = elem
                            continuing = False #This branch does not continue beyond this point
                            break
                        elif connected_elems_no == 2:
                            branch_end_elem = elem
                            #A check to see if a branch has looped back on itself or another existing branch. If it has we want to backtrack through the branch and break it in the middle. However, for now lets see what happens if we treat simply as a terminal.
                            if seen_elements[elems_at_node[first_node][1:connected_elems_no+1]].all():
                                continuing = False
                                cycle = True
                                connected_elems_no = 1 #trick into leaving main loop here

    if new_parents > 0:
       new_parent_list = new_parent_list[0:new_parents]    
    return new_parent_list,cycle,continuing,loop_parent,branch_end_elem
    
def fix_elem_direction(inlet_node,elems,nodes):

    # populate the elems_at_node array listing the elements connected to each node
    num_nodes = len(nodes)
    num_elems = len(elems)
    elems_at_node = np.zeros((num_nodes, 20), dtype=int)
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
    
def remove_small_radius(elems,radii,branch_id,branch_start,threshold):
    #Creating a list of elements with radii less than the threshold
    delete_list = []
    tmp_elems = elems[:,0]
    for nb in range(0,len(branch_start)):
        mean_radius = np.mean(radii[branch_id == nb +1])
        if mean_radius<=threshold:
            tmp_elems1 = tmp_elems[branch_id == nb+1]
            delete_list = np.append(delete_list,tmp_elems1)
    
    
    if delete_list != []:
        delete_list = delete_list.astype(int)
        elems = np.delete(elems,delete_list,axis=0)
        radii = np.delete(radii,delete_list,axis=0)

    #Renumber the elements now some have been deleted      
    for ne in range(0,len(elems)): #renumber elems
        elems[ne,0] = ne
    
    
    print('We have deleted based on radius, removing:',len(delete_list), " Elements")
    return elems, radii
    
    
def remove_order1(nodes,elems,branch_id, radii,length_threshold):
    num_elems = len(elems)
    delete_list = []
    length = np.zeros(len(elems))
    tmp_elems = elems[:,0]
    for ne in range(0,len(elems)):
        node_in = elems[ne,1]
        node_out = elems[ne,2]
        length[ne] = np.sqrt((nodes[node_in,1]-nodes[node_out,1])**2. + (nodes[node_in,2]-nodes[node_out,2])**2.+(nodes[node_in,3]-nodes[node_out,3])**2.)
    
    elem_cnct = pg_utilities.element_connectivity_1D(nodes, elems)
    print("Evaluating orders") 
    orders = analyse_tree.evaluate_orders(nodes, elems)
    print("Number of Strahler Orders",np.max(orders['strahler'])) 
    max_order = np.max(orders['strahler'])-1
    for ne in range(0,num_elems):
        num_upstream = elem_cnct['elem_up'][ne,0]
        if num_upstream == 1:
            upstream_order = orders['strahler'][elem_cnct['elem_up'][ne,1]]
        else:
            upstream_order = 0
        if orders['strahler'][ne]==1 and upstream_order>=max_order:
            branch_length = np.sum(length[branch_id == branch_id[ne]])
            if branch_length<length_threshold:
                tmp_elems1 = tmp_elems[branch_id == branch_id[ne]]
                delete_list = np.append(delete_list,tmp_elems1)
    
    
    if delete_list != []:
        delete_list = delete_list.astype(int)
        elems = np.delete(elems,delete_list,axis=0)
        radii = np.delete(radii,delete_list,axis=0)

    #Renumber the elements now some have been deleted      
    for ne in range(0,len(elems)): #renumber elems
        elems[ne,0] = ne
    print('We have deleted based on small order 1s connected to high orders, removing:',len(delete_list), " Elements")   
    return elems,radii
        
def sort_from_inlet(inlet_node,nodes,elems,branch_id,branch_start,branch_end):#
    
    num_elems = len(elems)
    num_nodes = len(nodes)
    elem_numbers = elems[:,0]
    
    elems_at_node = np.zeros((num_nodes, 20), dtype=int)
    for i in range(0, num_elems):
        elems_at_node[elems[i,1],0] = elems_at_node[elems[i,1],0] + 1
        j = elems_at_node[elems[i,1],0]
        elems_at_node[elems[i,1],j] = elems[i,0] #element number 1 at this node
        elems_at_node[elems[i,2],0] = elems_at_node[elems[i,2],0] + 1
        j = elems_at_node[elems[i,2],0]
        elems_at_node[elems[i,2],j] = elems[i,0]
    seen_elements = np.zeros((num_elems), dtype=bool)
    elem_map = np.zeros((num_elems),dtype=int)
    
    for i in range(0,num_nodes):
        if np.all(nodes[i,1:4]== inlet_node):
            first_node = i
            print("FOUND FIRST NODE",i)
            break
    for i in range(0,num_elems):
        seen_elements[i]=False
        if elems[i,1]==first_node:
            first_branch = np.where(branch_start==i)[0]
            print("FOUND FIRST BRANCH", first_branch)

            
    tmp_elems =  np.zeros((num_elems, 3), dtype=int) #empty array to contain seen elements 
    kount_elems = 0
        
    ne = int(branch_start[first_branch])
    tmp_elems[kount_elems,0] = kount_elems
    tmp_elems[kount_elems,1:3] = elems[ne,1:3]
    elem_map[kount_elems]=ne
    seen_elements[ne] = True
    kount_elems = kount_elems+1
    while ne != int(branch_end[first_branch]):
        if elems_at_node[elems[ne,2],1]==ne:
            ne = int(elems_at_node[elems[ne,2],2])
        else:
            ne = int(elems_at_node[elems[ne,2],1])
        tmp_elems[kount_elems,0] = kount_elems
        tmp_elems[kount_elems,1:3] = elems[ne,1:3]
        elem_map[kount_elems]=ne
        seen_elements[ne]=True
        kount_elems = kount_elems+1


    #subsequent branches can be any order but beneficial to order elems in correct order
    for nb in range(0,len(branch_start)):    
        if nb != first_branch:
            ne = int(branch_start[nb])
            #print(nb,seen_elements[ne],len(elems[branch_id == nb+1]))
            if not seen_elements[ne]:
                tmp_elems[kount_elems,0] = kount_elems
                tmp_elems[kount_elems,1:3] = elems[ne,1:3]
                elem_map[kount_elems]=ne
                kount_elems = kount_elems+1
                seen_elements[ne]=True
            internal_kount = 0
            while ne != int(branch_end[nb]) and internal_kount <len(elems[branch_id == nb+1]):
                if elems_at_node[elems[ne,2],1]==ne:
                    ne = int(elems_at_node[elems[ne,2],2])
                else:
                    ne = int(elems_at_node[elems[ne,2],1])
                internal_kount = internal_kount+1
                if not seen_elements[ne]:
                    tmp_elems[kount_elems,0] = kount_elems
                    tmp_elems[kount_elems,1:3] = elems[ne,1:3]
                    elem_map[kount_elems]=ne
                    kount_elems = kount_elems+1
                    seen_elements[ne]=True
            #print(nb,internal_kount,kount_elems,branch_start[nb],branch_end[nb])
    #print('kount_elems',kount_elems)
    
    return tmp_elems[0:kount_elems,:],elem_map[0:kount_elems]
    
    
######## included to test remove multifurcations functionality 
  
def remove_multiple_elements(geom, elem_connect, type):
    # unpackage information
    max_down = check_multiple(elem_connect)

    while max_down > 2:
        elem_down = elem_connect['elem_down']
        for ne in range(0, len(elem_down)):
            if elem_down[ne][0] > 2:  # more than 2 connected downstream
                if type == 'subtree':
                    geom['nodes'], node2 = extend_node_subtree(ne, geom)  # create new node
                    geom = update_elems(ne, node2, geom, elem_connect)  # create new element and update old
                else:
                    geom['nodes'], node2 = extend_node(ne, geom)  # create new node
                    geom = update_elems(ne, node2, geom, elem_connect)  # create new element and update old
        if type == 'subtree':
            elem_connect = element_connectivity_1D_subtree(geom['nodes'], geom['elems'], 6)
        else:
            elem_connect = element_connectivity_1D(geom['nodes'], geom['elems'], 6)

        max_down = check_multiple(elem_connect)
    num_elems = len(geom['elems'])
    elem_down = elem_connect['elem_down']
    elem_up = elem_connect['elem_up']
    geom['elem_down'] = elem_down[:,0:3]
    geom['elem_up'] = elem_up[:,0:3]


    return geom, elem_connect   
   
   
   
def check_multiple(elem_connect):
    up = elem_connect['elem_up']
    down = elem_connect['elem_down']
    for i in range(0, len(up)):
        if up[i][0] > 2:
            print ('element ', i, 'has >2 upstream elements')
    max_down=0
    count = 0
    for i in range(0, len(down)):
        if down[i][0] > 2:
            print ('element ', i, 'has ', down[i][0], ' downstream elements')
            if max_down < down[i][0]:
                max_down=down[i][0]
            count = count + 1

    print ('number of elements with >2 down stream: ', count)

    return max_down

def extend_node(elem_i, geom):
    nodes=geom['nodes']
    elems=geom['elems']  # nodes and elems initiated
    num_nodes = len(nodes)
    print('original nodes are:', nodes)
    dif = np.zeros(3)   #to store the difference between the two nodes(x,y,z) in a numpy array
    print('dif old:', dif)
    new_node = -1 * np.ones(3)  #new_node initiated
    norm = np.ones(3)    # newly added by VS
    print('new node and the dif array are:', new_node, dif )  #newly added by VS
    print('number of nodes:',num_nodes)      #newly added by VS
    node1 = int(elems[elem_i][1])
    node2 = int(elems[elem_i][2])  # node at other end of the element

    for i in range(0, 3):
        # assuming nodes starts index = node number (start at 0)
        #dif[i] = np.abs(nodes[node1][i] - nodes[node2][i])  # store difference of xyz
        dif[i] = (nodes[node2][i] - nodes[node1][i])
    print('storing the dif:',dif)    #newly added by VS

    #for i in range(0, 3):
    #norm[i] = np.linalg.norm(dif[i])  # newly added by VS
    norm = np.linalg.norm(dif)
    print('normalized vector:', norm)
    #print('dif and norm product:', (dif*norm))
    print('dif and norm division (unit vector):', (dif/norm))
    new_node = (nodes[node2] + ((dif/norm)*1e-3))
    print('new_node created at:', new_node)  # newly added by VS
    nodes = np.vstack((nodes, new_node))  # new node created in the same axis as the old one
    node2 = int(num_nodes)
    print('nodes, node1 and node2:', nodes, node1, node2) #newly added by VS
    return nodes, node2
            
   
def update_elems(elem_i, node2, geom, elem_connect):
    #unpack inputs
    elem_up = elem_connect['elem_up']
    elem_down = elem_connect['elem_down']
    elems=geom['elems']
    nodes=geom['nodes']

    num_elem = len(elems)
    new_elem = -1 * np.ones(3)
    node1=int(elems[elem_i][2]) #node other end of elem

    # create new elem connecting node1 and new node2
    new_elem[0] = num_elem  # elem numbering starts from 0; [0 1 2] num_elem = 3; new elem = 3
    new_elem[1] = node1
    new_elem[2] = node2

    # add new element to end
    elems = np.vstack((elems, new_elem))

    # update after second downstream element with new node
    for i in range(2,elem_down[elem_i][0]+1):
        old_elem = elem_down[elem_i][i]  # first down stream element
        elems[old_elem][1] = node2  # change starting node of old_elem to new node2

    geom['elems']=elems

    # add copy of node1 geom for node2 at end
    for item in geom.keys():
        current = geom[item]
        # print 'key:', item
        #print 'current', current
        # print 'current[ne]', current[ne]
        if item == 'nodes' or item == 'elems':
            continue #node and element already appended
        elif item == 'length': #radii 1D array
            new_length = find_length_single(nodes, node1, node2)
            current = np.hstack((current, new_length))
        else:
            current = np.hstack((current, current[elem_i]))
        geom[item]=current

    return geom     
   


   

