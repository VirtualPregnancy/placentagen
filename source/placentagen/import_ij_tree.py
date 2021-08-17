import numpy as np
from .pg_utilities import remove_rows,row_swap_2d,row_swap_1d
from ismember import ismember
from matplotlib import pyplot as plt
from tabulate import tabulate
import scipy

######
# Function: takes data from the csv and converts it to arrays
# Inputs: data_file - generated from the panadas read_csv function, containing results from imageJ image analysis
#         Arrays - a group of arrays each with length N for their first axis
# Outputs: nodes - an M x 3 array giving cartesian coordinates (x,y,z) for the node locations in the tree
#         elems - an N x 3 array, the first colum in the element number, the second two columns are the index of the start and end node
#         radii, length, euclidean_length - there are all an N x 1 array containing a property for each element
######

def sort_data(data_file):
    # get rid of any skeletons other than the main one
    data_file = data_file[data_file.SkeletonID == 1]

    # get skeleton properties as arrays
    euclid_length = data_file.Euclideandistance.values
    length = data_file.Branchlength.values
    radii = data_file.averageintensityinner3rd.values
    radii = radii

    #for i in range(0,len(radii)):
    #    if radii[i]==0:
    #        radii[i] = np.min(radii[radii>0])


    # get elem and node data
    data_file = data_file.drop(['SkeletonID', 'Branchlength', 'averageintensityinner3rd', 'Euclideandistance'], axis=1)
    data_file = data_file.values
    (elems, nodes) = sort_elements(data_file[:, 0:3], data_file[:, 3:6])

    # get rid of dud elements
    (elems, [length, euclid_length, radii]) = remove_rows(elems, [length, euclid_length, radii])


    print('Finalised nodes and elements', len(nodes), 'nodes', len(elems),'elements')
    return {'nodes': nodes, 'elems': elems, 'radii': radii, 'length': length, 'euclidean length': euclid_length}

def sort_elements(v1, v2):
    #v1 and v2 are the start and end nodes of each "braih"
    Nelem = len(v1) #number of elements
    elems = np.zeros([Nelem, 3])
    #nodes = np.zeros([Nelem * 2, 3])  # max number of nodes possible

    iN = 0  # node index
    combined_nodes = np.concatenate((v1,v2))
    unique_combined_nodes = np.unique(combined_nodes,axis=0)
    nodes = np.reshape(unique_combined_nodes,(len(unique_combined_nodes),3))
    Iloc_in, index_in = ismember(v1, nodes, 'rows')
    Iloc_out, index_out = ismember(v2, nodes, 'rows')

    for iE in range(0, Nelem):
        elems[iE, 0] = int(iE)
        elems[iE, 1] = int(index_in[iE])
        elems[iE, 2] = int(index_out[iE])

    return (elems, nodes)

######
# Function: rearranges elems (and corresponding properties) according to their strahler order, to be compatible with placentagen functions
# Inputs: geom - contains elems, nodes and other properties of the skeleton
#         inlet_loc - the coordinates of the parent node for the entire tree (if known)
#         find_inlet_loc - a boolean variable specifying whether to use inlet location provided (0) or to find the inlet location automatically (1)
# Outputs: geom - contains elems and properties, reordered according to strahler order so that no element can be higher in the element list than a higher order branch
######

def arrange_by_strahler_order(geom, find_inlet_loc, inlet_loc):
    # set up arrays
    nodes = geom['nodes']
    elem_properties = np.column_stack([geom['radii'], geom['length'], geom['euclidean length'], geom['elems']])
    elems = np.copy(geom['elems'])  # as elems is altered in this function
    elems = elems[:, 1:3]  # get rid of first column which means nothing
    radii = geom['radii']

    Ne = len(elems)
    Nn = len(nodes)
    elem_properties_new = np.zeros([Ne, 6])

    # find parent node
    (elems, elem_properties) = find_parent_node(find_inlet_loc, inlet_loc, nodes, radii, elems, elem_properties)

    # loop through by strahler order
    counter_new = 0
    counter = 1
    while (counter < Ne):

        # find elements which are terminal
        terminal_elems = np.zeros([Ne, 1])

        # go through each node
        for i in range(0, Nn + 1):

            # find number of occurrences of the node
            places = np.where(elems == i)
            ind1 = places[0]
            ind2 = places[1]

            if (len(ind1) == 1) and ((ind1[0]) != 0):  # if occurs once, then element is terminal (avoids root element)

                ind1 = ind1[0]
                ind2 = ind2[0]

                # swap to ensure element points right way
                if ind2 == 0:
                    elems[ind1, :] = row_swap_1d(np.squeeze(elems[ind1, :]), 1, 0)
                    elem_properties[ind1, 4:6] = row_swap_1d(np.squeeze(elem_properties[ind1, 4:6]), 1, 0)

                # assign element under the new element ordering scheme
                elem_properties_new[counter_new, :] = elem_properties[ind1, :]
                counter_new = counter_new + 1

                terminal_elems[ind1] = 1

                # join up element with upstream elements
                nodeNumNew = elems[ind1, 0]  # this is node number at other end of element
                nodeNum = i
                places = np.where(elems == nodeNumNew)  # find where the new node occurs
                ind1 = places[0]
                ind2 = places[1]

                counter2 = 1

                while ((len(ind1) == 2) & (counter2 < Ne)):  # as can only be present twice if a joining node

                    # see if branch joins to yet another branch, that we haven't yet encountered (i.e. not nodeNum)
                    if (elems[ind1[0], ~ind2[0]] == nodeNum):
                        k = 1
                    else:
                        k = 0
                    terminal_elems[ind1[k]] = 1  # label terminal_elems as joining elements

                    # switch the way element points
                    if (ind2[k] == 0):
                        elems[ind1[k], :] = row_swap_1d(np.squeeze(elems[ind1[k], :]), 1, 0)
                        elem_properties[ind1[k], 4:6] = row_swap_1d(np.squeeze(elem_properties[ind1[k], 4:6]), 1, 0)

                    nodeNum = nodeNumNew
                    nodeNumNew = elems[ind1[k], 0]

                    # assign new order
                    elem_properties_new[counter_new, :] = elem_properties[ind1[k], :]
                    counter_new = counter_new + 1

                    # update loop criteria
                    places = np.where(elems == nodeNumNew)
                    ind1 = places[0]
                    ind2 = places[1]
                    counter2 = counter2 + 1

        # update elems to 'get rid of' terminal elements from the list
        terminal_elems[0] = 0  # the root node can never be terminal
        terminal_elems_pair = np.column_stack([terminal_elems, terminal_elems])
        elems[terminal_elems_pair == 1] = -1

        # loop exit criteria
        places = np.where(terminal_elems == 1)
        places = places[1]
        if len(places) == 0:
            counter = Ne + 1
        counter = counter + 1

    # assign root element in new order systems
    elem_properties_new[Ne - 1, :] = elem_properties[0, :]

    # reduce size due to elements removed
    elem_properties_new = elem_properties_new[0:Ne, :]
    # reverse order
    elem_properties_new = np.flip(elem_properties_new, 0)

    elems = geom['elems']
    elems = elems[0:Ne, :]
    elems[:, 1:3] = elem_properties_new[:, 4:6]
    radii = elem_properties_new[:, 0]
    lengths = elem_properties_new[:, 1]
    euclid_lengths = elem_properties_new[:, 2]

    return {'elems': elems, 'radii': radii, 'length': lengths, 'euclidean length': euclid_lengths, 'nodes': nodes}


######
# Function: find parent node in array either from given coordinates or by finding the fattest terminal branch
# Inputs: nodes - an M x 3 list of node coordinates
#        radii - N x 1 array with radius of each element
#         elems - an Nx2(!!!) array with node indices for start and end of node
#         elem_properties - an N x K array, with each row containing various element properties (radii etc.)
#         inlet_loc - the coordinates of the parent node for the entire tree (if known)
#         find_inlet_loc - a boolean variable specifying whether to use inlet location provided (0) or to find the inlet location automatically (1)
# Outputs: elems and elem_properties updates so that inlet element is the first element in the list
######
def find_parent_node(find_inlet_loc, inlet_loc, nodes, radii, elems, elem_properties):
    # will define inlet as terminal element of largest radius
    if find_inlet_loc == 1:

        maxRad = -1
        # go through each node
        for i in range(0, len(nodes) + 1):

            # find number of occurrences of the node
            places = np.where(elems == i)
            ind1 = places[0]
            ind2 = places[1]

            if (len(ind1) == 1):  # if occurs once, then element is terminal (avoids root element)

                ind1 = ind1[0]
                ind2 = ind2[0]
                radius = radii[ind1]

                if radius > maxRad:
                    maxRad = radius
                    maxRadInd = i
                    Ne_root = ind1

        inlet_loc = np.squeeze(nodes[maxRadInd, :])
        Nn_root = maxRadInd
    # find root node and element from coordinates provided
    else:
        Nn_root = ismember(inlet_loc, nodes,'rows')
        if (Nn_root == -1):
            print("Warning, root node not located")

    print('Inlet Coordinates:' + str(inlet_loc) + ' ' + str(Ne_root))

    # find root element
    Ne_place = elems[Ne_root]
    #Ne_root = Ne_place[0]  # only need first index
    #if len(Ne_root) > 1:
    #    print("Warning, root node is associated with multiple elements")
    #if len(Ne_root) == 0:
    #    print("Warning, no root element located")
    #Ne_root = Ne_root[0]

    # make root element the first element
    elems = row_swap_2d(elems, 0, Ne_root)
    elem_properties = row_swap_2d(elem_properties, 0, Ne_root)

    # get element pointing right way
    if (np.squeeze(Ne_place[1]) != 0):
        elems[0, :] = row_swap_1d(np.squeeze(elems[0, :]), 1, 0)
        elem_properties[0, 4:6] = row_swap_1d(np.squeeze(elem_properties[0, 4:6]), 1, 0)

    return (elems, elem_properties)

def find_maximum_joins(elems):
    elems = np.concatenate([np.squeeze(elems[:, 1]), np.squeeze(elems[:, 2])])
    elems = elems.astype(int)
    result = np.bincount(elems)
    Nc = (max(result)) + 1

    # Warning if detect an unusual value
    if Nc > 12:
        print('Warning, large number of elements at one node: ' + str(Nc))
        Nc = 12

    return result, Nc

def element_connectivity_1D_rachel(node_loc, elems, Nc):
    # Initialise connectivity arrays
    num_elems = len(elems)
    elem_upstream = np.zeros((num_elems, Nc), dtype=int)
    elem_downstream = np.zeros((num_elems, Nc), dtype=int)

    num_nodes = len(node_loc)
    elems_at_node = np.zeros((num_nodes, Nc), dtype=int)

    # determine elements that are associated with each node
    for ne in range(0, num_elems):
        for nn in range(1, 3):
            nnod = int(elems[ne][nn])
            elems_at_node[nnod][0] = elems_at_node[nnod][0] + 1
            elems_at_node[nnod][elems_at_node[nnod][0]] = ne

    # assign connectivity
    for ne in range(0, num_elems):
        nnod2 = int(elems[ne][2])

        for noelem in range(1, elems_at_node[nnod2][0] + 1):
            ne2 = elems_at_node[nnod2][noelem]

            if ne2 != ne:
                elem_upstream[ne2][0] = elem_upstream[ne2][0] + 1
                elem_upstream[ne2][elem_upstream[ne2][0]] = ne
                elem_downstream[ne][0] = elem_downstream[ne][0] + 1
                elem_downstream[ne][elem_downstream[ne][0]] = ne2

    return {'elem_up': elem_upstream, 'elem_down': elem_downstream}

def evaluate_orders_rachel(elems, elem_connect):
    num_elems = len(elems)

    elem_upstream = elem_connect['elem_up']
    elem_downstream = elem_connect['elem_down']

    # Initialise order definition arrays
    strahler = np.zeros(len(elems), dtype=int)
    horsfield = np.zeros(len(elems), dtype=int)
    generation = np.zeros(len(elems), dtype=int)

    # Calculate generation of each element
    maxgen = 1  # Maximum possible generation
    for ne in range(0, num_elems):
        ne0 = elem_upstream[ne][1]
        if ne0 != 0:
            # Calculate parent generation
            n_generation = generation[ne0]
            if elem_downstream[ne0][0] == 1:
                # Continuation of previous element
                generation[ne] = n_generation
            elif elem_downstream[ne0][0] >= 2:
                # Bifurcation (or morefurcation)
                generation[ne] = n_generation + 1
        else:
            generation[ne] = 1  # Inlet
        maxgen = np.maximum(maxgen, generation[ne])

    # Now need to loop backwards to do ordering systems
    for ne in range(num_elems - 1, -1, -1):
        n_horsfield = np.maximum(horsfield[ne], 1)
        n_children = elem_downstream[ne][0]
        if n_children == 1:
            if generation[elem_downstream[ne][1]] == 0:
                n_children = 0
        temp_strahler = 0
        strahler_add = 1
        if n_children >= 2:  # Bifurcation downstream
            temp_strahler = strahler[elem_downstream[ne][1]]  # first daughter
            for noelem in range(1, n_children + 1):
                ne2 = elem_downstream[ne][noelem]
                temp_horsfield = horsfield[ne2]
                if temp_horsfield > n_horsfield:
                    n_horsfield = temp_horsfield
                if strahler[ne2] < temp_strahler:
                    strahler_add = 0
                elif strahler[ne2] > temp_strahler:
                    strahler_add = 0
                    temp_strahler = strahler[ne2]  # strahler of highest daughter
            n_horsfield = n_horsfield + 1
        elif n_children == 1:
            ne2 = elem_downstream[ne][1]  # element no of daughter
            n_horsfield = horsfield[ne2]
            strahler_add = strahler[ne2]
        horsfield[ne] = n_horsfield
        strahler[ne] = temp_strahler + strahler_add

    return {'strahler': strahler, 'horsfield': horsfield, 'generation': generation}

def remove_multifurcations(geom,elem_connect,num_con):
   count_corrected = 1000
   while count_corrected > 0:
        #loop though elements and work out which ones have more than 2
        count_corrected = 0
        for i in range(0,len(geom["elems"])):
            if elem_connect['elem_down'][i,0] > 3:
                geom["nodes"],node2 = extend_node(i,geom)
                geom = update_elems(i,node2,geom,elem_connect)
                count_corrected = count_corrected + 1
        elem_connect = element_connectivity_1D_rachel(geom['nodes'], geom['elems'], 12)

   return geom,elem_connect



######
# Function: adds another node to the end of node list that is the input node extended slightly in the longest axis of the
#           associated element
# Inputs: elem_i - index of element associated with >2 downstream elements
#		  nodes - node info in a skeleton where nodes=[x, y, z]
#		  elems - element info where elems=[elem_no, node1, node2] elements and nodes start from 0
# Output: updated nodes - new node appended to end of input nodes
#         node2 - the node number of the new node
######
def extend_node(elem_i, geom):
    nodes=geom['nodes']
    elems=geom['elems']
    num_nodes = len(nodes)
    dif = np.zeros(3)
    new_node = -1 * np.ones(3)

    node1 = int(elems[elem_i][1])
    node2 = int(elems[elem_i][2])  # node at other end of the element
    for i in range(0, 3):
        # assuming nodes starts index = node number (start at 0)
        dif[i] = np.abs(nodes[node1][i] - nodes[node2][i])  # store difference of xyz
    max_i = np.argmax(dif)  # longest axis (x, y or z)
    for i in range(0, 3):
        new_node[i] = nodes[node1][i]  # replicate old node
        if i == max_i:
            if nodes[node2][i] < 0:
                new_node[i] = nodes[node2][i] - 1e-10  # extend node slightly in longest axis
            else:
                new_node[i] = nodes[node2][i] + 1e-10
    # add new node to end
    nodes = np.vstack((nodes, new_node))
    node2 = int(num_nodes)

    return nodes, node2

####
# Function: creates new element connecting node1 and node2 and updates existing elements that start at node1
# Inputs: elem_i - elem that has more than 2 downstream elements
#		  node2 - new node number
#		  geom - array containing tree information eg. elems, nodes, radii, length, euclidean length
#         elem_connect - array containing element connectivity information of elem_down and elem_up
# Outputs: geom_new - updated input geom array with extra element and node
######
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


def find_length_single(nodes,node1,node2):

    d_x = nodes[node1][0] - nodes[node2][0]
    d_y = nodes[node1][1] - nodes[node2][1]
    d_z = nodes[node1][2] - nodes[node2][2]
    distance = np.sqrt(d_x ** 2 + d_y ** 2 + d_z ** 2)

    return distance


#####
# Function: finds properties of according to each Branch of the tree, where a branch is a set of elements with the
#          same Strahler order
# Inputs: geom - contains elems, and various element properties (length, radius etc.)
#         order - contains strahler order and generation of each element
#         elem_up - contains index of upstream elements for each element
# Outputs: branchGeom: contains the properties arrange in arrays according to each branch:
#           radius / length / euclidean length / strahler order: all M x 1 arrays where M is number of branches
#          branches: an N x 1 array where N is the number of elements, contains branch number of each element
######

def arrange_by_branches(geom, elem_up, order):

    # find branches
    Ne = len(order)
    branches = np.zeros(Ne)
    branchNum = 1

    for i in range(0, Ne):
        if order[i] != order[elem_up[i, 1]]:  # belongs to new branch
            branchNum = branchNum + 1
        branches[i] = branchNum

    Nb = int(max(branches))

    # sort results into branch groups
    lengths = geom['length']
    radii = geom['radii']
    nodes= geom['nodes']
    elems = geom['elems']

    branchRad = np.zeros(Nb)
    branchLen = np.zeros(Nb)
    branchEucLen = np.zeros(Nb)
    branchOrder = -1. * np.ones(Nb)

    for i in range(0, Nb):
        branchElements = np.where(branches == i+1) #find elements belonging to branch number
        branchElements = branchElements[0]

        for j in range(0, len(branchElements)): #go through all elements in branch
            ne = branchElements[j]

            branchOrder[i] = order[ne]
            branchLen[i] = branchLen[i] + lengths[ne]
            branchRad[i] = branchRad[i] + radii[ne]

        branchRad[i] = branchRad[i] / len(branchElements) # to get average radius

        startNode=nodes[int(elems[branchElements[0],1]),:]
        endNode=nodes[int(elems[branchElements[len(branchElements)-1],2]),:]

        branchEucLen[i]=np.sqrt(np.sum(np.square(startNode-endNode)))

    return {'radii': branchRad, 'length': branchLen, 'euclidean length': branchEucLen, 'order': branchOrder,
            'branches': branches}


######
# Function: find branch angles + L/LParent & D/Dparent
#          scale all results into mm and degrees
# Inputs: geom - contains elems, and various element properties (length, radius etc.)
#         orders - contains strahler order and generation of each element
#         elem_connect - contains upstream and downstream elements for each element
#         branchGeom - contains branch properties (length, radius, etc.)
#         voxelSize - for conversion to mm (must be isotropic)
#         conversionFactor - to scale radii correction, printed in log of ImageJ during MySkeletonizationProcess
# Outputs: geom and branchGeom are altered so all there arrays are in correct units (except nodes, and radii_unscaled, which remain in voxels) ##################
#          seg_angles - angle (radians) at each element junction in the tree assigned to each element according to how it branches from its parent
#          diam_ratio - ratio of length/diameter of each branch, accounting for multi-segment branches
#          length_ratio - ratio of parent / child lengths, accounting for multi-segment branches
#          diam_ratio / length_ratio / branch_angles are the same but for whole branches
######

def find_branch_angles(geom, orders, elem_connect, branchGeom, voxelSize, conversionFactor):

    # unpackage inputs
    nodes = geom['nodes']
    elems = geom['elems']
    elems = elems[:, 1:3]  # get rid of useless first column
    radii = geom['radii']
    lengths = geom['length']

    branches = branchGeom['branches']
    branchRad = branchGeom['radii']
    branchLen = branchGeom['length']

    strahler = orders['strahler']
    generations = orders['generation']

    elem_up = elem_connect['elem_up']

    # new arrays
    num_elems = len(elems)
    num_branches = len(branchRad)

    branch_angles = -1. * np.ones(num_branches) # results by branch (Strahler)
    diam_ratio_branch = -1. * np.ones(num_branches)
    length_ratio_branch = -1. * np.ones(num_branches)

    diam_ratio = -1. * np.ones(num_elems)  # results by generation
    length_ratio = -1. * np.ones(num_elems)
    seg_angles = -1. * np.ones(num_elems)

    # find results for each element (ignoring parent element)
    for ne in range(1, num_elems):

        neUp = elem_up[ne, 1] # find parent

        if (generations[neUp] < generations[ne]): # there is branching but not necessarily a new strahler branch

            # parent node
            endNode = int(elems[neUp, 0])
            startNode = int(elems[neUp, 1])
            v_parent = nodes[endNode, :] - nodes[startNode, :]
            v_parent = v_parent / np.linalg.norm(v_parent)

            d_parent = 2 * radii[neUp]
            L_parent = lengths[neUp]

            # daughter
            endNode = int(elems[ne, 1])
            startNode = int(elems[ne, 0])
            v_daughter = nodes[startNode, :] - nodes[endNode, :]
            v_daughter = v_daughter / np.linalg.norm(v_daughter)

            d_daughter = 2 * radii[ne]
            L_daughter = lengths[ne]

            # calculate angle
            dotProd = np.dot(v_parent, v_daughter)
            if abs(dotProd <= 1):
                angle=np.arccos(dotProd)
                seg_angles[ne] = angle
            else:
                angle=-1
                print('Angle Error, element: ' + str(ne))

            if d_parent != 0:
                diam_ratio[ne] = d_daughter/ d_parent
            if L_parent != 0:
                length_ratio[ne] = L_daughter / L_parent

            if (strahler[neUp] > strahler[ne]): #then this also is a new strahler branch

                # assign results
                branchNum = int(branches[ne])-1
                parentBranch = int(branches[neUp])-1

                branch_angles[branchNum] = angle

                if branchRad[parentBranch] != 0:
                    diam_ratio_branch[branchNum] = branchRad[branchNum] / branchRad[parentBranch]
                if branchLen[parentBranch] != 0:
                    length_ratio_branch[branchNum] = branchLen[branchNum] / branchLen[parentBranch]

    # scale results into mm and degrees & package them up
    geom['radii'] = geom['radii'] / conversionFactor * voxelSize
    geom['length'] = geom['length'] * voxelSize
    geom['nodes'] = geom['nodes'] * voxelSize
    geom['euclidean length'] = geom['euclidean length'] * voxelSize
    geom['branch angles'] = seg_angles * 180 / np.pi
    geom['diam_ratio'] = diam_ratio
    geom['length_ratio'] = length_ratio

    branchGeom['radii']= branchGeom['radii']/ conversionFactor
    branchGeom['radii'] = branchGeom['radii'] * voxelSize
    branchGeom['branch_angles'] = branch_angles * 180 / np.pi
    branchGeom['length'] = branchGeom['length'] * voxelSize
    branchGeom['euclidean length'] = branchGeom['euclidean length'] * voxelSize

    branchGeom['length ratio'] = length_ratio_branch
    branchGeom['diam ratio'] = diam_ratio_branch

    return (geom, branchGeom)


#######################
# Function: Find the Major/Minor ratios of length, diameter and branch angle
# Inputs:  geom - contains elements, and their radii, angles and lengths
#          elem_down - contains the index of the downstream elements at each element
# Outputs: grad - the diameter scaling coefficient
#######################

def major_minor(geom, elem_down):

    # extract data
    radii=geom['radii']
    angles=geom['branch angles']
    length=geom['length']

    # create arrays
    Ne=len(elem_down)

    Minor_angle=-1*np.ones(Ne)
    Major_angle = -1*np.ones(Ne)

    D_Major_Minor = -1 * np.ones(Ne)
    D_min_parent = -1 * np.ones(Ne)
    D_maj_parent = -1 * np.ones(Ne)

    L_Major_Minor = -1 * np.ones(Ne)
    L_min_parent = -1 * np.ones(Ne)
    L_maj_parent = -1 * np.ones(Ne)

    for i in range(0, Ne):
        numDown=elem_down[i, 0]

        if numDown>1: # then this element has multiple children, find minor / major child

            d_min=100000
            d_max=0
            for j in range(1, numDown+1): #look throigh children and find widest & thinnest one
                child=np.int(elem_down[i, j])
                d_child=radii[child]

                if d_child>d_max:
                    d_max=d_child
                    daughter_max=child
                if d_child<d_min:
                    d_min = d_child
                    daughter_min = child

            if daughter_max!=daughter_min: # ensure two distinct daughters

                Minor_angle[i]=angles[daughter_min]
                Major_angle[i]=angles[daughter_max]

                if radii[daughter_min]!=0: # avoid divide by zero errors
                    D_Major_Minor[i]=radii[daughter_max]/radii[daughter_min]
                if radii[i] != 0:
                    D_min_parent[i]=radii[daughter_min]/radii[i]
                    D_maj_parent[i]=radii[daughter_max]/radii[i]

                if length[daughter_min] != 0:
                    L_Major_Minor[i] = length[daughter_max] / length[daughter_min]
                if length[i] != 0:
                    L_min_parent[i] = length[daughter_min] / length[i]
                    L_maj_parent[i] = length[daughter_max] / length[i]
    return {'Minor_angle': Minor_angle, 'Major_angle': Major_angle, 'D_maj_min': D_Major_Minor, 'D_min_P': D_min_parent,'D_maj_P': D_maj_parent, 'L_maj_min': L_Major_Minor, 'L_min_P': L_min_parent,'L_maj_P': L_maj_parent}


#######################
# Function: Find & print vascular depth, span & volume
# Inputs:  nodes - M x 3 array with node coordinates
#          elems - N x 3 array containing elements
#          orders - N x 1 array with order of each element
#          vascVol - vascularVolume, in mm^3 (number of voxels in the volume image)
# Outputs: depth- Vascular depth, mm
#          span, vascular span, mm
#######################

def overall_shape_parameters(nodes, elems,orders, vascVol):

    # find umbilical insertion node
    inds=np.where(orders==np.max(orders))
    inds=inds[0]
    inds=inds[len(inds)-1] #last umbilical elemnt should be further down
    umbNode=np.int(elems[inds,2]) #take second node as elems point downstream
    umbEnd=nodes[umbNode,:]

    # extract termimal nodes
    inds = np.where(orders == 1)
    inds = inds[0]
    endNodes=elems[inds,2] #take second node as elems point downstream
    endNodes=np.squeeze(endNodes.astype(int))
    endPoints=nodes[endNodes,:]

    # Vascular Span
    dists=scipy.spatial.distance.pdist(endPoints, 'euclidean') # pairwise distance between points
    span=(np.max(dists))

    # Get placenta volume
    ET = EllipsoidTool()
    (center, radii, rotation) = ET.getMinVolEllipse(endPoints, .01)
    ellipseVol = ET.getEllipsoidVolume(radii)
    ET.plotEllipsoid(center, radii, rotation, ax=None, plotAxes=False, cageColor='b', cageAlpha=0.2)

    print('\nOverall Placenta Shape')
    print('-----------------------')
    print("Vascular Span = " + str(span) + ' mm')
    print('Vascular Volume = ' + str(vascVol) + ' mm^3' )
    print('Placenta Volume' + str(ellipseVol)+ ' mm^3' )
    print("Vascular Density = " + str(vascVol/ellipseVol))
    print('\n')
    return span



#######################
# Function: Finds diameter scaling coefficent and creates a plot of diameters
# Inputs:  diam - an N x 1 array containing the diameter of all N branches
#          cutoff - the diameter below which will be excluded when calculating the diameter scaling coefficient
# Outputs: grad - the diameter scaling coefficient
#######################

def diam_log_cdf(diam, cutoff):

    plt.figure()

    # Reversed cumulative histogram.
    plt.subplot(1, 2, 1)
    plt.title('Reverse CDF')
    plt.ylabel('Proportion of Segments')
    plt.xlabel('Diameter')
    n_bins = np.int(np.round(len(diam) / 25)) # bin spacing
    n, bins, patches = plt.hist(diam, n_bins, density=True, histtype='step', cumulative=-1, label='Reverse cdf') # cdf

    # Add cutoff.
    plt.plot([cutoff,cutoff], [0,1],label ='cutoff')
    plt.legend()

    # Log Log plot
    plt.subplot(1, 2, 2)
    plt.ylabel('log(Number Segments)')
    plt.xlabel('log(Diameter)')

    # only take bins above cutoff
    bins=bins[0:len(bins)-1] # so same size as n
    n=n[bins>cutoff]
    bins=bins[bins>cutoff]

    # log log plot
    x=np.log(bins)
    yData=np.log(n)
    plt.plot(x,yData , 'k--', linewidth=1.5, label='Data')

    # fit line to data
    xFit=np.unique(x)
    yFit=np.poly1d(np.polyfit(x, yData, 1))(np.unique(x))
    plt.plot(np.unique(x), np.poly1d(np.polyfit(x, yData, 1))(np.unique(x)),label='linear fit')

    # Scaling Coefficient is gradient
    grad=(yFit[len(yFit)-1]-yFit[0])/(xFit[len(xFit)-1]-xFit[0])
    heading=('Diam Coeff = ' + str(grad))
    plt.title(heading)
    plt.legend()
    plt.show()

    # R^2 value
    yMean = [np.mean(yData) for y in yData]
    r2=1 - (sum((yFit - yData) * (yFit - yData))/sum((yMean-yData) * (yMean - yData)))
    print('Diameter Scaling Coefficient = ' + str(grad) + ' Rsquared = ' + str(r2))
    return grad, r2


####
# Function: find statistics on branching tree and display as table, sorting my generations in the tree
# Inputs: geom - contains various element properties (length, radius etc.) by element
#         orders - contains strahler order and generation of each element
# Outputs: table of information according to generation prints to screen
######

def generation_summary_statistics_rachel(geom, orders, major_minor_results):

    # unpack inputs
    generation = orders['generation']

    diam = 2 * geom['radii']
    length = geom['length']
    euclid_length = geom['euclidean length']
    angles = geom['branch angles']

    diam_ratio = geom['diam_ratio']
    length_ratio = geom['length_ratio']

    Minor_angle = major_minor_results['Minor_angle']
    Major_angle = major_minor_results['Major_angle']

    D_Major_Minor = major_minor_results['D_maj_min']
    D_min_parent = major_minor_results['D_min_P']
    D_maj_parent = major_minor_results['D_maj_P']

    L_Major_Minor = major_minor_results['L_maj_min']
    L_min_parent = major_minor_results['L_min_P']
    L_maj_parent = major_minor_results['L_maj_P']

    # statisitcs by generation
    num_gens= int(max(generation))
    values_by_gen = np.zeros([num_gens, 34])

    for n_gen in range(0, num_gens):

        element_list = (generation == n_gen + 1)

        diam_list = np.extract(element_list, diam)
        len_list = np.extract(element_list, length)

        # account for zero diameters
        diam_bool = diam_list > 0
        len_bool = len_list > 0
        list = np.logical_and(diam_bool, len_bool)
        diam_list = diam_list[list]
        len_list = len_list[list]

        # assign stats for each order
        values_by_gen[n_gen, 0] = n_gen + 1  # order
        values_by_gen[n_gen, 1] = len(np.extract(element_list, element_list))  # number of branches

        values_by_gen[n_gen, 2] = np.mean(np.extract(element_list, length))  # length
        values_by_gen[n_gen, 3] = np.std(np.extract(element_list, length))  # length std

        values_by_gen[n_gen, 4] = np.mean(diam_list)  # diameter
        values_by_gen[n_gen, 5] = np.std(diam_list)  # diameter std

        values_by_gen[n_gen, 6] = np.mean(np.extract(element_list, euclid_length))  # euclidean length
        values_by_gen[n_gen, 7] = np.std(np.extract(element_list, euclid_length))  # euclidean length std

        values_by_gen[n_gen, 8] = np.mean(len_list / diam_list)  # length / diameter
        values_by_gen[n_gen, 9] = np.std(len_list / diam_list)  # length / diameter std

        values_by_gen[n_gen, 10] = np.mean(
            np.extract(element_list, length) / np.extract(element_list, euclid_length))  # tortuosity
        values_by_gen[n_gen, 11] = np.std(
            np.extract(element_list, length) / np.extract(element_list, euclid_length))  # tortuosity


        if n_gen > 0:


            angle_list = np.extract(element_list, angles)
            angle_list = angle_list[angle_list > 0]
            if len(angle_list)>0:
                values_by_gen[n_gen, 12] = np.mean(angle_list)  # angles
                values_by_gen[n_gen, 13] = np.std(angle_list)  # angles std

            Minor_angle_list = np.extract(element_list, Minor_angle)
            Minor_angle_list = Minor_angle_list[Minor_angle_list > 0]
            Major_angle_list = np.extract(element_list, Major_angle)
            Major_angle_list = Major_angle_list[Major_angle_list > 0]
            if len(Minor_angle_list) > 0:
                values_by_gen[n_gen, 14] = np.mean(Minor_angle_list)  # minor angles
                values_by_gen[n_gen, 15] = np.std(Minor_angle_list)
                values_by_gen[n_gen, 16] = np.mean(Major_angle_list)  # major angles
                values_by_gen[n_gen, 17] = np.std(Major_angle_list)

            lengthRatio = np.extract(element_list, length_ratio)
            lengthRatio = lengthRatio[lengthRatio > 0]

            L_min_parent_list = np.extract(element_list, L_min_parent)
            L_min_parent_list = L_min_parent_list[L_min_parent_list > 0]
            L_maj_parent_list = np.extract(element_list, L_maj_parent)
            L_maj_parent_list = L_maj_parent_list[L_maj_parent_list > 0]
            L_Major_Minor_list = np.extract(element_list, L_Major_Minor)
            L_Major_Minor_list = L_Major_Minor_list[L_Major_Minor_list > 0]
            if len(L_min_parent_list) > 0:
                values_by_gen[n_gen, 18] = np.mean(lengthRatio)  # len ratio
                values_by_gen[n_gen, 19] = np.std(lengthRatio)  # len ratio
                values_by_gen[n_gen, 20] = np.mean(L_min_parent_list)
                values_by_gen[n_gen, 21] = np.std(L_min_parent_list)
                values_by_gen[n_gen, 22] = np.mean(L_maj_parent_list)
                values_by_gen[n_gen, 23] = np.std(L_maj_parent_list)
                values_by_gen[n_gen, 24] = np.mean(L_Major_Minor_list)
                values_by_gen[n_gen, 25] = np.std(L_Major_Minor_list)

            diamRatio = np.extract(element_list, diam_ratio)
            diamRatio = diamRatio[diamRatio > 0]
            D_min_parent_list = np.extract(element_list, D_min_parent)
            D_min_parent_list = D_min_parent_list[D_min_parent_list > 0]
            D_maj_parent_list = np.extract(element_list, D_maj_parent)
            D_maj_parent_list = D_maj_parent_list[D_maj_parent_list > 0]
            D_Major_Minor_list = np.extract(element_list, D_Major_Minor)
            D_Major_Minor_list = D_Major_Minor_list[D_Major_Minor_list > 0]
            if len(D_min_parent_list) > 0:
                values_by_gen[n_gen, 26] = np.mean(diamRatio)  # diam ratio
                values_by_gen[n_gen, 27] = np.std(diamRatio)  # diam std
                values_by_gen[n_gen, 28] = np.mean(D_min_parent_list)
                values_by_gen[n_gen, 29] = np.std(D_min_parent_list)
                values_by_gen[n_gen, 30] = np.mean(D_maj_parent_list)
                values_by_gen[n_gen, 31] = np.std(D_maj_parent_list)
                values_by_gen[n_gen, 32] = np.mean(D_Major_Minor_list)
                values_by_gen[n_gen, 33] = np.std(D_Major_Minor_list)

    # print table
    header = ['Gen', 'NumBranches', 'Length(mm)', 'std', 'Diameter(mm)', 'std', 'Euclidean Length(mm)', 'std',
              'Len/Diam', 'std', 'Tortuosity', 'std', 'Angles', 'std','Minor Angle','std','Major Angle','std', 'LLparent', 'std', 'LminLparent', 'std', 'LmajLparent', 'std', 'LminLmaj', 'std', 'DDparent', 'std','DminDparent', 'std','DmajDparent', 'std','DminDmaj', 'std']
    print('\n')
    print('Statistics By Generation: ')
    print('..................')
    print(tabulate(values_by_gen, headers=header))

    # statistics independent of order
    values_overall = np.zeros([1, 34])

    element_list = (generation > 0)
    diam_list = np.extract(element_list, diam)

    len_list = np.extract(element_list, length)
    len_list = len_list[diam_list > 0]
    diam_list = diam_list[diam_list > 0]

    angle_list = np.extract(element_list, angles)
    angle_list = angle_list[angle_list > 0]

    Minor_angle_list = np.extract(element_list, Minor_angle)
    Minor_angle_list = Minor_angle_list[Minor_angle_list > 0]
    Major_angle_list = np.extract(element_list, Major_angle)
    Major_angle_list = Major_angle_list[Major_angle_list > 0]
    L_min_parent_list = np.extract(element_list, L_min_parent)
    L_min_parent_list = L_min_parent_list[L_min_parent_list > 0]
    L_maj_parent_list = np.extract(element_list, L_maj_parent)
    L_maj_parent_list = L_maj_parent_list[L_maj_parent_list > 0]
    L_Major_Minor_list = np.extract(element_list, L_Major_Minor)
    L_Major_Minor_list = L_Major_Minor_list[L_Major_Minor_list > 0]
    D_min_parent_list = np.extract(element_list, D_min_parent)
    D_min_parent_list = D_min_parent_list[D_min_parent_list > 0]
    D_maj_parent_list = np.extract(element_list, D_maj_parent)
    D_maj_parent_list = D_maj_parent_list[D_maj_parent_list > 0]
    D_Major_Minor_list = np.extract(element_list, D_Major_Minor)
    D_Major_Minor_list = D_Major_Minor_list[D_Major_Minor_list > 0]

    # assign stats for each order
    values_overall[0, 0] = -1
    values_overall[0, 1] = len(np.extract(element_list, element_list))  # number of branches

    values_overall[0, 2] = np.mean(len_list)  # length
    values_overall[0, 3] = np.std(len_list)  # length std

    values_overall[0, 4] = np.mean(diam_list)  # diameter
    values_overall[0, 5] = np.std(diam_list)  # diameter std

    values_overall[0, 6] = np.mean(np.extract(element_list, euclid_length))  # euclidean length
    values_overall[0, 7] = np.std(np.extract(element_list, euclid_length))  # euclidean length std

    values_overall[0, 8] = np.mean(len_list / diam_list)  # length / diameter
    values_overall[0, 9] = np.std(len_list / diam_list)  # length / diameter std

    values_overall[0, 10] = np.mean(
        np.extract(element_list, length) / np.extract(element_list, euclid_length))  # tortuosity
    values_overall[0, 11] = np.std(
        np.extract(element_list, length) / np.extract(element_list, euclid_length))  # tortuosity

    values_overall[0, 12] = np.mean(angle_list)  # angles
    values_overall[0, 13] = np.std(angle_list)  # angles std
    values_overall[0, 14] = np.mean(Minor_angle_list)  # minor angles
    values_overall[0, 15] = np.std(Minor_angle_list)
    values_overall[0, 16] = np.mean(Major_angle_list)  # major angles
    values_overall[0, 17] = np.std(Major_angle_list)

    lengthRatio = np.extract(element_list, length_ratio)
    lengthRatio = lengthRatio[lengthRatio > 0]
    values_overall[0, 18] = np.mean(lengthRatio)  # len ratio
    values_overall[0, 19] = np.std(lengthRatio)  # len ratio
    values_overall[0, 20] = np.mean(L_min_parent_list)
    values_overall[0, 21] = np.std(L_min_parent_list)
    values_overall[0, 22] = np.mean(L_maj_parent_list)
    values_overall[0, 23] = np.std(L_maj_parent_list)
    values_overall[0, 24] = np.mean(L_Major_Minor_list)
    values_overall[0, 25] = np.std(L_Major_Minor_list)

    diamRatio = np.extract(element_list, diam_ratio)
    diamRatio = diamRatio[diamRatio > 0]
    values_overall[0, 26] = np.mean(diamRatio)  # diam ratio
    values_overall[0, 27] = np.std(diamRatio)  # diam std
    values_overall[0, 28] = np.mean(D_min_parent_list)
    values_overall[0, 29] = np.std(D_min_parent_list)
    values_overall[0, 30] = np.mean(D_maj_parent_list)
    values_overall[0, 31] = np.std(D_maj_parent_list)
    values_overall[0, 32] = np.mean(D_Major_Minor_list)
    values_overall[0, 33] = np.std(D_Major_Minor_list)

    # print table
    header = ['Gen', 'NumBranches', 'Length(mm)', 'std', 'Diameter(mm)', 'std', 'Euclidean Length(mm)', 'std',
              'Len/Diam', 'std', 'Tortuosity', 'std', 'Angles', 'std','Minor Angle','std','Major Angle','std', 'LLparent', 'std', 'LminLparent', 'std', 'LmajLparent', 'std', 'LminLmaj', 'std', 'DDparent', 'std','DminDparent', 'std','DmajDparent', 'std','DminDmaj', 'std']

    print(tabulate(values_overall, headers=header))
    print('\n')

    return np.concatenate((values_by_gen, values_overall),0)

####
# Function: find statistics on branching tree (by Strahler order) and display as table
# Inputs: branchGeom - contains properties (length, radius etc.) by Strahler branch
#         geom - contains various element properties (length, radius etc.) by element
#         orders - contains strahler order and generation of each element
# Outputs: table of information according to order and other information printed to screen
######

def summary_statistics_rachel(branchGeom, geom, orders, major_minor_results):

    # branch inputs
    branchDiam = 2 * branchGeom['radii']
    branchLen = branchGeom['length']
    branchEucLen = branchGeom['euclidean length']
    branchOrder = branchGeom['order']
    branchAngles = branchGeom['branch_angles']
    branchLenRatio = branchGeom['length ratio']
    branchDiamRatio = branchGeom['diam ratio']

    # statisitcs by order
    num_orders = int(max(branchOrder))
    values_by_order = np.zeros([num_orders, 20])

    for n_ord in range(0, num_orders):

        branch_list = (branchOrder == n_ord + 1)

        diam_list = np.extract(branch_list, branchDiam)
        len_list = np.extract(branch_list, branchLen)

        # account for zero diameters
        diam_bool = diam_list > 0
        len_bool = len_list > 0
        list = np.logical_and(diam_bool, len_bool)
        diam_list = diam_list[list]
        len_list = len_list[list]

        # assign stats for each order
        values_by_order[n_ord, 0] = n_ord + 1  # order
        values_by_order[n_ord, 1] = len(np.extract(branch_list, branch_list))  # number of branches

        values_by_order[n_ord, 2] = np.mean(np.extract(branch_list, branchLen))  # length
        values_by_order[n_ord, 3] = np.std(np.extract(branch_list, branchLen))  # length std

        values_by_order[n_ord, 4] = np.mean(diam_list)  # diameter
        values_by_order[n_ord, 5] = np.std(diam_list)  # diameter std

        values_by_order[n_ord, 6] = np.mean(np.extract(branch_list, branchEucLen))  # euclidean length
        values_by_order[n_ord, 7] = np.std(np.extract(branch_list, branchEucLen))  # euclidean length std

        values_by_order[n_ord, 8] = np.mean(len_list / diam_list)  # length / diameter
        values_by_order[n_ord, 9] = np.std(len_list / diam_list)  # length / diameter std

        values_by_order[n_ord, 10] = np.mean(
            np.extract(branch_list, branchLen) / np.extract(branch_list, branchEucLen))  # tortuosity
        values_by_order[n_ord, 11] = np.std(
            np.extract(branch_list, branchLen) / np.extract(branch_list, branchEucLen))  # tortuosity


        if n_ord < num_orders - 1:


            angle_list = np.extract(branch_list, branchAngles)
            angle_list = angle_list[angle_list > 0]

            values_by_order[n_ord, 12] = np.mean(angle_list)  # angles
            values_by_order[n_ord, 13] = np.std(angle_list)  # angles std

            lengthRatio = np.extract(branch_list, branchLenRatio)
            lengthRatio = lengthRatio[lengthRatio > 0]

            values_by_order[n_ord, 14] = np.mean(lengthRatio)  # len ratio
            values_by_order[n_ord, 15] = np.std(lengthRatio)  # len ratio

            diamRatio = np.extract(branch_list, branchDiamRatio)
            diamRatio = diamRatio[diamRatio > 0]

            values_by_order[n_ord, 16] = np.mean(diamRatio)  # diam ratio
            values_by_order[n_ord, 17] = np.std(diamRatio)  # diam std

        values_by_order[n_ord, 18] = values_by_order[n_ord-1, 1]/values_by_order[n_ord, 1]  # Bifurcation ratio
        values_by_order[n_ord, 19] = np.sum(np.square(diam_list)*np.pi/4 ) # Total CSA

    # print table
    header = ['Order', 'NumBranches', 'Length(mm)', 'std', 'Diameter(mm)', 'std', 'Euclidean Length(mm)', 'std',
              'Len/Diam', 'std', 'Tortuosity', 'std', 'Angles', 'std', 'LenRatio', 'std', 'DiamRatio', 'std','Bifurcation Ratio','TotalCSA']
    print('\n')
    print('Statistics By Order: ')
    print('..................')
    print(tabulate(values_by_order, headers=header))

    # statistics independent of order
    values_overall = np.zeros([1, 20])

    branch_list = (branchOrder > 0)
    diam_list = np.extract(branch_list, branchDiam)

    len_list = np.extract(branch_list, branchLen)
    len_list = len_list[diam_list > 0]
    diam_list = diam_list[diam_list > 0]

    angle_list = np.extract(branch_list, branchAngles)
    angle_list = angle_list[angle_list > 0]

    # assign stats for each order
    values_overall[0, 0] = -1
    values_overall[0, 1] = len(np.extract(branch_list, branch_list))  # number of branches

    values_overall[0, 2] = np.mean(len_list)  # length
    values_overall[0, 3] = np.std(len_list)  # length std

    values_overall[0, 4] = np.mean(diam_list)  # diameter
    values_overall[0, 5] = np.std(diam_list)  # diameter std

    values_overall[0, 6] = np.mean(np.extract(branch_list, branchEucLen))  # euclidean length
    values_overall[0, 7] = np.std(np.extract(branch_list, branchEucLen))  # euclidean length std

    values_overall[0, 8] = np.mean(len_list / diam_list)  # length / diameter
    values_overall[0, 9] = np.std(len_list / diam_list)  # length / diameter std

    values_overall[0, 10] = np.mean(
        np.extract(branch_list, branchLen) / np.extract(branch_list, branchEucLen))  # tortuosity
    values_overall[0, 11] = np.std(
        np.extract(branch_list, branchLen) / np.extract(branch_list, branchEucLen))  # tortuosity

    values_overall[0, 12] = np.mean(angle_list)  # angles
    values_overall[0, 13] = np.std(angle_list)  # angles std

    lengthRatio = np.extract(branch_list, branchLenRatio)
    lengthRatio = lengthRatio[lengthRatio > 0]
    values_overall[0, 14] = np.mean(lengthRatio)  # len ratio
    values_overall[0, 15] = np.std(lengthRatio)  # len ratio

    diamRatio = np.extract(branch_list, branchDiamRatio)
    diamRatio = diamRatio[diamRatio > 0]
    values_overall[0, 16] = np.mean(diamRatio)  # diam ratio
    values_overall[0, 17] = np.std(diamRatio)  # diam std

    values_overall[0, 18] = np.mean(values_by_order[1:num_orders, 18])  # Bifurcation ratio



    # print table
    header = ['     ', '           ', '          ', '   ', '            ', '   ', '                     ', '   ',
              '        ', '   ', '           ', '   ', '      ', '   ', '        ', '   ', '         ', '   ','                 ']

    print(tabulate(values_overall, headers=header))
    print('\n')

    # unpack inputs
    strahler = orders['strahler']
    generation = orders['generation']

    diam = 2*geom['radii']
    length = geom['length']
    length2 = length[(diam > 0)]
    diam = diam[(diam > 0)]
    euclid_length = geom['euclidean length']

    angles = geom['branch angles']
    angles = angles[angles > 0]  # get rid of first elem

    diam_ratio = geom['diam_ratio']
    diam_ratio = diam_ratio[(diam_ratio > 0)]

    length_ratio = geom['length_ratio']
    length_ratio = length_ratio[(length_ratio > 0)]

    # unpack inputs
    Minor_angle = major_minor_results['Minor_angle']
    Minor_angle = Minor_angle[Minor_angle > 0]

    Major_angle = major_minor_results['Major_angle']
    Major_angle = Major_angle[Major_angle > 0]

    D_Major_Minor = major_minor_results['D_maj_min']
    D_Major_Minor = D_Major_Minor[D_Major_Minor > 0]

    D_min_parent = major_minor_results['D_min_P']
    D_min_parent = D_min_parent[(D_min_parent > 0)]

    D_maj_parent = major_minor_results['D_maj_P']
    D_maj_parent = D_maj_parent[(D_maj_parent > 0)]

    L_Major_Minor = major_minor_results['L_maj_min']
    L_Major_Minor = L_Major_Minor[L_Major_Minor > 0]

    L_min_parent = major_minor_results['L_min_P']
    L_min_parent = L_min_parent[(L_min_parent > 0)]

    L_maj_parent = major_minor_results['L_maj_P']
    L_maj_parent = L_maj_parent[(L_maj_parent > 0)]

    # Segment statistics
    print('Segment statistics: ')
    print('..................')
    print('Num Segments = ' + str(len(strahler)))
    print('Total length = ' + str(np.sum(branchGeom['length'])) + ' mm')
    print('Num generations = ' + str(max(generation)))
    terminalGen = generation[(strahler == 1)]
    print('Average Terminal generation (std) = ' + str(np.mean(terminalGen)) + ' (' + str(np.std(terminalGen)) + ')')
    print('Segment Tortuosity = ' + str(np.mean(length / euclid_length)) + ' (' + str(
        np.std(length / euclid_length)) + ')')
    print('Average Length (std) = ' + str(np.mean(length)) + ' (' + str(np.std(length)) + ')')
    print('Average Euclidean Length (std) = ' + str(np.mean(euclid_length)) + ' (' + str(np.std(euclid_length)) + ')')
    print('Average Diameter (std) = ' + str(np.mean(diam)) + ' (' + str(np.std(diam)) + ')')
    print('Average L/D (std) = ' + str(np.mean(length2/diam)) + ' (' + str(np.std(length2/diam)) + ')') ########

    print('Segment Angles = ' + str(np.mean(angles)) + ' (' + str(np.std(angles)) + ')')
    print('    Minor Angle = ' + str(np.mean(Minor_angle)) + ' (' + str(np.std(Minor_angle)) + ')')
    print('    Major Angle = ' + str(np.mean(Major_angle)) + ' (' + str(np.std(Major_angle)) + ')')
    print('D/Dparent = ' + str(np.mean(diam_ratio)) + ' (' + str(np.std(diam_ratio)) + ')')
    print('    Dmin/Dparent = ' + str(np.mean(D_min_parent)) + ' (' + str(np.std(D_min_parent)) + ')')
    print('    Dmaj/Dparent = ' + str(np.mean(D_maj_parent)) + ' (' + str(np.std(D_maj_parent)) + ')')
    print('    Dmaj/Dmin = ' + str(np.mean(D_Major_Minor)) + ' (' + str(np.std(D_Major_Minor)) + ')')
    print('L/Lparent = ' + str(np.mean(length_ratio)) + ' (' + str(np.std(length_ratio)) + ')')
    print('    Lmin/Lparent = ' + str(np.mean(L_min_parent)) + ' (' + str(np.std(L_min_parent)) + ')')
    print('    Lmaj/Lparent = ' + str(np.mean(L_maj_parent)) + ' (' + str(np.std(L_maj_parent)) + ')')
    print('    Lmaj/Lmin = ' + str(np.mean(L_Major_Minor)) + ' (' + str(np.std(L_Major_Minor)) + ')')
    print('\n')

    # Find  Strahler Ratios: Rb, Rl, Rd
    Num_Branches = values_by_order[:, 1]
    Diameter_strahler = values_by_order[:, 4]
    Length_strahler = values_by_order[:, 2]
    Orders_strahler = values_by_order[:, 0]

    [Rb, r2] = find_strahler_ratio(Orders_strahler, Num_Branches)
    print('Rb = ' + str(Rb) + ' Rsq = ' + str(r2))
    [Rd, r2] = find_strahler_ratio(Orders_strahler, Diameter_strahler)
    print('Rd = ' + str(Rd) + ' Rsq = ' + str(r2))
    [Rl, r2] = find_strahler_ratio(Orders_strahler, Length_strahler)
    print('Rl = ' + str(Rl) + ' Rsq = ' + str(r2))

    return np.concatenate((values_by_order, values_overall),0)

######
# Function: Finds Strahler ratio of variable
#     Inputs: Orders- an array containing the orders of the vascular tree
#             Factor - an array containing a value for each order of the tree e.g. number of branches at each order
#     Outputs: Strahler ratio e.g. Rb, Rd for that factor and the R^2 value of the linear fit used to produce it
######

def find_strahler_ratio(Orders, Factor):

    x = Orders
    yData = np.log(Factor)
    plt.plot(x, yData, 'k--', linewidth=1.5, label='Data')

    # fit line to data
    xFit = np.unique(Orders)
    yFit = np.poly1d(np.polyfit(x, yData, 1))(np.unique(x))
    plt.plot(np.unique(x), yFit, label='linear fit')

    # Scaling Coefficient is gradient of the fit
    grad = (yFit[len(yFit) - 1] - yFit[0]) / (xFit[len(xFit) - 1] - xFit[0])
    grad=np.abs(grad)
    grad=np.exp(grad)

    # R^2 value
    yMean = [np.mean(yData) for y in yData]
    r2 = 1 - (sum((yFit - yData) * (yFit - yData)) / sum((yMean - yData) * (yMean - yData)))

    heading = ('Strahler Ratio = ' + str(grad))
    plt.title(heading)
    plt.legend()
    plt.show()

    return grad, r2

def export_solution_2(data, groupname, filename, name):
    # Write header
    type = "exelem"
    data_num = len(data)
    filename = filename + '.' + type
    f = open(filename, 'w')
    f.write(" Group name: %s\n" % groupname)
    f.write("Shape. Dimension=1\n")
    f.write("#Scale factor sets=0\n")
    f.write("#Nodes=0\n")
    f.write(" #Fields=1\n")
    f.write("1) " + name + ", field, rectangular cartesian, #Components=1\n")
    f.write(name + ".  l.Lagrange, no modify, grid based.\n")
    f.write(" #xi1=1\n")

    # Write element values
    for x in range(0, data_num):
        f.write(" Element:            %s 0 0\n" % int(x + 1))
        f.write("   Values:\n")
        f.write("          %s" % np.squeeze(data[x]))
        f.write("   %s \n" % np.squeeze(data[x]))
    f.close()

    return

def renumber_elems(elems):
    for i in range(0,len(elems)):
        elems[i,0] = i

    return elems

def find_radius_normal_projection(SkeletonImage, VolumeImage, elems, nodes, euclid_radii):

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
            print('large neighbourhood in', large_neighbourhood)
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

    # Compare radii to euclidean radii
    normal_radii[normal_radii < 0] = euclid_radii[normal_radii<0]
    euclid_radii[euclid_radii == 0] = normal_radii[euclid_radii == 0]
    euclid_radii[euclid_radii == 0]=np.min(euclid_radii[euclid_radii>0]) # so no chance of div 0
    difference = abs(normal_radii - euclid_radii)/ euclid_radii
    cutoff = 0.33333 # distances that are larger than cutoff are not used, 1 means that the distance is the same magnitude as the euclidean distance
    normal_radii[difference > cutoff] = euclid_radii[difference > cutoff]

    return (normal_radii)

######
# Function: Find radius for a single point on an element using normal projects
# Inputs: coord1, coord2 - 3d coordinates of two points on the centreline of the element
#         VolumeImage - a logical matrix with volume image, Nx x Ny x Nz
# Outputs: distances - an M x 1 array of various distances estimated at this slice
#######

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

             step = step + 0.2 # step update by 1/5 of a voxel (could increase in order to speed up)
             counter = counter + 1
             currentPosition = np.double(coord1) + step*normal # take step in direction of normal vector
             currentPosition = np.round(currentPosition)
             currentValue = VolumeImage[int(currentPosition[0]), int(currentPosition[1]), int(currentPosition[2])]
        distances[i] = step - 0.2

    return distances

######
# Function: Find next voxel up the skeleton element
# Inputs: coord - current position on element
#         endCoord - final position on element
#         large neighbourhood - the 5 x 5 x 5 area surrounding coord, a section of skeleton image
# Outputs: nextCoord - next position along element
#######

def find_connected_voxels(large_neighbourhood, coord, endCoord):
    print(large_neighbourhood)
    error=0 # default
    large_neighbourhood[2, 2, 2] = 0  # so can't stay at same spot
    neighbourhood = np.copy(large_neighbourhood[1:4, 1:4, 1:4])  # 3 x 3 x 3 area surrounding the coord

    inds = np.where(neighbourhood == 1) # new places to step to
    numInds=len(inds[0])
    subs = np.squeeze(np.column_stack([inds[0],inds[1],inds[2]]))

    # Case of one choice of where to step
    if numInds == 1:
        nextCoord = [coord[0] + subs[0] - 1, coord[1] + subs[1] - 1, coord[2] + subs[2] - 1] # current voxel is at [1,1,1] in neighbourhood

    # Case of no options of where to go
    else: #Either numInds==0 (nowhere to go in direct neighbour hood) OR numInds>1 (too many places to choose in direct neighbourhood

        large_inds = np.where(large_neighbourhood == 1) # places to jump to
        large_subs=np.squeeze(np.column_stack([large_inds[0],large_inds[1],large_inds[2]]))
        numLargeInds = len(large_inds[0])

        if numLargeInds == 1: # only one place to jump, go there
            nextCoord = [coord[0] + large_subs[0] - 2, coord[1] + large_subs[1] - 2, coord[2] + large_subs[2] - 2] # current voxel is at [2, 2, 2] in large neighbourhood

        elif numLargeInds > 1:  # choose where to jump using distance criteria
            endDist = np.zeros(numLargeInds)
            done = 0

            for i in range(0,numLargeInds):

                nextCoord = [coord[0] + large_subs[i,0] - 2, coord[1] + large_subs[i,1] - 2, coord[2] + large_subs[i,2] - 2]  # current voxel is at [2, 2, 2] in large neighbourhood

                # fast track to end if it is in the large neighbourhood (most of the time this section solves it)
                if (np.prod(nextCoord == endCoord)):
                    done = 1
                    break

                # distance to end coord
                endDist[i] = np.sqrt(np.sum(np.square(nextCoord - endCoord)))

            if done == 0: # choose point closest to end point
                # print('Connectivity Error')
                error = 1  # have failed to get to end coord of this branch
                nextCoord = endCoord  # fast track to end

        else: # can't go anywhere, error
            #print('Connectivity Error')
            error = 1 # have failed to get to end coord of this branch
            nextCoord = endCoord # fast track to end

    next_coord=np.array(nextCoord) #check type is correct

    return (next_coord, error)