#!/usr/bin/env python
import numpy as np
from . import pg_utilities
from . import imports_and_exports
from . import skeleton_to_tree
import sys
from numpy import matlib
from scipy import spatial as sp


"""
.. module:: analyse_tree
  :synopsis: One sentence synopis (brief) could appear in module index.

:synopsis:A longer synopsis that could appear on the home page for that module in documentation.

"""

def analyse_branching(geom,branch_geom,ordering_system,conversionFactor,voxelSize):
    """ Does a branching analysis on the tree defined by 'geom'
     Inputs:
       -  geom:  A geometry structure consisting of element list, node location and radii/lengths
       -  ordering_system: the ordering system to be used in analysis (e.g. 'strahler', 'horsfield'

    Returns: Prints to screen a  table of branching properties (one per generation, one per order) and overall summary statistics

    """

    elem_cnct = pg_utilities.element_connectivity_1D(geom['nodes'], geom['elems'])
    geom['order'] = evaluate_orders(geom['nodes'], geom['elems'])
    print("Full geom has max order", np.max(geom['order']['strahler']))
    geom['elem_up'] = elem_cnct['elem_up']
    geom['elem_down'] = elem_cnct['elem_down']

    elem_cnct_branch = pg_utilities.element_connectivity_1D(geom['nodes'], branch_geom['elems'])
    branch_geom['order'] = evaluate_orders(geom['nodes'], branch_geom['elems'])
    print("branch geom has max order", np.max(branch_geom['order']['strahler']))
    branch_geom['elem_up'] = elem_cnct_branch['elem_up']
    branch_geom['elem_down'] = elem_cnct_branch['elem_down']

    ## Find Results
    branch_geom = branch_properties(geom,branch_geom)
    major_minor_results=major_minor(branch_geom, branch_geom['elem_down']) #major/minor child stuff

    ## tabulate data
    
    generation_table,bs = summary_statistics(branch_geom,  major_minor_results,'generation')
    strahler_table,bs = summary_statistics(branch_geom,  major_minor_results,'strahler')

    
    return geom, branch_geom, generation_table,strahler_table,bs

def arrange_by_branches(geom, elem_up, order,generation):
    """ Finds properties of according to each Branch of the tree, where a branch is a set of elements with the
              same order. Ordering system can be any defined in 'evaluate_ordering'
     Inputs:
        - geom: contains elems, and various element properties (length, radius etc.)
        - elem_up - contains index of upstream elements for each element
        - order:  contains order of each element
        - generation: contains generation of each element
     Outputs:
        branchGeom: contains the properties arrange in arrays according to each branch:
               radius / length / euclidean length / strahler order: all M x 1 arrays where M is number of branches
        branches: an N x 1 array where N is the number of elements, contains branch number of each element
    """

    # find branches, which are branches with the same 'generation' as one another
    num_elems = len(order)
    branches = np.zeros(num_elems,dtype=int)
    branchNum = 0

    for i in range(0, num_elems):
        if generation[i] != generation[elem_up[i, 1]]:  # does not belong with upstream branch
            branchNum = branchNum + 1
        else:
            branchNum = branches[elem_up[i, 1]]
        branches[i] = branchNum

    num_branches = int(max(branches)) +1 #including inlet

    # sort results into branch groups
    lengths = geom['length']
    radii = geom['radii']
    nodes= geom['nodes']
    elems = geom['elems']

    branchRad = np.zeros(num_branches)
    branchLen = np.zeros(num_branches)
    branchEucLen = np.zeros(num_branches)
    branchOrder = -1. * np.ones(num_branches)

    for i in range(0, num_branches):
        branchElements = np.where(branches == i) #find elements belonging to branch number
        branchElements = branchElements[0]

        for j in range(0, len(branchElements)): #go through all elements in branch
            ne = branchElements[j]
            branchOrder[i] = order[ne]
            branchLen[i] = branchLen[i] + lengths[ne]
            branchRad[i] = branchRad[i] + radii[ne]

        branchRad[i] = branchRad[i] / len(branchElements) # to get average radius

        startNode=nodes[int(elems[branchElements[0],1]),:]
        endNode=nodes[int(elems[branchElements[len(branchElements)-1],2]),:]

        branchEucLen[i]=np.sqrt(np.sum(np.square(startNode[1:4]-endNode[1:4])))

    return {'radii': branchRad, 'length': branchLen, 'euclidean length': branchEucLen, 'order': branchOrder,
            'branches': branches}

def arrange_by_strahler_order(geom, find_inlet_loc, inlet_loc):
    """ Rearranges elems (and corresponding properties) according to their strahler order

    Inputs:
       - geom: A geometry structure consisting of element list, node location and radii/lengths
       - inlet_loc: the coordinates of the parent node for the entire tree (if known)
       - find_inlet_loc: boolean variable specifying whether to use inlet location provided (0) or to find the inlet location automatically (1)


    Returns:
       - geom: contains elems and properties, reordered according to strahler order so that no element can be higher in the element list than a higher order branch

    """
    # set up arrays
    nodes = geom['nodes']
    elem_properties = np.column_stack([geom['radii'], geom['length'], geom['euclidean length'], geom['elems']])
    elems = np.copy(geom['elems'])  # as elems is altered in this function
    elems = elems[:, 1:3]  # get rid of first column which contains element numbef
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

def branch_properties(geom,branch_geom):
    num_branch = len(branch_geom['elems'])
    num_elems = len(geom['elems'])
    branch_geom['length'] = np.zeros(num_branch)
    branch_geom['radii']=np.zeros(num_branch)
    branch_geom['branch angles'] = -1.*np.ones(num_branch)# seg_angles * 180 / np.pi
    branch_geom['diam ratio'] = -1.*np.ones(num_branch)#diam_ratio
    branch_geom['length ratio'] = -1.*np.ones(num_branch)#length_ratio
    branch_geom['angles'] = -1.*np.ones(num_branch)# seg_angles * 180 / np.pi
    branch_geom['inlet radius'] = 0.
    my_elems = geom['elems'][:,0]
    for nb in range(0,num_branch):
       tmp_elems = my_elems[geom['branch id']==nb+1]
       branch_geom['length'][nb] = np.sum(geom['length'][tmp_elems])
       branch_geom['radii'][nb] = np.sum(geom['radii'][tmp_elems])/float(len(tmp_elems))
       if nb<2:
           print(geom['length'][tmp_elems],len(tmp_elems))
       
    for nb in range(0,num_branch):
        if branch_geom['elem_up'][nb,0]>0:
           nb_up = branch_geom['elem_up'][nb,1]
           branch_geom['diam ratio'][nb] = branch_geom['radii'][nb]/branch_geom['radii'][nb_up]
           branch_geom['length ratio'][nb] = branch_geom['length'][nb]/branch_geom['length'][nb_up]   
        else: #Umbilical artery  (inlet)
           branch_geom['inlet radius'] = branch_geom['radii'][nb]
        if branch_geom['elem_down'][nb,0]>1:
           nbUp = nb
           endNode = int(branch_geom['elems'][nbUp, 2])
           startNode = int(branch_geom['elems'][nbUp, 1])
           v_parent = geom['nodes'][endNode, 1:4] - geom['nodes'][startNode, 1:4]
           v_parent = v_parent / np.linalg.norm(v_parent)
           for b in range(0,branch_geom['elem_down'][nb,0]):
               nbDown = branch_geom['elem_down'][nb,b+1]
               endNode = int(branch_geom['elems'][nbDown, 2])
               startNode = int(branch_geom['elems'][nbDown, 1])
               v_daughter = geom['nodes'][endNode, 1:4] - geom['nodes'][startNode, 1:4]
               v_daughter = v_daughter / np.linalg.norm(v_parent)
               angle = pg_utilities.angle_two_vectors(v_parent, v_daughter)*180./np.pi
               #print(nb_daughter)
               branch_geom['branch angles'][nbDown]=angle
               #print(branch_geom['branch angles'][nb_daughter])
        

   # kb=0
   # for ne in range(0,num_elems):
   #     if geom['elem_down'][ne,0]>1:
   #        neUp = ne
   #        endNode = int(geom['elems'][neUp, 2])
   #        startNode = int(geom['elems'][neUp, 1])
   #        v_parent = geom['nodes'][endNode, 1:4] - geom['nodes'][startNode, 1:4]
   #        v_parent = v_parent / np.linalg.norm(v_parent)
   #        for b in range(0,geom['elem_down'][ne,0]):
   #            neDown = geom['elem_down'][ne,b+1]
   #            endNode = int(geom['elems'][neDown, 2])
   #            startNode = int(geom['elems'][neDown, 1])
   #            v_daughter = geom['nodes'][endNode, 1:4] - geom['nodes'][startNode, 1:4]
   #            v_daughter = v_daughter / np.linalg.norm(v_parent)
   #            angle = pg_utilities.angle_two_vectors(v_parent, v_daughter)*180./np.pi
   #            nb_daughter = int(geom['branch id'][neDown-1])
   #            #print(nb_daughter)
   #            branch_geom['angles'][nb_daughter]=angle
   #            #print(branch_geom['branch angles'][nb_daughter])
   #     
   #        kb= kb+geom['elem_down'][ne,0]
   #        #print('branch_point',ne,kb,geom['branch id'][ne])
    return branch_geom

def calc_terminal_branch(node_loc, elems):
    """ Generates a list of terminal nodes associated with a branching geometry based on element connectivity.

    Inputs:
       - node_loc: array of coordinates (locations) of nodes of tree branches
       - elems: array of elements showing element connectivity
       
    Returns:
       - terminal_elems: array of terminal element number
       - terminal_nodes: array of terminal nodes
       - total_terminals: total number of terminal branches in the whole tree
            
    A way you might want to use me is:    

    >>> node_loc =np.array([[ 0.,0.,0.,-1.,2.,0.,0.], [1.,0.,0.,-0.5,2.,0.,0.],[2.,0.,-0.5,0.,1.31578947,0.,0.],[3.,0.,0.5,0.,0.,0.,0.]])
    >>> elems = np.array([[0 ,0 ,1], [1 ,1 ,2], [2 ,1 ,3]])
    >>> calc_terminal_branch(node_loc,elems)

    This will return:  

    >>> terminal_elems: [1 2]
    >>> terminal_nodes: [2 3]
    >>> total_terminals: 2
    """
    # This function generates a list of terminal nodes associated with a branching geometry
    # inputs are node locations and elements
    num_elems = len(elems)
    num_nodes = len(node_loc)
    elem_cnct = pg_utilities.element_connectivity_1D(node_loc, elems)

    terminal_branches = np.zeros(num_elems, dtype=int)
    terminal_nodes = np.zeros(num_nodes, dtype=int)

    num_term = 0
    for ne in range(0, num_elems):
        if elem_cnct['elem_down'][ne][0] == 0:  # no downstream element
            terminal_branches[num_term] = ne
            terminal_nodes[num_term] = elems[ne][2]  # node that leaves the terminal element
            num_term = num_term + 1

    terminal_branches = np.resize(terminal_branches, num_term)
    terminal_nodes = np.resize(terminal_nodes, num_term)

    print('Total number of terminals assessed, num_terminals =  ' + str(num_term))

    return {'terminal_elems': terminal_branches, 'terminal_nodes': terminal_nodes, 'total_terminals': num_term}

def cal_br_vol_samp_grid(rectangular_mesh, branch_nodes, branch_elems, branch_radius, volume, thickness, ellipticity,
                         start_elem):
    """ Calculate total volume and diameter of branches in each sampling grid element

    Inputs are:
    - rectangular_mesh: rectangular sampling grid
    - branch_nodes: array of coordinates (locations) of nodes of tree branches
    - branch_elems: array of element showing element connectivity
    - branch_radius: array of branch radius
    - volume: volume of placenta
    - thickness: thickness of placenta
    - ellipticity: ellipticity of placenta
    - start_elem: number of element to start calculating tissue volume

    Return:
    - br_vol_in_grid: array of total tissue volume in each sampling grid element
    - br_diameter_in_grid: array of total diameter*volume in each sampling grid element

    A way you might want to use me is:

    >>> thickness =  2.1  #mm
    >>> ellipticity = 1.00  #no unit
    >>> volume=5    #mm3
    >>> rectangular_mesh = {}
    >>> rectangular_mesh['nodes'] = np.array([[-0.5, -0.5, -1.5],[ 0.5, -0.5,-1.5],[-0.5,  0.5 ,-1.5],[ 0.5 , 0.5, -1.5],[-0.5 ,-0.5, -0.5],[ 0.5 ,-0.5 ,-0.5],[-0.5 , 0.5 ,-0.5],[ 0.5 , 0.5 ,-0.5],[-0.5, -0.5 , 0.5],[ 0.5, -0.5 , 0.5],[-0.5  ,0.5 , 0.5],[ 0.5 , 0.5  ,0.5]])
    >>> rectangular_mesh['elems'] = [[ 0,  0,  1,  2,  3,  4, 5, 6, 7],[1,4,5,6,7,8,9,10,11]]
    >>> rectangular_mesh['total_elems'] = 2
    >>> branch_elems={}
    >>> branch_elems['elems']=[[0 ,0, 1]]
    >>> branch_nodes={}
    >>> branch_nodes['nodes']=np.array([[ 0.,0.,0., -1., 2.,0.,0.],[ 1.,0.,0.,-0.5 ,2.,0.,0.]])
    >>> branch_radius=[0.1]
    >>> start_elem=0
    >>> cal_br_vol_samp_grid(rectangular_mesh,  branch_nodes['nodes'], branch_elems['elems'],branch_radius, volume, thickness,ellipticity, start_elem)

    This will return:

    >>> br_vol_in_grid[0]: 0.01396263
    >>> br_diameter_in_grid[0]: 0.00279253
    """

    # Define the resolution of cylinder for analysis
    num_points_xy = 8
    num_points_z = 8
    # Define information about sampling grid required to place data points in correct locations
    total_sample_elems = rectangular_mesh['total_elems']
    gr = pg_utilities.samp_gr_for_node_loc(rectangular_mesh)
    # Define the placental ellipsoid
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)  # calculate radii of ellipsoid
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    unit_cyl_points = np.zeros((num_points_xy * num_points_xy * num_points_z, 3))
    # Define a cylinder of points of radius 1 and length 1
    x = np.linspace(-1, 1, num_points_xy)
    y = np.linspace(-1, 1, num_points_xy)
    num_accepted = 0
    for k in range(0, num_points_z + 1):
        for i in range(0, num_points_xy):
            for j in range(0, num_points_xy):
                if (x[i] ** 2 + y[j] ** 2) <= 1:
                    new_z = 1 / np.double(num_points_z) * k
                    unit_cyl_points[num_accepted][0] = x[i]
                    unit_cyl_points[num_accepted][1] = y[j]
                    unit_cyl_points[num_accepted][2] = new_z
                    num_accepted = num_accepted + 1
    unit_cyl_points.resize(num_accepted, 3, refcheck=False)
    cyl_points = np.copy(unit_cyl_points)
    cylindervector = np.array([0.0, 0.0, 1.0])

    ###Define and initialise arrays to be populated

    # The volume of each branch
    vol_each_br = np.zeros(len(branch_elems))
    # Array for total volume of sampling grid in each element
    total_vol_samp_gr = np.zeros(total_sample_elems)
    # Array for diameter variable of sampling grid in each element (this variable is to be used for weighted diameter calculation)
    total_diameter_samp_gr = np.zeros(total_sample_elems)
    # initialise counters
    branch_count = 0
    volume_outside_ellipsoid = 0.0
    volume_inside_ellipsoid = 0.0

    for ne in range(start_elem, len(branch_elems)):  # len(branch_elems)):  # looping for all branchs in tree

        node1 = branch_nodes[branch_elems[ne][1]][1:4]  # coor of start node of a branch element
        node2 = branch_nodes[branch_elems[ne][2]][1:4]  # coor of end node of a branch element
        node1in = pg_utilities.check_in_on_ellipsoid(node1[0], node1[1], node1[2], x_radius, y_radius, z_radius)
        node2in = pg_utilities.check_in_on_ellipsoid(node2[0], node2[1], node2[2], x_radius, y_radius, z_radius)

        if not node1in and not node2in:
            print('Warning, element ' + str(ne) + 'is not in ellipsoid, if this is not expected check your geometry')
            print('Skipping this element from analysis')
            continue
        elif not node1in or not node2in:
            print('Warning, element ' + str(ne) + 'has one node not in the ellipsoid.')
            print('The first node ' + str(node1) + ' is ' + str(node1in) + ' (True means inside).')
            print('The second node ' + str(node2) + ' is ' + str(node2in) + ' (True means inside).')
            print('Skipping this element from analysis')
            continue

        branch_vector = node2 - node1
        r = branch_radius[ne]
        length = np.linalg.norm(branch_vector)
        vol_each_br[ne] = np.pi * length * r ** 2.0
        vol_per_point = vol_each_br[ne] / (np.double(num_accepted))

        cyl_points[:, 0:2] = unit_cyl_points[:, 0:2] * r
        cyl_points[:, 2] = unit_cyl_points[:, 2] * length

        desiredvector = branch_vector / np.linalg.norm(branch_vector)

        rotation_axis = np.cross(desiredvector, cylindervector)

        if np.linalg.norm(rotation_axis) == 0:  # aligned
            if node2[2] - node1[2] < 0:
                cyl_points[:, 2] = -1.0 * cyl_points[:, 2]
        else:
            angle = pg_utilities.angle_two_vectors(cylindervector, desiredvector)
            rotation_mat = pg_utilities.rotation_matrix_3d(rotation_axis, angle)
            cyl_points = np.array(np.matrix(cyl_points) * np.matrix(rotation_mat))

        cyl_points[:, 0] = cyl_points[:, 0] + node1[0]
        cyl_points[:, 1] = cyl_points[:, 1] + node1[1]
        cyl_points[:, 2] = cyl_points[:, 2] + node1[2]

        # Array for vol distribution of inidvidual branch (not total)
        vol_distribution_each_br = np.zeros(total_sample_elems, dtype=float)

        for nt in range(0, num_accepted):

            coord_point = cyl_points[nt][0:3]
            inside = pg_utilities.check_in_on_ellipsoid(coord_point[0], coord_point[1], coord_point[2], x_radius,
                                                        y_radius, z_radius)
            if inside:
                nelem = pg_utilities.locate_node(gr[0], gr[1], gr[2], gr[3], gr[4], gr[5], gr[6], gr[7], gr[8],
                                                 coord_point)
                total_vol_samp_gr[nelem] = total_vol_samp_gr[nelem] + vol_per_point
                vol_distribution_each_br[nelem] = vol_distribution_each_br[nelem] + vol_per_point
                volume_inside_ellipsoid = volume_inside_ellipsoid + vol_per_point
            else:
                # Data points lie outside the ellipsoid - this is OK in some cases, so the code shouldn't exit. However,
                # users should be able to check how much is outside of ellipsoid if they believe their branching geometry
                # is set up NOT to go outside the ellipsoid at all.
                volume_outside_ellipsoid = volume_outside_ellipsoid + vol_per_point

        total_diameter_samp_gr = total_diameter_samp_gr + vol_distribution_each_br * 2.0 * r  # this variable is calculated as summation of diameter * vol of branch in grid (to be used for weight_diam)

    percent_outside = volume_outside_ellipsoid / np.sum(total_vol_samp_gr) * 100.0

    total_vol_ml = (volume_outside_ellipsoid + np.sum(total_vol_samp_gr))/1000.0
    sum_branch_ml = np.sum(vol_each_br)/1000.0

    print('Analysis complete ' + str(percent_outside) + '% of analysed points lie outside the ellipsoid.')
    print('Total branch volume analysed ' + str(total_vol_ml) + ' (compared with summed branch vol ' + str(
        sum_branch_ml) + ')')

    return {'br_vol_in_grid': total_vol_samp_gr, 'br_diameter_in_grid': total_diameter_samp_gr}


def conductivity_samp_gr(vol_frac, weighted_diameter, elem_list):
    """Calculate conductivity of sampling grid element where villous branches are located

    Inputs are:
    - vol_frac: tissue volume fraction of sampling grid element
    - weighted_diameter: weighted diameter of sampling grid element
    - elem_list: list of elements to assess

    Return:
    - conductivity: conductivity of sampling grid element where the placental tissue are located
    will be in the same units as the weighted diameter (typically mm)

    A way you might want to use me:

    >>> vol_frac= [0.72401065]
    >>> weighted_diameter=[0.17988357]
    >>> non_empties=[0]
    >>> conductivity_samp_gr(vol_frac,weighted_diameter,non_empties)

    This will return:

    >>> conductivity: 7.20937313e-06"""
    max_cond = 0.52
    conductivity = np.zeros(len(vol_frac))
    for i in range(0, len(elem_list)):
        ne = elem_list[i]
        if vol_frac[ne] != 0.0:
            conductivity[ne] = weighted_diameter[ne] ** 2 * (1 - vol_frac[ne]) ** 3 / (180.0 * vol_frac[ne] ** 2)
        elif vol_frac[ne] == 0.0:  # see mabelles thesis
            conductivity[ne] = max_cond
        if conductivity[ne] > max_cond:
            conductivity[ne] = max_cond

    return conductivity


def ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, num_test_points):
    """ Calculates the placental volume associated with each element in a samplling grid

    Inputs are:
    - rectangular_mesh: the rectangular sampling grid
    - volume: placental volume
    - thickness: placental thickness
    - ellipiticity: placental ellipticity
    - num_test_points: resolution of integration quadrature

    Return:
    - pl_vol_in_grid: array of placental volume in each sampling grid element
    - non_empty_rects: array of sampling grid elements that are occupied by placental tissue

    A way you might want to use me is:

    >>> thickness =  (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0  #mm
    >>> ellipticity = 1.00  #no unit
    >>> spacing = 1.0 #no unit
    >>> volume=1 #mm3
    >>> rectangular_mesh = {}
    >>> rectangular_mesh['nodes'] = [[0., 0., 0.], [ thickness/2.0, 0., 0.],[0., thickness/2.0, 0.],[ thickness/2.0, thickness/2.0, 0.],[0., 0., thickness/2.0], [ thickness/2.0, 0., thickness/2.0],[0., thickness/2.0,thickness/2.0],[ thickness/2.0, thickness/2.0, thickness/2.0]]
    >>> rectangular_mesh['elems'] = [[ 0,  0,  1,  2,  3,  4, 5, 6, 7]]
    >>> rectangular_mesh['total_nodes'] =8
    >>> rectangular_mesh['total_elems'] = 1
    >>> num_test_points=25
    >>> ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, num_test_points)

    This will return:

    >>> pl_vol_in_grid: 0.12485807941
    >>> non_empty_rects: 0
    """
    total_elems = rectangular_mesh['total_elems']
    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    # Initialise the array that defines the volume of placenta in each grid element
    pl_vol_in_grid = np.zeros(total_elems)
    non_empty_loc = np.zeros(total_elems, dtype=int)
    non_empty_count = 0

    for ne in range(0, len(elems)):  # looping through elements
        count_in_range = 0
        nod_in_range = np.zeros(8, dtype=int)
        # define range of x, y , and z in the element
        startx = nodes[elems[ne][1]][0]
        endx = nodes[elems[ne][8]][0]
        starty = nodes[elems[ne][1]][1]
        endy = nodes[elems[ne][8]][1]
        startz = nodes[elems[ne][1]][2]
        endz = nodes[elems[ne][8]][2]
        for nod in range(1, 9):
            check_in_range = pg_utilities.check_in_ellipsoid(nodes[elems[ne][nod]][0], nodes[elems[ne][nod]][1],
                                                             nodes[elems[ne][nod]][2], x_radius, y_radius, z_radius)
            check_on_range = pg_utilities.check_on_ellipsoid(nodes[elems[ne][nod]][0], nodes[elems[ne][nod]][1],
                                                             nodes[elems[ne][nod]][2], x_radius, y_radius, z_radius)
            if check_in_range or check_on_range:
                count_in_range = count_in_range + 1
                nod_in_range[nod - 1] = 1
        if count_in_range == 8:  # if all 8 nodes are inside the ellipsoid
            non_empty_loc[non_empty_count] = ne
            non_empty_count = non_empty_count + 1
            pl_vol_in_grid[ne] = (endx - startx) * (endy - starty) * (
                    endz - startz)  # the placental vol in that samp_grid_el is same as vol of samp_grid_el
        elif count_in_range == 0:  # if all 8 nodes are outside the ellpsiod
            # since this samp_grid_el is completely outside, the placental vol is zero
            pl_vol_in_grid[ne] = 0
        else:  # if some nodes in and some nodes out, the samp_grid_el is at the edge of ellipsoid
            # Use trapezoidal quadrature to caculate the volume under the surface of the ellipsoid in each element
            non_empty_loc[non_empty_count] = ne
            non_empty_count = non_empty_count + 1
            # need to map to positive quadrant
            repeat = False
            if (startz < 0 and endz <= 0):
                # need to project to positive z axis
                startz = abs(nodes[elems[ne][8]][2])
                endz = abs(nodes[elems[ne][1]][2])
            elif (startz < 0 and endz > 0):
                # Need to split into components above and below the axis and sum the two
                startz = 0
                endz = abs(nodes[elems[ne][1]][2])
                startz_2 = 0
                endz_2 = nodes[elems[ne][8]][2]
                repeat = True
            xVector = np.linspace(startx, endx, num_test_points)
            yVector = np.linspace(starty, endy, num_test_points)
            xv, yv = np.meshgrid(xVector, yVector)
            zv = z_radius ** 2 * (1 - (xv / x_radius) ** 2 - (yv / y_radius) ** 2)
            for i in range(num_test_points):
                for j in range(num_test_points):
                    if zv[i, j] <= startz ** 2:
                        zv[i, j] = startz ** 2
                    zv[i, j] = np.sqrt(zv[i, j])
                    if zv[i, j] > endz:
                        zv[i, j] = endz
                    elif zv[i, j] < startz:
                        zv[i, j] = startz
            intermediate = np.zeros(num_test_points)
            for i in range(0, num_test_points):
                intermediate[i] = np.trapz(zv[:, i], xVector)
            Value1 = np.trapz(intermediate, yVector)
            pl_vol_in_grid[ne] = (Value1 - startz * (endx - startx) * (endy - starty))
            if repeat:
                xVector = np.linspace(startx, endx, num_test_points)
                yVector = np.linspace(starty, endy, num_test_points)
                xv, yv = np.meshgrid(xVector, yVector)
                zv = z_radius ** 2 * (1 - (xv / x_radius) ** 2 - (yv / y_radius) ** 2)
                for i in range(num_test_points):
                    for j in range(num_test_points):
                        if zv[i, j] <= startz_2 ** 2:
                            zv[i, j] = startz_2 ** 2
                        zv[i, j] = np.sqrt(zv[i, j])
                        if zv[i, j] > endz_2:
                            zv[i, j] = endz_2
                        elif zv[i, j] < startz_2:
                            zv[i, j] = startz_2
                intermediate = np.zeros(num_test_points)
                for i in range(0, num_test_points):
                    intermediate[i] = np.trapz(zv[:, i], xVector)
                Value1 = np.trapz(intermediate, yVector)
                pl_vol_in_grid[ne] = pl_vol_in_grid[ne] + (Value1 - startz_2 * (endx - startx) * (
                        endy - starty))

    print('Number of Non-empty cells: ' + str(non_empty_count))
    print('Total number of cells: ' + str(total_elems))
    non_empty_loc = np.resize(non_empty_loc, non_empty_count)

    return {'pl_vol_in_grid': pl_vol_in_grid, 'non_empty_rects': non_empty_loc}

def evaluate_orders(node_loc, elems):
    """Calculates generations, Horsfield orders, Strahler orders for a given tree
       Works for diverging trees only, but accounts for more than three elements joining at a node
       Inputs:
          node_loc = array with location of nodes
          elems = array with location of elements
    """
    num_elems = len(elems)
    # Calculate connectivity of elements
    elem_connect = pg_utilities.element_connectivity_1D(node_loc, elems)
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
        if elem_upstream[ne][0] != 0:
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


def define_radius_by_order(node_loc, elems, system, inlet_elem, inlet_radius, radius_ratio):
    """ This function defines radii in a branching tree by 'order' of the vessel

     Inputs are:
     - node_loc: The nodes in the branching tree
     - elems: The elements in the branching tree
     - system: 'strahler','horsfield' or 'generation' to define vessel order
     - inlet_elem: element number that you want to define as having inlet_radius
     - inlet_radius: the radius of your inlet vessel
     - radius ratio: Strahler or Horsfield type ratio, defines the slope of log(order) vs log(radius)

     Returns: 
     -radius of each branch

     A way you might want to use me is: 

     >>> node_loc =np.array([[ 0.,0.,0.,-1.,2.,0.,0.], [1.,0.,0.,-0.5,2.,0.,0.],[2.,0.,-0.5,0.,1.31578947,0.,0.],[3.,0.,0.5,0.,0.,0.,0.]])
     >>> elems = np.array([[0 ,0 ,1], [1 ,1 ,2], [2 ,1 ,3]])
     >>> system='strahler'
     >>> inlet_elem=0
     >>> inlet_radius=0.1
     >>> radius_ratio=1.53
     >>> define_radius_by_order(node_loc, elems, system, inlet_elem, inlet_radius, radius_ratio)

    This will return:

    >> radius: [ 0.1, 0.06535948 , 0.06535948]"""
    num_elems = len(elems)
    radius = np.zeros(num_elems)  # initialise radius array
    # Evaluate orders in the system
    orders = evaluate_orders(node_loc, elems)
    elem_order = orders[system]
    ne = inlet_elem
    n_max_ord = elem_order[ne]
    radius[ne] = inlet_radius

    for ne in range(0, num_elems):
        radius[ne] = 10. ** (np.log10(radius_ratio) * (elem_order[ne] - n_max_ord) + np.log10(inlet_radius))

    return radius


def define_radius_by_order_stem(node_loc, elems, system, filename_stem, inlet_radius, radius_ratio):
    """ This function defines radii in a branching tree by 'order' of the vessel

     Inputs are:
     - node_loc: The nodes in the branching tree
     - elems: The elements in the branching tree
     - system: 'strahler','horsfield' or 'generation' to define vessel order
     - filename_stem: filename that includes list of stem villi location and element number
     - inlet_radius: the radius of your inlet vessel
     - radius ratio: Strahler or Horsfield type ratio, defines the slope of log(order) vs log(radius)

     Returns:
     -radius of each branch

     A way you might want to use me is:

    """

    num_elems = len(elems)
    radius = np.zeros(num_elems)  # initialise radius array
    #define stem elems and connectivity
    stem_elems = imports_and_exports.import_stemxy(filename_stem)['elem']
    elem_cnct = pg_utilities.element_connectivity_1D(node_loc, elems)
    # Evaluate orders in the system
    orders = evaluate_orders(node_loc, elems)
    elem_order = orders[system]
    for stem in range(0,len(stem_elems)):
        #For each stem need to form a list of elems that are children of that element
        ne = stem_elems[stem]
        elem_list = pg_utilities.group_elem_parent(ne,elem_cnct['elem_down'])
        n_max_ord = elem_order[ne]
        radius[ne] = inlet_radius

        for noelem in range(0, len(elem_list)):
            ne = elem_list[noelem]
            radius[ne] = 10. ** (np.log10(radius_ratio) * (elem_order[ne] - n_max_ord) + np.log10(inlet_radius))


    return radius



def define_elem_lengths(node_loc, elems):
    """ This function defines element length in a branching tree

         Inputs are:
         - node_loc: The nodes in the branching tree
         - elems: The elements in the branching tree

         Returns:
         -length of each branch

         A way you might want to use me is:

        """

    num_elems = len(elems)
    # length array
    lengths = np.zeros(num_elems)

    for ne in range(0, num_elems):
        np1 = elems[ne][1]
        np2 = elems[ne][2]
        point1 = node_loc[np1][1:4]
        point2 = node_loc[np2][1:4]
        lengths[ne] = np.linalg.norm(point1 - point2)

    return lengths

def find_branch_angles(geom, orders, elem_connect, branchGeom, voxelSize, conversionFactor):
    """Finds branch angles + L/LParent & D/Dparent and scale all results into desired units and degrees
       Inputs:
        - geom: contains elems, and various element properties (length, radius etc.)
        - orders: contains strahler order and generation of each element
        -  elem_connect: contains upstream and downstream elements for each element
        - branchGeom: contains branch properties (length, radius, etc.)
        - voxelSize: for conversion to mm (must be isotropic)
        - conversionFactor: to scale radii correction, printed in log of ImageJ during MySkeletonizationProcess
       Outputs:
        - geom and branchGeom are altered so all there arrays are in correct units (except nodes, and radii_unscaled, which remain in voxels) ##################
        - seg_angles: angle (radians) at each element junction in the tree assigned to each element according to how it branches from its parent
        - diam_ratio:  ratio of length/diameter of each branch, accounting for multi-segment branches
        - length_ratio:  ratio of parent / child lengths, accounting for multi-segment branches
        - diam_ratio: length_ratio / branch_angles are the same but for whole branches
    """

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
    branchGeom['branch angles'] = branch_angles * 180 / np.pi
    branchGeom['length'] = branchGeom['length'] * voxelSize
    branchGeom['euclidean length'] = branchGeom['euclidean length'] * voxelSize

    branchGeom['length ratio'] = length_ratio_branch
    branchGeom['diam ratio'] = diam_ratio_branch

    return (geom, branchGeom)


def find_parent_node(find_inlet_loc, inlet_loc, nodes, radii, elems, elem_properties):
    """Finds the parent node in array either from given coordinates or by finding the terminal branch with the largest radius
       Inputs:
         - nodes: an list of node coordinates in with structure [node num, coord1,coord2,coord3,...]
         - radii: N x 1 array with radius of each element
         - elems: an Nx2(!!!) array with node indices for start and end of node
         - elem_properties:an N x K array, with each row containing various element properties (radii etc.)
         - inlet_loc : the coordinates of the parent node for the entire tree (if known)
         = find_inlet_loc - a boolean variable specifying whether to use inlet location provided (0) or to find the inlet location automatically (1)
       Returns: elems and elem_properties updates so that inlet element is the first element in the list
    """
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

        inlet_loc = np.squeeze(nodes[maxRadInd, 1:4])
        Nn_root = maxRadInd
    # find root node and element from coordinates provided
    else:
        Nn_root = pg_utilities.is_member(inlet_loc, nodes[:,1:4])
        if (Nn_root == -1):
            print("Warning, root node not located")

    print('Inlet Coordinates:' + str(inlet_loc))

    # find root element
    Ne_place = np.where(elems == Nn_root)
    Ne_root = Ne_place[0]  # only need first index
    if len(Ne_root) > 1:
        print("Warning, root node is associated with multiple elements")
    if len(Ne_root) == 0:
        print("Warning, no root element located")
    Ne_root = Ne_root[0]

    # make root element the first element
    elems = pg_utilities.row_swap_2d(elems, 0, Ne_root)
    elem_properties = pg_utilities.row_swap_2d(elem_properties, 0, Ne_root)

    # get element pointing right way
    if (np.squeeze(Ne_place[1]) != 0):
        elems[0, :] = pg_utilities.row_swap_1d(np.squeeze(elems[0, :]), 1, 0)
        elem_properties[0, 4:6] = pg_utilities.row_swap_1d(np.squeeze(elem_properties[0, 4:6]), 1, 0)

    return (elems, elem_properties)


def generation_summary_statistics(geom, orders, major_minor_results):
    """Calculates statistics on branching tree and display as table, sorting my generations in the tree
     Inputs:
       - geom: contains various element properties (length, radius etc.) by element
       - orders: contains strahler order and generation of each element
     Outputs: table of information according to generation prints to screen
    """
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

    # 'LLparent', 'std', 'LminLparent', 'std', 'LmajLparent', 'std', 'LminLmaj', 'std', 'DDparent', 'std','DminDparent', 'std','DmajDparent', 'std','DminDmaj', 'std']
    print('\n')
    print('Statistics By Generation: ')
    print('..................')
    print(
        '   Gen   |   Num   |    L    |  L(std) |    D    |  D(std) |   LEuc  |LEuc(std)|   L_D   | L_D(std)|   Tort  |Tort(std)|   Ang   | Ang(std)|   Amin  |Amin(std)|   Amaj  |Amaj(std)|')
    for n_gen in range(0, num_gens):
        print (
                    ' %7i | %7i | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f |' % (
            values_by_gen[n_gen, 0], values_by_gen[n_gen, 1],
            values_by_gen[n_gen, 2], values_by_gen[n_gen, 3],
            values_by_gen[n_gen, 4], values_by_gen[n_gen, 5],
            values_by_gen[n_gen, 6], values_by_gen[n_gen, 7],
            values_by_gen[n_gen, 8], values_by_gen[n_gen, 9],
            values_by_gen[n_gen, 10], values_by_gen[n_gen, 11],
            values_by_gen[n_gen, 12], values_by_gen[n_gen, 13],
            values_by_gen[n_gen, 14], values_by_gen[n_gen, 15],
            values_by_gen[n_gen, 16], values_by_gen[n_gen, 17]))

    print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
    print (' OVERALL | %7i | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f |' % (
        values_overall[0, 1], values_overall[0, 2], values_overall[0, 3],
        values_overall[0, 4], values_overall[0, 5],
        values_overall[0, 6], values_overall[0, 7],
        values_overall[0, 8], values_overall[0, 9],
        values_overall[0, 10], values_overall[0, 11],
        values_overall[0, 12], values_overall[0, 13],
        values_overall[0, 14], values_overall[0, 15],
        values_overall[0, 16], values_overall[0, 17]))

    print('..................')

    #   'DDparent', 'std','DminDparent', 'std','DmajDparent', 'std','DminDmaj', 'std']
    print('\n')
    print('Statistics By Generation: ')
    print('..................')
    print(
        '   Gen   |   L_Lp  |L_Lp(std)| Lmin_Lp |   std   | Lmaj_Lp |   std   |Lmin_Lmaj|   std   |   D_Dp  |   std   | Dmin_Dp |   std   | Dmaj_Dp |   std   |Dmin_Dmaj|   std   |')
    for n_gen in range(0, num_gens):
        print (
                ' %7i | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f |' % (
            values_by_gen[n_gen, 0], values_by_gen[n_gen, 18],
            values_by_gen[n_gen, 19], values_by_gen[n_gen, 20],
            values_by_gen[n_gen, 21], values_by_gen[n_gen, 22],
            values_by_gen[n_gen, 23], values_by_gen[n_gen, 24],
            values_by_gen[n_gen, 25], values_by_gen[n_gen, 26],
            values_by_gen[n_gen, 27], values_by_gen[n_gen, 28],
            values_by_gen[n_gen, 29], values_by_gen[n_gen, 30],
            values_by_gen[n_gen, 31], values_by_gen[n_gen, 32],
            values_by_gen[n_gen, 33]))

    print(
        '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
    print (
                ' OVERALL | %7i | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f |' % (
            values_overall[0, 18], values_overall[0, 19], values_overall[0, 20],
            values_overall[0, 21], values_overall[0, 22],
            values_overall[0, 23], values_overall[0, 24],
            values_overall[0, 25], values_overall[0, 26],
            values_overall[0, 27], values_overall[0, 28],
            values_overall[0, 29], values_overall[0, 30],
            values_overall[0, 31], values_overall[0, 32],
            values_overall[0, 33]))

    print('-------------')
    print('     |||||   ')
    print('   \ (   )   ')
    print('    ---|---  ')
    print('       |   \ ')
    print('       |     ')
    print('      / \    ')
    print('     /   \   ')
    print('-------------')


    return np.concatenate((values_by_gen, values_overall),0)

def major_minor(geom, elem_down):
    """
     Find the Major/Minor ratios of length, diameter and branch angle
       Inputs:
       - geom: contains elements, and their radii, angles and lengths
       - elem_down: contains the index of the downstream elements at each element
       Outputs:
       - major and minor angle info for each element
    """

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
                child=int(elem_down[i, j])
                d_child=radii[child]

                if d_child>=d_max:
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


#Unused function, no documentation
#def mapping_fields_from_data(datapoints,rectangular_mesh,field1, field2, field3):
#    data_elems = np.zeros(len(datapoints), dtype=int)
#    data_fields = np.zeros((len(datapoints),3))
#    gr = pg_utilities.samp_gr_for_node_loc(rectangular_mesh)
#    for nt in range(0,len(datapoints)):
#        data_elems[nt] = pg_utilities.locate_node(gr[0], gr[1], gr[2], gr[3], gr[4], gr[5], gr[6], gr[7], gr[8],
#                                               datapoints[nt][:])
#        data_fields[nt,0]= field1[data_elems[nt]]
#        data_fields[nt,1] = field2[data_elems[nt]]
#        data_fields[nt, 2] = field3[data_elems[nt]]
#
#
#    return data_fields




def mapping_mesh_sampl_gr(mesh_node_elems, non_empty_rects, conductivity, porosity, export, exportfile):
    """Map the conductivity and porosity value of mesh node with sampling grid element

      Inputs are:
       - mesg_node_elems: array showing where darcy nodes are located inside the sampling grid
       - non_empty_rects: array of non empty sampling grid element
       - conductiviy: conductivity of non-empty sampling grid element
       - porosity: porosity of non-empty sampling grid element

     Return:
       - mapped_con_por: mapped value of conductivity and porosity of each darcy mesh node"""

    mapped_con_por = np.zeros((len(mesh_node_elems), 3)).astype(object)
    mapped_con_por[:, 0] = mapped_con_por[:, 0].astype(int)

    if (export):
        f = open(exportfile, 'w')

    for el in range(0, len(mesh_node_elems)):
        mapped_con_por[el, 0] = el + 1
        print(non_empty_rects, mesh_node_elems)
        if (np.argwhere(non_empty_rects == mesh_node_elems[el,1])):
            mapped_con_por[el, 1] = conductivity[np.argwhere(non_empty_rects == mesh_node_elems[el][1])][0, 0]
            mapped_con_por[el, 2] = porosity[np.where(non_empty_rects == mesh_node_elems[el][1])][0]
        else:  # node sits right on surface, assume empty
            # print('surface node',mesh_node_elems[el][1])
            mapped_con_por[el, 1] = 0.52
            mapped_con_por[el, 2] = 1.0
        if (export):
            f.write("%s %s %s\n" % (mesh_node_elems[el][0], mapped_con_por[el, 1], mapped_con_por[el, 2]))

    if (export):
        f.close()

    return mapped_con_por

##Unused function, no documentation
#def map_mesh_terminals(mesh_nodes, terminal_nodes, branch_nodes, export, exportfile):
#    node_info = np.zeros((len(mesh_nodes), 2), dtype=int)
#    for nnod in terminal_nodes:
#        min_distance = 10000
#        for i in range(0, len(mesh_nodes)):
#            distance = np.sqrt((mesh_nodes[i][1] - branch_nodes[nnod][1]) ** 2.0 + (
#                        mesh_nodes[i][2] - branch_nodes[nnod][2]) ** 2.0 + (
#                                           mesh_nodes[i][3] - branch_nodes[nnod][3]) ** 2.0)
#            if (distance < min_distance):
#                min_distance = distance
#                close_node = int(mesh_nodes[i][0])
#        node_info[close_node - 1][1] = node_info[close_node - 1][1] + 1
#    if (export):
#        f = open(exportfile, 'w')
#    for i in range(0, len(mesh_nodes)):
#        node_info[i][0] = int(mesh_nodes[i][0])
#        if (export):
#            f.write("%s %s\n" % (node_info[i][0], node_info[i][1]))
#    if (export):
#        f.close()


def node_in_sampling_grid(rectangular_mesh, mesh_node_loc):
    """Locate where the 3D mesh nodes are located inside the sampling grid mesh

     Inputs are:
      - rectangular mesh: rectangular sampling grid mesh
      - mesh_node_loc: node locations of mesh

     Return:
      - mesh_node_elems: array which shows the sampling grid element where the mesh nodes are located

    """
    mesh_node_elems = np.zeros((len(mesh_node_loc), 2), dtype=int)
    gr = pg_utilities.samp_gr_for_node_loc(rectangular_mesh)
    for nt in range(0, len(mesh_node_loc)):
        coord_node = mesh_node_loc[nt][1:4]
        nelem = pg_utilities.locate_node(gr[0], gr[1], gr[2], gr[3], gr[4], gr[5], gr[6], gr[7], gr[8], coord_node)
        mesh_node_elems[nt][0] = int(mesh_node_loc[nt][0])
        mesh_node_elems[nt][1] = nelem  # record what element the darcy node is in
        # print(mesh_node_elems[nt])
    return mesh_node_elems


def porosity(vol_frac):
    """ Calculate porosity

    Input is:
     - vol_frac: volume fraction of element

    Return:
     - porosity: porosity of element

    """
    porosity = np.zeros(len(vol_frac))
    porosity = 1 - vol_frac
    return porosity

def summary_statistics(branchGeom, major_minor_results,ordering_system):

    # branch inputs
    branchDiam = np.multiply(2.,branchGeom['radii'])
    branchLen = branchGeom['length']
    branchEucLen = branchGeom['euclidean length']
    branchOrder = branchGeom['order'][ordering_system]
    branchAngles = branchGeom['branch angles']
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
    header = ['LenRatio', 'std', 'DiamRatio', 'std','Bifurcation Ratio','TotalCSA']
    print('\n')
    print('Statistics By Order: ')
    print(ordering_system)
    print('..................')

    print(
        '   Order |   num   |    L    |  L(std) |    D    |  D(std) |   LEuc  |LEuc(std)|   L_D   | L_D(std)|   Tort  |Tort(std)|   Ang   | Ang(std)|    Rl   | Rl(std) |    Rd   | Rd(std)|   Rb    |   CSA   ')
    for n_ord in range(0, num_orders):
        print(' %7i | %7i | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f |%7.4f | %7.4f | %7.4f |'%(values_by_order[n_ord,0],
            values_by_order[n_ord, 1],values_by_order[n_ord,2],
            values_by_order[n_ord, 3],values_by_order[n_ord,4],
            values_by_order[n_ord, 5], values_by_order[n_ord,6],
            values_by_order[n_ord, 7], values_by_order[n_ord, 8],
            values_by_order[n_ord, 9], values_by_order[n_ord, 10],
            values_by_order[n_ord, 11], values_by_order[n_ord, 12],
            values_by_order[n_ord, 13], values_by_order[n_ord, 14],
            values_by_order[n_ord, 15], values_by_order[n_ord, 16],
            values_by_order[n_ord, 17],values_by_order[n_ord, 18],
            values_by_order[n_ord, 19]))

    print('-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
    #        ' %7i | %7.4f | %7.4f
    #        print(tabulate(values_by_order, headers=header))


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



    print (
                ' OVERALL | %7i | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f |' % (
            values_overall[0, 1], values_overall[0, 2], values_overall[0, 3],
            values_overall[0, 4], values_overall[0, 5],
            values_overall[0, 6], values_overall[0, 7],
            values_overall[0, 8], values_overall[0, 9],
            values_overall[0, 10], values_overall[0, 11],
            values_overall[0, 12], values_overall[0, 13],
            values_overall[0, 14], values_overall[0, 15],
            values_overall[0, 16],values_overall[0, 17],values_overall[0, 18]))

    print('-------------')
    print('     |||||   ')
    print('     (   ) / ')
    print('    ---|---  ')
    print('   /   |     ')
    print('       |     ')
    print('      / \    ')
    print('     /   \   ')
    print('-------------')
    
    ## unpack inputs
    strahler = branchGeom['order']['strahler']
    generation = branchGeom['order']['generation']

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
    
    branchVols = np.multiply(branchGeom['radii'],branchGeom['radii'])*np.pi*branchLen
    if (branchGeom['nodes'][:,3]!=branchGeom['nodes'][0,3]).all():
        hull = sp.ConvexHull(branchGeom['nodes'][:,1:4])
    
        #find maxim
        distance_max = 0.
        for i in range(0,len(hull.vertices)):
          for j in range(0,len(hull.vertices)):
              if i != j:
                 distance = np.linalg.norm(branchGeom['nodes'][hull.vertices[i],1:4]-branchGeom['nodes'][hull.vertices[j],1:4])
                 if distance > distance_max:
                     distance_max = distance
    else:
        hull = sp.ConvexHull(branchGeom['nodes'][:, 1:3])
        #hull.volume = 0.
        # find maxim
        distance_max = 0.
        for i in range(0, len(hull.vertices)):
            for j in range(0, len(hull.vertices)):
                if i != j:
                    distance = np.linalg.norm(
                        branchGeom['nodes'][hull.vertices[i], 1:4] - branchGeom['nodes'][hull.vertices[j], 1:4])
                    if distance > distance_max:
                        distance_max = distance

    branch_statistics = np.zeros((48,1))



    # Segment statistics
    print('Branch statistics: ')
    print('..................')
    print('Num Segments = ' + str(len(branchGeom['elems'])))
    branch_statistics[0]=len(branchGeom['elems'])

    print('Total length = ' + str(np.sum(branchGeom['length'])) + ' mm')
    branch_statistics[1]=np.sum(branchGeom['length'])
    print('Total volume of  vessels identified = ' + str(np.sum(branchVols)) + ' mm3')
    branch_statistics[2] =np.sum(branchVols)
    print('Total volume of complex hull representing vessel volume = ' + str(hull.volume) + ' mm3')
    branch_statistics[3] = hull.volume
    print('Max vascular span = ' + str(distance_max) + ' mm')
    branch_statistics[4] = distance_max
    print('Inlet diameter = ' + str(branchGeom['inlet radius']*2.) + ' mm')
    branch_statistics[5] = branchGeom['inlet radius']*2.
    print('Num generations = ' + str(max(branchGeom['order']['generation'])))
    branch_statistics[6]=max(branchGeom['order']['generation'])
    print('Num Strahler Orders = ' + str(max(branchGeom['order']['strahler'])))
    branch_statistics[7] =max(branchGeom['order']['strahler'])

    terminalGen = generation[(strahler == 1)]
    print('Average Terminal generation (std) = ' + str(np.mean(terminalGen)) + ' (' + str(np.std(terminalGen)) + ')')   
    branch_statistics[8] = np.mean(terminalGen)
    branch_statistics[9] = np.std(terminalGen)
    print('Branch Tortuosity = ' + str(np.mean(branchLen / branchEucLen)) + ' (' + str(
        np.std(branchLen / branchEucLen)) + ')')
    branch_statistics[10] = np.mean(branchLen / branchEucLen)
    branch_statistics[11] = np.std(branchLen / branchEucLen)
    print('Average Length mm (std) = ' + str(np.mean(branchLen)) + ' (' + str(np.std(branchLen)) + ')')
    branch_statistics[12] = np.mean(branchLen)
    branch_statistics[13] =np.std(branchLen)
    print('Average Euclidean Length mm (std) = ' + str(np.mean(branchEucLen)) + ' (' + str(np.std(branchEucLen)) + ')')
    branch_statistics[14] = np.mean(branchEucLen)
    branch_statistics[15] = np.std(branchEucLen)
    print('Average Diameter (std) = ' + str(np.mean(branchDiam)) + ' (' + str(np.std(branchDiam)) + ')')
    branch_statistics[16] =np.mean(branchDiam)
    branch_statistics[17] = np.std(branchDiam)
    print('Average L/D (std) = ' + str(np.mean(branchLen/branchDiam)) + ' (' + str(np.std(branchLen/branchDiam)) + ')') ########
    branch_statistics[18] =np.mean(branchLen/branchDiam)
    branch_statistics[19] =np.std(branchLen/branchDiam)
    print('Branch Angles = ' + str(np.mean(branchAngles)) + ' (' + str(np.std(branchAngles)) + ')')
    branch_statistics[20] = np.mean(branchAngles)
    branch_statistics[21] = np.std(branchAngles)
    print('    Minor Angle = ' + str(np.mean(Minor_angle)) + ' (' + str(np.std(Minor_angle)) + ')')
    branch_statistics[22] = np.mean(Minor_angle)
    branch_statistics[23] = np.std(Minor_angle)
    print('    Major Angle = ' + str(np.mean(Major_angle)) + ' (' + str(np.std(Major_angle)) + ')')
    branch_statistics[24] = np.mean(Major_angle)
    branch_statistics[25] = np.std(Major_angle)
    print('D/Dparent = ' + str(np.mean(branchDiamRatio)) + ' (' + str(np.std(branchDiamRatio)) + ')')
    branch_statistics[26] = np.mean(branchDiamRatio)
    branch_statistics[27] = np.std(branchDiamRatio)
    print('    Dmin/Dparent = ' + str(np.mean(D_min_parent)) + ' (' + str(np.std(D_min_parent)) + ')')
    branch_statistics[28] = np.mean(D_min_parent)
    branch_statistics[29] = np.std(D_min_parent)
    print('    Dmaj/Dparent = ' + str(np.mean(D_maj_parent)) + ' (' + str(np.std(D_maj_parent)) + ')')
    branch_statistics[30] = np.mean(D_maj_parent)
    branch_statistics[31] = np.std(D_maj_parent)
    print('    Dmaj/Dmin = ' + str(np.mean(D_Major_Minor)) + ' (' + str(np.std(D_Major_Minor)) + ')')
    branch_statistics[32] = np.mean(D_Major_Minor)
    branch_statistics[33] = np.std(D_Major_Minor)
    print('L/Lparent = ' + str(np.mean(branchLenRatio)) + ' (' + str(np.std(branchLenRatio)) + ')')
    branch_statistics[34] = np.mean(branchLenRatio)
    branch_statistics[35] = np.std(branchLenRatio)
    print('    Lmin/Lparent = ' + str(np.mean(L_min_parent)) + ' (' + str(np.std(L_min_parent)) + ')')
    branch_statistics[36] = np.mean(L_min_parent)
    branch_statistics[37] = np.std(L_min_parent)
    print('    Lmaj/Lparent = ' + str(np.mean(L_maj_parent)) + ' (' + str(np.std(L_maj_parent)) + ')')
    branch_statistics[38] = np.mean(L_maj_parent)
    branch_statistics[39] = np.std(L_maj_parent)
    print('    Lmaj/Lmin = ' + str(np.mean(L_Major_Minor)) + ' (' + str(np.std(L_Major_Minor)) + ')')
    branch_statistics[40] = np.mean(L_Major_Minor)
    branch_statistics[41] = np.std(L_Major_Minor)
    print('\n')


    if  ordering_system == 'strahler':
        # Find  Strahler Ratios: Rb, Rl, Rd
        Num_Branches = values_by_order[0:num_orders, 1]
        Diameter_strahler = values_by_order[0:num_orders, 4]
        Length_strahler = values_by_order[0:num_orders, 2]
        Orders_strahler = values_by_order[0:num_orders, 0]

        print('Branching/length/diameter ratios: ')
        print('..................................')

        [Rb, r2] = pg_utilities.find_strahler_ratio(Orders_strahler, Num_Branches)
        print('Rb = ' + str(Rb) + ' Rsq = ' + str(r2))
        branch_statistics[42] = Rb
        branch_statistics[43] = r2
        [Rd, r2] = pg_utilities.find_strahler_ratio(Orders_strahler[Diameter_strahler>0], Diameter_strahler[Diameter_strahler>0])
        print('Rd = ' + str(Rd) + ' Rsq = ' + str(r2))
        branch_statistics[44] = Rd
        branch_statistics[45] = r2
        [Rl, r2] = pg_utilities.find_strahler_ratio(Orders_strahler[Length_strahler>0], Length_strahler[Length_strahler>0])
        print('Rl = ' + str(Rl) + ' Rsq = ' + str(r2))
        branch_statistics[46] = Rl
        branch_statistics[47] = r2

    print('-------------')
    print('     |||||   ')
    print('   \ (   ) / ')
    print('    ---|---  ')
    print('       |     ')
    print('       |     ')
    print('      / \    ')
    print('     /   \   ')
    print('-------------')
    
    branch_statistics = np.transpose(branch_statistics)

    return np.concatenate((values_by_order, values_overall),0),branch_statistics


def smooth_on_sg(rectangular_mesh, non_empties, field):
    node_field = np.zeros((rectangular_mesh['total_nodes'], 2))

    for i in range(0, len(non_empties)):
        ne = non_empties[i]
        for j in range(1, 9):
            nnod = rectangular_mesh['elems'][ne][j]
            node_field[nnod][0] = node_field[nnod][0] + field[ne]
            node_field[nnod][1] = node_field[nnod][1] + 1.0

    for i in range(0, rectangular_mesh['total_nodes']):
        if (node_field[i][1] != 0.0):
            node_field[i][0] = node_field[i][0] / node_field[i][1]

    for i in range(0, len(non_empties)):
        ne = non_empties[i]
        elem_field = 0.0
        for j in range(1, 9):
            nnod = rectangular_mesh['elems'][ne][j]
            elem_field = elem_field + node_field[nnod][0]
        elem_field = elem_field / 8.0
        field[ne] = elem_field

    return field


def terminal_villous_diameter(num_int_gens, num_convolutes, len_int, rad_int, len_convolute, rad_convolute):
    """ The concept to calculate terminal villous diameter follows the same as terminal_villous_volume calculation.
    Multiply vol of each branch with diameter of each branch and summation of them to be able to calculate the weighted_diameter in the next subroutine

    Inputs:
       - num_int_gens: Number of generations of intermediate villous per terminal 'stem' villus
       - num_convolutes: Number of terminal convolutes per intermediate villous
       - len_int: Length of a typical intermediate villous
       - rad_int: Radius of a typical intermediate villous
       - len_convolute: Length of a typical terminal convolute
       - rad_convolute: Radius of a typical terminal convolute

    Return:
    - term_vill_diameter: diameter value of terminal conduits

    A way you might want to use me is:

    >>> num_int_gens = 3
    >>> num_convolutes = 10
    >>> len_int = 1.5 #mm
    >>> rad_int = 0.03 #mm
    >>> len_convolute = 3.0 #mm
    >>> rad_convolute = 0.025 #mm

    This will return:

    >>> term_vill_diameter: 0.09
    """
    num_ints = 1
    term_vill_diameter = 0.0
    for i in range(0, num_int_gens + 2):
        num_ints = num_ints * 2.0
        diameter_ints = num_ints * (np.pi * len_int * rad_int ** 2.0) * 2 * rad_int
        if i > 0:
            diameter_convolutes = num_ints * num_convolutes * (
                    np.pi * len_convolute * rad_convolute ** 2.0) * 2 * rad_convolute
        else:
            diameter_convolutes = 0.0
        term_vill_diameter = term_vill_diameter + diameter_ints + diameter_convolutes

    return term_vill_diameter


def terminals_in_sampling_grid_fast(rectangular_mesh, terminal_list, node_loc):
    """ Counts the number of terminals in a sampling grid element, will only work with
    rectangular mesh created as in generate_shapes.gen_rectangular_mesh

    Inputs are:
     - Rectangular mesh: the rectangular sampling grid
     - terminal_list: a list of terminal branch
     - node_loc: array of coordinates (locations) of nodes of tree branches

    Return:
     - terminals_in_grid: array showing how many terminal branches are in each sampling grid element
     - terminal_elems: array showing the number of sampling grid element where terminal branches are located

    A way you might want to use me is:

    >>> terminal_list={}
    >>> terminal_list['terminal_nodes']=[3]
    >>> terminal_list['total_terminals']=1
    >>> rectangular_mesh = {}
    >>> rectangular_mesh['nodes'] =np.array( [[ 0.,  0.,  0.],[ 1.,  0. , 0.],[ 0.,  1. , 0.],[ 1. , 1. , 0.],[ 0.,  0. , 1.],[ 1.,  0. , 1.],[ 0. , 1. , 1.],[ 1. , 1. , 1.]])
    >>> rectangular_mesh['elems']=[[0, 0, 1, 2, 3, 4, 5, 6, 7]]
    >>> node_loc =np.array([[ 0.,0.,0.,-1.,2.,0.,0.], [1.,0.,0.,-0.5,2.,0.,0.],[2.,0.,-0.5,0.,1.31578947,0.,0.],[3.,0.,0.5,0.,0.,0.,0.]])   
    >>> terminals_in_sampling_grid_fast(rectangular_mesh, terminal_list, node_loc)

    This will return:

    >>> terminals_in_grid: 1
    >>> terminal_elems: 0
    """
    num_terminals = terminal_list['total_terminals']
    terminals_in_grid = np.zeros(len(rectangular_mesh['elems']), dtype=int)
    terminal_elems = np.zeros(num_terminals, dtype=int)
    gr = pg_utilities.samp_gr_for_node_loc(rectangular_mesh)

    for nt in range(0, num_terminals):
        coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
        nelem = pg_utilities.locate_node(gr[0], gr[1], gr[2], gr[3], gr[4], gr[5], gr[6], gr[7], gr[8], coord_terminal)
        terminals_in_grid[nelem] = terminals_in_grid[nelem] + 1
        terminal_elems[nt] = nelem  # record what element the terminal is in
    return {'terminals_in_grid': terminals_in_grid, 'terminal_elems': terminal_elems}



def terminals_in_sampling_grid(rectangular_mesh, placenta_list, terminal_list, node_loc):
    """ Counts the number of terminals in a sampling grid element for general mesh

    Inputs are:
     - Rectangular mesh: the rectangular sampling grid
     - placenta_list: array of sampling grid element that are located inside the ellipsoid
     - terminal_list: a list of terminal branch
     - node_loc: array of coordinates (locations) of nodes of tree branches

    Return:
     - terminals_in_grid: array showing how many terminal branches are in each sampling grid element
     - terminal_elems: array showing the number of sampling grid element where terminal branches are located
    
    A way you might want to use me is:

    >>> terminal_list = {}
    >>> terminal_list['terminal_nodes'] = [3]
    >>> terminal_list['total_terminals'] = 1
    >>> placenta_list = [7]
    >>> rectangular_mesh = {}
    >>> rectangular_mesh['elems'] = np.zeros((8, 9), dtype=int)
    >>> rectangular_mesh['nodes'] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.], [0., 0., 1.], [1., 0., 1.],
                                     [0., 1., 1.], [1., 1., 1.]]
    >>> rectangular_mesh['elems'][7] = [0, 0, 1, 2, 3, 4, 5, 6, 7]
    >>> node_loc =np.array([[ 0.,0.,0.,-1.,2.,0.,0.], [1.,0.,0.,-0.5,2.,0.,0.],[2.,0.,-0.5,0.,1.31578947,0.,0.],[3.,0.,0.5,0.,0.,0.,0.]])  
    >>> terminals_in_sampling_grid(rectangular_mesh, placenta_list, terminal_list, node_loc)

    This will return:

    >>> terminals_in_grid[7]: 1
    >>> terminal_elems[0]: 7
    """
    num_sample_elems = len(placenta_list)
    num_terminals = terminal_list['total_terminals']
    terminals_in_grid = np.zeros(len(rectangular_mesh['elems']), dtype=int)
    terminal_mapped = np.zeros(num_terminals, dtype=int)
    terminal_elems = np.zeros(num_terminals, dtype=int)

    for ne_i in range(0, num_sample_elems):
        # First node has min x,y,z and last node has max x,y,z
        ne = placenta_list[ne_i]
        if placenta_list[ne_i] > 0:  # There is some placenta in this element (assuming none in el 0)
            first_node = rectangular_mesh['elems'][ne][1]
            last_node = rectangular_mesh['elems'][ne][8]
            min_coords = rectangular_mesh['nodes'][first_node][0:3]
            max_coords = rectangular_mesh['nodes'][last_node][0:3]
            for nt in range(0, num_terminals):
                if terminal_mapped[nt] == 0:
                    in_element = False
                    coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
                    if coord_terminal[0] >= min_coords[0]:
                        if coord_terminal[0] < max_coords[0]:
                            if coord_terminal[1] >= min_coords[1]:
                                if coord_terminal[1] < max_coords[1]:
                                    if coord_terminal[2] >= min_coords[2]:
                                        if coord_terminal[2] < max_coords[2]:
                                            in_element = True
                    if in_element:
                        terminals_in_grid[ne] = terminals_in_grid[ne] + 1
                        terminal_mapped[nt] = 1
                        terminal_elems[nt] = ne

    return {'terminals_in_grid': terminals_in_grid, 'terminal_elems': terminal_elems}



def terminal_volume_to_grid(rectangular_mesh, terminal_list, node_loc, volume, thickness, ellipticity, term_total_vol,
                            term_tissue_vol, term_tissue_diam):
    """ Calculates the volume of terminal unit associated with each sampling grid element

    Inputs are:
     - Rectangular mesh: the rectangular sampling grid
     - terminal_list: a list of terminal branch
     - node_loc: array of coordinates (locations) of nodes of tree branches 
     - volume: volume of placental ellipsoid
     - thickness: thickness of placental ellipsoid
     - ellipticity: ellipticity of placental ellipsoid
     - term_total_vol: total volume of terminal villus
     - term_tissue_vol: volume of terminal conduits
     - term_tissue_diameter: weighted diameter of terminal conduits

    Return:

     - term_vol_in_grid
     - term_diameter_in_grid

    A way you might want to use me is:

      >>> node_loc =np.array([[ 0.,0.,0.,-1.,2.,0.,0.], [1.,0.,0.,-0.5,2.,0.,0.],[2.,0.,-0.5,0.,1.31578947,0.,0.],[3.,0.,0.5,0.,0.,0.,0.]])
      >>> rectangular_mesh = {}
      >>> rectangular_mesh['nodes'] = np.array([[-2., -2. ,-2.],[ 0. ,-2. ,-2.],[ 2. ,-2. ,-2.],[-2. , 0., -2.],[ 0. , 0. ,-2.],[ 2. , 0. ,-2.],[-2. ,-2. , 0.],[ 0. ,-2. , 0.],[ 2. ,-2. , 0.],[-2. , 0. ,0.],[ 0. , 0.,  0.],[ 2.,  0. , 0.],[-2. ,-2. , 2.],[ 0. ,-2. , 2.],[ 2., -2.,  2.],[-2. , 0. , 2.],[ 0.,  0. , 2.],[ 2. , 0.,  2.]])
      >>> rectangular_mesh['elems'] = [[ 0,0,1,3,4,6,7,9,10],[ 1,  1,2,4,5,7,8,10,11],[2,6,7,9,10,12,13,15,16],[4,7,8,10,11,13,14,16,17]]
      >>> rectangular_mesh['total_elems'] = 4
      >>> terminal_list={}
      >>> terminal_list['total_terminals']=1
      >>> terminal_list['terminal_nodes']=[2]
      >>> volume=5
      >>> thickness=2.1
      >>> ellipticity=1
      >>> term_total_vol=0.04#artificial value to match with a smaller ellipsoid
      >>> term_tissue_vol=1.77657064561
      >>> term_tissue_diam=0.090100877305
      >>> terminal_volume_to_grid(rectangular_mesh, terminal_list, node_loc,volume, thickness, ellipticity,term_total_vol, term_tissue_vol, term_tissue_diam)
     
    This will return:
      >>> term_vol_in_grid[0]: 0.44414266
      >>> term_diameter_in_grid[0]: 0.08003529"""

    # Define the resolution of block for analysis
    num_points_xyz = 8
    # number of terminals to assess
    num_terminals = terminal_list['total_terminals']

    # Define information about sampling grid required to place data points in correct locations
    total_sample_elems = rectangular_mesh['total_elems']
    gr = pg_utilities.samp_gr_for_node_loc(rectangular_mesh)
    # Array for total volume  and diameter of sampling grid in each element
    total_vol_samp_gr = np.zeros(total_sample_elems)
    total_diameter_samp_gr = np.zeros(total_sample_elems)

    # Define the placental ellipsoid
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)  # calculate radii of ellipsoid
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    term_vol_points = np.zeros((num_points_xyz * num_points_xyz * num_points_xyz, 3))
    # Define a cylinder of points of radius 1 and length 1 #ARC is this a cube
    x = np.linspace(-1, 1, num_points_xyz)
    y = np.linspace(-1, 1, num_points_xyz)
    zlist = np.linspace(-1, 1, num_points_xyz)
    num_accepted = 0
    for k in range(0, num_points_xyz):
        for i in range(0, num_points_xyz):
            for j in range(0, num_points_xyz):
                new_z = zlist[k]
                term_vol_points[num_accepted][0] = x[i]
                term_vol_points[num_accepted][1] = y[j]
                term_vol_points[num_accepted][2] = new_z
                num_accepted = num_accepted + 1
    term_vol_points.resize(num_accepted, 3, refcheck=False)
    term_vol_points = term_vol_points * term_total_vol ** (1.0 / 3.0)
    vol_per_point = term_tissue_vol / (num_points_xyz * num_points_xyz * num_points_xyz)
    total_volume = 0.0
    for nt in range(0, num_terminals):
        coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
        local_term_points = np.copy(term_vol_points)

        local_term_points[:, 0] = local_term_points[:, 0] + coord_terminal[0]
        local_term_points[:, 1] = local_term_points[:, 1] + coord_terminal[1]
        local_term_points[:, 2] = local_term_points[:, 2] + coord_terminal[2]

        # Array for vol distribution of inidvidual branch (not total)
        vol_distribution_each_br = np.zeros(total_sample_elems, dtype=float)

        for npoint in range(0, num_accepted):

            coord_point = local_term_points[npoint][0:3]
            inside = pg_utilities.check_in_on_ellipsoid(coord_point[0], coord_point[1], coord_point[2], x_radius,
                                                        y_radius, z_radius)
            if inside:
                nelem = pg_utilities.locate_node(gr[0], gr[1], gr[2], gr[3], gr[4], gr[5], gr[6], gr[7], gr[8],
                                                 coord_point)
                total_vol_samp_gr[nelem] = total_vol_samp_gr[nelem] + vol_per_point
                total_volume = total_volume + vol_per_point
                vol_distribution_each_br[nelem] = vol_distribution_each_br[nelem] + vol_per_point

        total_diameter_samp_gr = total_diameter_samp_gr + vol_distribution_each_br * 2 * term_tissue_diam

    return {'term_vol_in_grid': total_vol_samp_gr, 'term_diameter_in_grid': total_diameter_samp_gr}



def terminal_villous_volume(num_int_gens, num_convolutes, len_int, rad_int, len_convolute, rad_convolute,
                            smallest_radius):
    """ This function calculates the average volume of a terminal villous based on structural
    characteristics measured in the literature.

    Inputs:
       - num_int_gens: Number of generations of intermediate villous per terminal 'stem' villois
       - num_convolutes: Number of terminal convolutes per intermediate villous
       - len_int: Length of a typical intermediate villous
       - rad_int: Radius of a typical intermediate villous
       - len_convolute: Length of a typical terminal convolute
       - rad_convolute: Radius of a typical terminal convolute
       - smallest_radius: Minimum radius of a branch in your villoous tree

    Returns:
       - term_vill_volume: Typical volume of a terminal villous


    A way you might want to use me is:

    >>> num_int_gens = 3
    >>> num_convolutes = 10
    >>> len_int = 1.5 #mm
    >>> rad_int = 0.03 #mm
    >>> len_convolute = 3.0 #mm
    >>> rad_convolute = 0.025 #mm
    >>> smallest radius = 0.03 mm
    >>> terminal_villous_volume(num_int_gens,num_convolutes,len_int,rad_int,len_convulute,rad_convolute,smallest_radius)

    This will take the normal average data from Leiser et al (1990, IBBN:3805554680) and calculate
    average volume of terminal villi to be ~1.77 mm^3
    """

    # Each terminal stem villous branches to two immature intermediate villi
    # and then to three generations of mature intermediate villi each with ~10 terminal conduits
    num_ints = 1
    term_vill_volume = 0.0
    term_vill_diameter = 0.0
    sum_vol_ints = 0.0
    sum_vol_conv = 0.0
    radius_step = (smallest_radius - rad_int)/(num_int_gens)
    for i in range(0, num_int_gens + 2):
        num_ints = num_ints * 2.0
        vol_ints = num_ints * np.pi * len_int * (smallest_radius - i*radius_step) ** 2.0
        diameter_ints = vol_ints * 2 * rad_int
        sum_vol_ints = sum_vol_ints + vol_ints
        if i > 0:
            vol_convolutes = num_ints * num_convolutes * np.pi * len_convolute * rad_convolute ** 2.0
            diameter_convolutes = vol_convolutes* 2 * rad_convolute
            sum_vol_conv = sum_vol_conv + vol_convolutes
        else:
            vol_convolutes = 0.0
            diameter_convolutes = 0.0
        term_vill_volume = term_vill_volume + vol_ints + vol_convolutes
        term_vill_diameter = term_vill_diameter + diameter_ints + diameter_convolutes

    proportion_terminal = sum_vol_conv/(sum_vol_conv+sum_vol_ints)

    return {'volume': term_vill_volume, 'diameter': term_vill_diameter, 'propterm': proportion_terminal}


def tissue_vol_in_samp_gr(term_vol_in_grid, br_vol_in_grid):
    """Calculate the total tissue volume (i.e. including terminal conduits) of tree branches

       Inputs are:
       - term_vol_in_grid:total volume of terminal conduits
       - br_vol_in_grid:total volume of branches before terminal conduits

       Return:
       tissue_vol: total tissue volume of whole tree
      
       A way you might want to use me is:

       >>> term_vol_in_grid=0.444
       >>> br_vol_in_grid=0.008

       This will return:
       >>> tissue_vol: 0.452"""

    tissue_vol = br_vol_in_grid + term_vol_in_grid

    return tissue_vol


def vol_frac_in_samp_gr(tissue_vol, sampling_grid_vol,max_allowed,min_allowed):
    """Calculate volume fraction of sampling grid mesh where the villous branches are located

       Inputs are: 
       - tissue_vol: tissue volume in each sampling grid
       - sampling_grid_vol:volume of sampling grid element where placental tissue are located

       Return:
       - vol_frac: volume fraction of sampling grid element where the placental tissue are located

       A way you might want to use me:

       >>> tissue_vol=[0.453]
       >>> sampling_grid_vol={}
       >>> sampling_grid_vol['non_empty_rects']=[0]
       >>> sampling_grid_vol['pl_vol_in_grid']=[0.625]
       vol_frac_in_samp_gr(tissue_vol,sampling_grid_vol)
      
       This will return:
       >>> vol_frac: 0.7248"""

    volumes = sampling_grid_vol['pl_vol_in_grid']
    non_empties = sampling_grid_vol['non_empty_rects']
    vol_frac = np.zeros(len(volumes))

    for i in range(0, len(non_empties)):
        ne = non_empties[i]
        vol_frac[ne] = tissue_vol[ne] / volumes[ne]
        if vol_frac[ne] > max_allowed:
            vol_frac[ne] = max_allowed
        elif vol_frac[ne] <min_allowed:
            vol_frac[ne] = min_allowed

    return vol_frac



def weighted_diameter_in_samp_gr(term_diameter_in_grid, br_diameter_in_grid, tissue_vol):
    """ Calculated weighted_diameter. 
    Weighted_diameter each sampling grid = (d1*v1+d2*v2+d3*v3+...+dn*vn)/(v1+v2+v2+...+vn)

    Inputs are:
    - term_vill_diameter: diameter of terminal conduits
    - br_diameter_in_grid: diameter of branches in each sampling grid
    - terminals_in_grid: number of terminals in each sampling grid
    - tissue_vol: tissue volume in each sampling grid
     
    Return:
    - weighted_diameter: weighted diameter of each sampling grid element
    """

    tissue_diameter = br_diameter_in_grid + term_diameter_in_grid
    np.seterr(divide='ignore', invalid='ignore')
    weighted_diameter = np.nan_to_num(tissue_diameter / tissue_vol)

    return weighted_diameter





