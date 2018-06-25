#!/usr/bin/env python
import numpy as np
from . import pg_utilities
import sys
from numpy import matlib

"""
.. module:: analyse_tree
  :synopsis: One sentence synopis (brief) could appear in module index.

:synopsis:A longer synopsis that could appear on the home page for that module in documentation.

"""


def calc_terminal_branch(node_loc, elems):
    """ Generates a list of terminal nodes associated with a branching geometry

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


def evaluate_orders(node_loc, elems):
    # calculates generations, Horsfield orders, Strahler orders for a given tree
    # Works for diverging trees only
    # Inputs are:
    # node_loc = array with location of nodes
    # elems = array with location of elements
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


def tree_statistics(node_loc, elems, radius, orders):
    # Caclulates tree statistics for a given tree and prints to terminal
    # Inputs are:
    # node_loc: The nodes in the branching tree
    # elems: The elements in the branching tree
    # radius: per element radius
    # orders: per element order

    num_elems = len(elems)
    diameters = 2.0 * radius
    connectivity = pg_utilities.element_connectivity_1D(node_loc, elems)
    elem_upstream = connectivity['elem_up']
    elem_downstream = connectivity['elem_down']
    num_schemes = 3
    generation = orders['generation']
    horsfield = orders['horsfield']
    strahler = orders['strahler']

    index = np.zeros(3, dtype=int)

    # local arrays
    # length array
    lengths = np.zeros(num_elems)
    # ratios: index i is order, index j is ratios of branching (j=1), length (j=2), and diameter (j=3)

    nbranches = np.zeros((num_schemes + 1, num_elems))

    # j = 0 length, j = 1 diameter, j = 4 L/D
    branches = np.zeros((5, num_elems))

    for ne in range(0, num_elems):
        np1 = elems[ne][1]
        np2 = elems[ne][2]
        point1 = node_loc[np1][1:4]
        point2 = node_loc[np2][1:4]
        lengths[ne] = np.linalg.norm(point1 - point2)

    ntotal = 0
    num_dpp = 0
    num_llp = 0
    N = 1
    for ne in range(0, num_elems):
        num_upstream = elem_upstream[ne][0]
        ne0 = elem_upstream[ne][1]
        index[0] = generation[ne]
        index[1] = horsfield[ne]
        index[2] = strahler[ne]
        add = False
        if (num_upstream == 0):  # nothing upstream
            add = True
        elif generation[ne0] != index[0]:  # something upstream and not the
            add = True
        if add:
            N = N + 1
            for i in range(0, num_schemes):
                nbranches[i][N - 1] = index[i]
            if num_upstream != 0:
                nbranches[num_schemes][N - 1] = strahler[ne0]  # strahler order of parent
            else:
                nbranches[num_schemes][N - 1] = 0

        # Add length of all segments along branch, and calculate mean diameter
        n_segments = 1

        mean_diameter = diameter[ne]
        branches[0][N - 1] = lengths[ne]
        ne_next = ne
        while elem_downstream[ne_next][0] == 1:
            ne_next = elem_downstream[ne_next][1]
            branches[0][N - 1] = branches[0][N - 1] + lengths[ne_next]
            mean_diameter = mean_diameter + diameters[ne_next]
            n_segments = n_segments + 1


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
    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']
    startx = np.min(nodes[:, 0])
    xside = nodes[elems[0][8]][0] - nodes[elems[0][1]][0]
    endx = np.max(nodes[:, 0])
    nelem_x = (endx - startx) / xside
    starty = np.min(nodes[:, 1])
    yside = nodes[elems[0][8]][1] - nodes[elems[0][1]][1]
    endy = np.max(nodes[:, 1])
    nelem_y = (endy - starty) / yside
    startz = np.min(nodes[:, 2])
    zside = nodes[elems[0][8]][2] - nodes[elems[0][1]][2]
    endz = np.max(nodes[:, 2])
    nelem_z = (endz - startz) / zside

    for nt in range(0, num_terminals):
        coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
        xelem_num = np.floor((coord_terminal[0] - startx) / xside)
        yelem_num = np.floor((coord_terminal[1] - starty) / yside)
        zelem_num = np.floor((coord_terminal[2] - startz) / zside)
        nelem = int(xelem_num + (yelem_num) * nelem_x + (zelem_num) * (nelem_x * nelem_y))
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


def terminal_volume_to_grid(rectangular_mesh, terminal_list, node_loc,volume, thickness, ellipticity,term_total_vol, term_tissue_vol, term_tissue_diam):
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
    #number of terminals to assess
    num_terminals = terminal_list['total_terminals']

    # Define information about sampling grid required to place data points in correct locations
    total_sample_elems = rectangular_mesh['total_elems']
    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']
    startx = np.min(nodes[:, 0])
    xside = nodes[elems[0][8]][0] - nodes[elems[0][1]][0]
    endx = np.max(nodes[:, 0])
    nelem_x = (endx - startx) / xside
    starty = np.min(nodes[:, 1])
    yside = nodes[elems[0][8]][1] - nodes[elems[0][1]][1]
    endy = np.max(nodes[:, 1])
    nelem_y = (endy - starty) / yside
    startz = np.min(nodes[:, 2])
    zside = nodes[elems[0][8]][2] - nodes[elems[0][1]][2]
    endz = np.max(nodes[:, 2])

    # Array for total volume  and diameter of sampling grid in each element
    total_vol_samp_gr = np.zeros(total_sample_elems)
    total_diameter_samp_gr = np.zeros(total_sample_elems)

    # Define the placental ellipsoid
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)  # calculate radii of ellipsoid
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    term_vol_points = np.zeros((num_points_xyz * num_points_xyz * num_points_xyz, 3))
    # Define a cylinder of points of radius 1 and length 1
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
    term_vol_points = term_vol_points*term_total_vol**(1.0/3.0)
    vol_per_point = term_tissue_vol/(num_points_xyz * num_points_xyz *num_points_xyz)
    total_volume = 0.0
    for nt in range(0, num_terminals):
        coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
        local_term_points = np.copy(term_vol_points)

        local_term_points[:, 0] = local_term_points[:, 0] + coord_terminal[0]
        local_term_points[:, 1] = local_term_points[:, 1] + coord_terminal[1]
        local_term_points[:, 2] = local_term_points[:, 2] + coord_terminal[2]

        # Array for vol distribution of inidvidual branch (not total)
        vol_distribution_each_br=np.zeros(total_sample_elems, dtype=float)

        for npoint in range(0, num_accepted):

            coord_point = local_term_points[npoint][0:3]
            inside = pg_utilities.check_in_on_ellipsoid(coord_point[0], coord_point[1], coord_point[2], x_radius,
                                                        y_radius, z_radius)
            if inside:
                xelem_num = np.floor((coord_point[0] - startx) / xside)
                yelem_num = np.floor((coord_point[1] - starty) / yside)
                zelem_num = np.floor((coord_point[2] - startz) / zside)
                nelem = int(xelem_num + (yelem_num) * nelem_x + (zelem_num) * (nelem_x * nelem_y))
                total_vol_samp_gr[nelem] = total_vol_samp_gr[nelem] + vol_per_point
                total_volume = total_volume + vol_per_point
                vol_distribution_each_br[nelem] = vol_distribution_each_br[nelem] + vol_per_point

        total_diameter_samp_gr = total_diameter_samp_gr + vol_distribution_each_br * 2 * term_tissue_diam

    return {'term_vol_in_grid': total_vol_samp_gr,'term_diameter_in_grid':total_diameter_samp_gr}


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


def cal_br_vol_samp_grid(rectangular_mesh, branch_nodes, branch_elems,branch_radius, volume, thickness, ellipticity, start_elem):
    """ Calculate total volume and diameter of branches in each samp_grid_el 

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

    #Define the resolution of cylinder for analysis
    num_points_xy = 8
    num_points_z = 8
    # Define information about sampling grid required to place data points in correct locations
    total_sample_elems = rectangular_mesh['total_elems']
    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']
    startx = np.min(nodes[:, 0])
    xside = nodes[elems[0][8]][0] - nodes[elems[0][1]][0]
    endx = np.max(nodes[:, 0])
    nelem_x = (endx - startx) / xside
    starty = np.min(nodes[:, 1])
    yside = nodes[elems[0][8]][1] - nodes[elems[0][1]][1]
    endy = np.max(nodes[:, 1])
    nelem_y = (endy - starty) / yside
    startz = np.min(nodes[:, 2])
    zside = nodes[elems[0][8]][2] - nodes[elems[0][1]][2]
    endz = np.max(nodes[:, 2])
    
    #Define the placental ellipsoid
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)  # calculate radii of ellipsoid
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    
    unit_cyl_points = np.zeros((num_points_xy*num_points_xy*num_points_z,3))
    #Define a cylinder of points of radius 1 and length 1
    x = np.linspace(-1,1,num_points_xy)
    y = np.linspace(-1,1,num_points_xy)
    num_accepted = 0
    for k in range(0,num_points_z+1):
        for i in range(0,num_points_xy):
            for j in range(0,num_points_xy):
                if(x[i]**2 + y[j]**2)<=1:
                    new_z = 1 / np.double(num_points_z) * k
                    unit_cyl_points[num_accepted][0] = x[i]
                    unit_cyl_points[num_accepted][1] = y[j]
                    unit_cyl_points[num_accepted][2] = new_z
                    num_accepted = num_accepted+1
    unit_cyl_points.resize(num_accepted, 3, refcheck=False)
    cyl_points = np.copy(unit_cyl_points)
    cylindervector = np.array([0.0,0.0,1.0])
    
    ###Define and initialise arrays to be populated

    #The volume of each branch
    vol_each_br = np.zeros(len(branch_elems))
    # Array for total volume of sampling grid in each element
    total_vol_samp_gr = np.zeros(total_sample_elems)
    # Array for diameter variable of sampling grid in each element (this variable is to be used for weighted diameter calculation)
    total_diameter_samp_gr = np.zeros(total_sample_elems)
    #initialise counters
    branch_count = 0
    volume_outside_ellipsoid = 0.0
    volume_inside_ellipsoid = 0.0

    for ne in range(start_elem,len(branch_elems)):#len(branch_elems)):  # looping for all branchs in tree
        
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
            print('The first node ' + str(node1) + ' is ' + srt(node1in) + ' (True means inside).')
            print('The second node ' + str(node2) + ' is ' + srt(node2in) + ' (True means inside).')
            print('Skipping this element from analysis')
            continue

        branch_vector = node2-node1
        r = branch_radius[ne]
        length = np.linalg.norm(branch_vector)
        vol_each_br[ne] = np.pi*length*r**2.0
        vol_per_point = vol_each_br[ne]/(np.double(num_accepted))
        
        cyl_points[:,0:2]=unit_cyl_points[:,0:2]*r
        cyl_points[:,2] = unit_cyl_points[:,2]*length

        desiredvector = branch_vector / np.linalg.norm(branch_vector)

        rotation_axis = np.cross(desiredvector,cylindervector)

        if np.linalg.norm(rotation_axis) == 0:#aligned
            if node2[2]-node1[2] < 0:
                cyl_points[:, 2] = -1.0*cyl_points[:, 2]
        else:
            angle = pg_utilities.angle_two_vectors(cylindervector,desiredvector)
            rotation_mat = pg_utilities.rotation_matrix_3d(rotation_axis,angle)
            cyl_points = np.array(np.matrix(cyl_points)*np.matrix(rotation_mat))

        cyl_points[:, 0] = cyl_points[:, 0] + node1[0]
        cyl_points[:, 1] = cyl_points[:, 1] + node1[1]
        cyl_points[:, 2] = cyl_points[:, 2] + node1[2]

        # Array for vol distribution of inidvidual branch (not total)
        vol_distribution_each_br=np.zeros(total_sample_elems, dtype=float)
               
        for nt in range(0, num_accepted):
        
            coord_point = cyl_points[nt][0:3]
            inside=pg_utilities.check_in_on_ellipsoid(coord_point[0], coord_point[1], coord_point[2], x_radius, y_radius, z_radius)
            if inside:
                xelem_num = np.floor((coord_point[0] - startx) / xside)
                yelem_num = np.floor((coord_point[1] - starty) / yside)
                zelem_num = np.floor((coord_point[2] - startz) / zside)
                nelem = int(xelem_num + (yelem_num) * nelem_x + (zelem_num) * (nelem_x * nelem_y))
                total_vol_samp_gr[nelem] = total_vol_samp_gr[nelem] + vol_per_point
                vol_distribution_each_br[nelem]=vol_distribution_each_br[nelem]+vol_per_point
                volume_inside_ellipsoid = volume_inside_ellipsoid + vol_per_point
            else:
                #Data points lie outside the ellipsoid - this is OK in some cases, so the code shouldn't exit. However,
                #users should be able to check how much is outside of ellipsoid if they believe their branching geometry
                #is set up NOT to go outside the ellipsoid at all.
                volume_outside_ellipsoid = volume_outside_ellipsoid + vol_per_point
        
        total_diameter_samp_gr = total_diameter_samp_gr + vol_distribution_each_br*2*r#this variable is calculated as summation of diameter * vol of branch in grid (to be used for weight_diam)
    
    percent_outside = volume_outside_ellipsoid/np.sum(total_vol_samp_gr)*100.0
    
    print('Analysis complete ' + str(percent_outside) + '% of analysed points lie outside the ellipsoid.')
    print('Total branch volume analysed ' + str(volume_outside_ellipsoid + np.sum(total_vol_samp_gr)) + ' (' + str(np.sum(vol_each_br)) +')')
    
    return {'br_vol_in_grid': total_vol_samp_gr,'br_diameter_in_grid':total_diameter_samp_gr}



def terminal_villous_volume(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute):
    """ This function calculates the average volume of a terminal villous based on structural
    characteristics measured in the literature.

    Inputs:
       - num_int_gens: Number of generations of intermediate villous per terminal 'stem' villois
       - num_convolutes: Number of terminal convolutes per intermediate villous
       - len_int: Length of a typical intermediate villous
       - rad_int: Radius of a typical intermediate villous
       - len_convolute: Length of a typical terminal convolute
       - rad_convolute: Radius of a typical terminal convolute

    Returns:
       - term_vill_volume: Typical volume of a terminal villous


    A way you might want to use me is:

    >>> num_int_gens = 3
    >>> num_convolutes = 10
    >>> len_int = 1.5 #mm
    >>> rad_int = 0.03 #mm
    >>> len_convolute = 3.0 #mm
    >>> rad_convolute = 0.025 #mm
    >>> terminal_villous_volume(num_int_gens,num_convolutes,len_int,rad_int,len_convulute,rad_convolute)

    This will take the normal average data from Leiser et al (1990, IBBN:3805554680) and calculate
    average volume of terminal villi to be ~1.77 mm^3
    """

    #Each terminal stem villous branches to two immature intermediate villi
    #and then to three generations of mature intermediate villi each with ~10 terminal conduits
    num_ints = 1
    term_vill_volume = 0.0
    for i in range(0,4):
        num_ints = num_ints*2.0
        vol_ints = num_ints*np.pi*len_int*rad_int**2.0
        if i > 0:
            vol_convolutes = num_ints*num_convolutes*np.pi*len_convolute*rad_convolute**2.0
        else:
            vol_convolutes = 0.0
        term_vill_volume = term_vill_volume + vol_ints + vol_convolutes
    
    return term_vill_volume

def tissue_vol_in_samp_gr(term_vol_in_grid,br_vol_in_grid):
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


def vol_frac_in_samp_gr(tissue_vol,sampling_grid_vol):

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

    for i in range(0,len(non_empties)):
        ne = non_empties[i]
        vol_frac[ne] = tissue_vol[ne]/volumes[ne]
        if vol_frac[ne] > 1.0:
            vol_frac[ne] = 1.0


    return vol_frac


def conductivity_samp_gr(vol_frac,weighted_diameter,non_empties):
    """Calculate conductivity of sampling grid element where villous branches are located

    Inputs are: 
    - vol_frac: tissue volume fraction of sampling grid element
    - weighted_diameter: weighted diameter of sampling grid element
    - non_empties: volume of sampling grid element where placenta tissue are located

    Return:
    - conductivity: conductivity of sampling grid element where the placental tissue are located

    A way you might want to use me:

    >>> vol_frac= [0.72401065]
    >>> weighted_diameter=[0.17988357]
    >>> non_empties=[0]
    >>> conductivity_samp_gr(vol_frac,weighted_diameter,non_empties)

    This will return:

    >>> conductivity: 7.20937313e-06"""
    max_cond = 0.52
    conductivity = np.zeros(len(vol_frac))
    for i in range(0,len(non_empties)):
        ne = non_empties[i]
        if vol_frac[ne] != 0.0:
            conductivity[ne] = weighted_diameter[ne]**2*(1-vol_frac[ne])**3/(180.0*vol_frac[ne]**2)
        elif vol_frac[ne] == 0.0:#see mabelles thesis
            conductivity[ne] = max_cond
        if conductivity[ne] > max_cond:
            conductivity[ne] = max_cond

    return conductivity

def terminal_villous_diameter(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute):
  
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
    for i in range(0,4):
        num_ints = num_ints*2.0
        diameter_ints = num_ints*(np.pi*len_int*rad_int**2.0)*2*rad_int
        if i > 0:
            diameter_convolutes = num_ints*num_convolutes*(np.pi*len_convolute*rad_convolute**2.0)*2*rad_convolute
        else:
            diameter_convolutes = 0.0
        term_vill_diameter = term_vill_diameter + diameter_ints + diameter_convolutes
      
    return term_vill_diameter


def weighted_diameter_in_samp_gr(term_diameter_in_grid,br_diameter_in_grid,tissue_vol):
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
    weighted_diameter = np.nan_to_num(tissue_diameter/tissue_vol)

    return weighted_diameter


def porosity(vol_frac):
    """ Calculate porosity

    Input is: 
     - vol_frac: volume fraction of element

    Return: 
     - porosity: porosity of element
    
    """
    porosity= np.zeros(len(vol_frac))
    porosity=1-vol_frac
    return porosity
