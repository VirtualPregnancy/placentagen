#!/usr/bin/env python
import numpy as np
from scipy.spatial import Delaunay
from . import pg_utilities

"""
.. module:: generate_shapes
  :synopsis: Contains code to generate placental shapes for generic placental models.

:synopsis:Contains code to generate placental shapes for generic placental models \n
 (i.e. from literature measures without specific data from an individual
 
"""

def equispaced_data_in_ellipsoid(n, volume, thickness, ellipticity):
    """ Generates equally spaced data points in an ellipsoid.

    Inputs:
       - n: number of data points which we aim to generate
       - volume: volume of ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions

    Returns:
       - Edata: A nx3 array of datapoints, with each point being defined by its x-,y-, and z- coordinates

    A way you might want to use me is:

    >>> n = 100
    >>> volume = 10
    >>> thickness = 3
    >>> ellipticity = 1.1
    >>> equispaced_data_in_ellipsoid(n, volume, thickness, ellipticity)

   This will return 100 data points in an ellipse with z-axis thickness 3, volume 10, and with the y-axis dimension 1.1 times the x-axis dimension.

    """
    # Generates equally spaced data points in an ellipsoid with the following inputs
    # n=number of data points which we aim to generate
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
    data_spacing = (volume / n) ** (1.0 / 3.0)
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    # Aiming to generate seed points that fill a cuboid encompasing the placental volume then remove seed points that
    # are external to the ellipsoid

    num_data = 0  # zero the total number of data points

    # Calculate the number of points that should lie in each dimension in a cube
    nd_x = np.floor(2.0 * (x_radius + data_spacing) / data_spacing)
    nd_y = np.floor(2.0 * (y_radius + data_spacing) / data_spacing)
    nd_z = np.floor(2.0 * (z_radius + data_spacing) / data_spacing)
    nd_x = int(nd_x)
    nd_y = int(nd_y)
    nd_z = int(nd_z)
    # Set up edge node coordinates
    x_coord = np.linspace(-x_radius - data_spacing / 2.0, x_radius + data_spacing / 2.0, nd_x)
    y_coord = np.linspace(-y_radius - data_spacing / 2.0, y_radius + data_spacing / 2.0, nd_y)
    z_coord = np.linspace(-z_radius - data_spacing / 2.0, z_radius + data_spacing / 2.0, nd_z)

    # Use these vectors to form a unifromly spaced grid
    data_coords = np.vstack(np.meshgrid(x_coord, y_coord, z_coord)).reshape(3, -1).T

    # Store nodes that lie within ellipsoid
    Edata = np.zeros((nd_x * nd_y * nd_z, 3))
    for i in range(len(data_coords)):  # Loop through grid
        coord_check = pg_utilities.check_in_ellipsoid(data_coords[i][0], data_coords[i][1], data_coords[i][2], x_radius,
                                                      y_radius, z_radius)

        if coord_check is True:  # Has to be strictly in the ellipsoid
            Edata[num_data, :] = data_coords[i, :]  # add to data array
            num_data = num_data + 1
    Edata.resize(num_data, 3)  # resize data array to correct size

    print('Data points within ellipsoid allocated. Total = ' + str(len(Edata)))

    return Edata


def uniform_data_on_ellipsoid(n, volume, thickness, ellipticity, random_seed):
    """ Generates equally spaced data points on the positive z-surface of an ellipsoid

    Inputs:
       - n: number of data points which we aim to generate
       - volume: volume of ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions

    Returns:
       - chorion_data: A nx3 array of datapoints, with each point being defined by its x-,y-, and z- coordinates

    A way you might want to use me is:

    >>> n = 100
    >>> volume = 10
    >>> thickness = 3
    >>> ellipticity = 1.1
    >>> equispaced_data_on_ellipsoid(n, volume, thickness, ellipticity)

   This will return 100 data points on the positive z-surface ellipse with z-axis thickness 3, volume 10, and with the y-axis dimension 1.1 times the x-axis dimension.

    """
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    area_estimate = np.pi * x_radius * y_radius
    data_spacing = 0.85 * np.sqrt(area_estimate / n)

    chorion_data = np.zeros((n, 3))
    np.random.seed(random_seed)
    generated_seed = 0
    acceptable_attempts = n * 1000  # try not to have too many failures
    attempts = 0

    while generated_seed < n and attempts < acceptable_attempts:
        # generate random x-y coordinates between negative and positive radii
        new_x = np.random.uniform(-x_radius, x_radius)
        new_y = np.random.uniform(-y_radius, y_radius)
        # check if new coordinate is on the ellipse
        if ((new_x / x_radius) ** 2 + (new_y / y_radius) ** 2) < 1:  # on the surface
            if generated_seed == 0:
                generated_seed = generated_seed + 1
                new_z = pg_utilities.z_from_xy(new_x, new_y, x_radius, y_radius, z_radius)
                chorion_data[generated_seed - 1][:] = [new_x, new_y, new_z]
            else:
                reject = False
                for j in range(0, generated_seed + 1):
                    distance = (chorion_data[j - 1][0] - new_x) ** 2 + (chorion_data[j - 1][1] - new_y) ** 2
                    distance = np.sqrt(distance)
                    if distance <= data_spacing:
                        reject = True
                        break
                if reject is False:
                    generated_seed = generated_seed + 1
                    new_z = pg_utilities.z_from_xy(new_x, new_y, x_radius, y_radius, z_radius)
                    chorion_data[generated_seed - 1][:] = [new_x, new_y, new_z]

        attempts = attempts + 1
    chorion_data.resize(generated_seed, 3)  # resize data array to correct size
    print('Data points on ellipsoid allocated. Total = ' + str(len(chorion_data)))

    return chorion_data


def gen_rect_cover_ellipsoid(volume, thickness, ellipticity, x_spacing, y_spacing, z_spacing):
    # Generates equally spaced data nodes and elements and constructs a rectangular 'mesh' that covers the space that is
    # made up of an ellipsoidal placenta
    # volume=volume of ellipsoid
    # thickness = placental thickness (z-dimension)
    # ellipticity = ratio of y to x axis dimensions
    # X,Y,Z spacing is the number of elements required in each of the x, y z directions

    # Calculate the dimensions of the ellipsoid
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    
    # z height of ellipsoid is 2* zradius
    # We want number of nodes to cover height and have prescribed spaing
    nnod_x = int(np.ceil(x_radius * 2.0 / x_spacing)) + 1
    x_width = x_spacing * (nnod_x - 1)
    nnod_y = int(np.ceil(y_radius * 2.0 / y_spacing)) + 1
    y_width = y_spacing * (nnod_y - 1)
    nnod_z = int(np.ceil(z_radius * 2.0 / z_spacing)) + 1
    z_width = z_spacing * (nnod_z - 1)
    
    node_loc=gen_rectangular_node(x_width, y_width, z_width, nnod_x, nnod_y, nnod_z)
    # Generating the element connectivity of each cube element, 8 nodes for each 3D cube element
    elems=cube_mesh_connectivity(nnod_x,nnod_y,nnod_z)

    return {'nodes': node_loc, 'elems': elems, 'total_nodes': nnod_x * nnod_y * nnod_z,
            'total_elems': (nnod_x - 1) * (nnod_y - 1) * (nnod_z - 1)}


def gen_ellip_mesh_tet(volume,thickness,ellipticity,n):
    """ Generates ellipsoid tetrahedral mesh for 3D problems

    Inputs:
       - volume: volume of placental ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions
       - n: number of datapoints generated to create the mesh

    Returns:
       - nodes: nodes location of mesh
       - elems: element connectivity of mesh (tetrahedral element)
       - node_array: array of nodes
       - element_array: array of elements
    """

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']
    
    nodeSpacing = (n/(2*x_radius*2*y_radius*2*z_radius)) **(1./3)
        
    nnod_x=2*x_radius*nodeSpacing
    nnod_y=2*y_radius*nodeSpacing
    nnod_z=2*z_radius*nodeSpacing
    nodes=gen_rectangular_node(x_radius*2, y_radius*2, z_radius*2, nnod_x, nnod_y, nnod_z)
    
    #nodes inside the ellipsoid    
    ellipsoid_node=np.zeros((len(nodes),3))
    count=0
    for nnode in range (0, len(nodes)):
       coord_point = nodes[nnode][0:3]
       inside=pg_utilities.check_in_on_ellipsoid(coord_point[0], coord_point[1], coord_point[2], x_radius, y_radius, z_radius)
       if inside:
          ellipsoid_node[count,:]=coord_point[:]
          count=count+1
    ellipsoid_node.resize(count,3)     
    xyList = ellipsoid_node[:,[0,1]]
    xyListUnique = np.vstack({tuple(row) for row in xyList})
    #looking for z_coordinate of surface nodes
    for xyColumn in xyListUnique:
        
        xyNodes = np.where(np.all(xyList == xyColumn, axis = 1))[0]
        if len(xyNodes) > 1:
           x_coord=ellipsoid_node[xyNodes[0],0]
           y_coord=ellipsoid_node[xyNodes[0],1]
           ellipsoid_node[xyNodes[len(xyNodes) - 1],2] = pg_utilities.z_from_xy(x_coord, y_coord, x_radius, y_radius, z_radius)   
           ellipsoid_node[xyNodes[0],2] =-1*( pg_utilities.z_from_xy(x_coord, y_coord, x_radius, y_radius, z_radius))

    #generate tetrahedral mesh
    pyMesh = Delaunay(ellipsoid_node)

    #Build arrays to pass into openCMISS conversion:
    node_loc = pyMesh.points
    temp_elems = pyMesh.simplices
    #CHECK ELEMENTS FOR 0 VOLUME:
    min_vol = 0.00001
    index = 0
    indexArr = []

    for element in temp_elems:
     x_coor = []
     y_coor = []
     z_coor = []
     for node in element:
        x_coor.append(node_loc[node][0])
        y_coor.append(node_loc[node][1])
        z_coor.append(node_loc[node][2])
          
     vmat = np.vstack((x_coor,y_coor,z_coor,[1.0,1.0,1.0,1.0]))#matrix of coor of element
     elem_volume = (1/6.0) * abs(np.linalg.det(vmat))#volume of each tetrahedral element
     
     #if volume is not zero
     if elem_volume > min_vol:
       
       indexArr.append(index)
     index = index+1
     
    #update arrays without 0 volume elements, to pass into openCMISS
    elems = temp_elems[indexArr,:]
    for i in range(len(elems)):
       elems[i] = [x+1 for x in elems[i]]
    element_array = range(1, len(elems)+1)
    node_array = range(1, len(node_loc)+1)

    return {'nodes': node_loc, 'elems': elems, 'element_array':element_array,'node_array': node_array,'nodeSpacing':nodeSpacing}


def gen_rectangular_node(x_width, y_width, z_width, nnod_x, nnod_y, nnod_z):
      
    # Create linspaces for x y and z coordinates
    x = np.linspace(-x_width / 2.0, x_width / 2.0, nnod_x)  # linspace for x axis
    y = np.linspace(-y_width / 2.0, y_width / 2.0, nnod_y)  # linspace for y axis
    z = np.linspace(-z_width / 2.0, z_width / 2.0, nnod_z)  # linspace for z axis
    node_loc_temp = np.vstack(np.meshgrid(y, z, x)).reshape(3, -1).T  # generate nodes for rectangular mesh
    node_loc = np.zeros((len(node_loc_temp),3))
    for i in range(0,len(node_loc)):
        node_loc[i][0] = node_loc_temp[i][2]
        node_loc[i][1] = node_loc_temp[i][0]
        node_loc[i][2] = node_loc_temp[i][1]
    
    return node_loc

def gen_rectangular_mesh2(nel_x,nel_y,nel_z,xdim,ydim,zdim,element_type):
    #generates a rectangular mesh of defined dimenions using either linear or quadratic elements
    if element_type==1: #linear element
       nnod_x= int(nel_x+1)
       nnod_y = int(nel_y+1)
       nnod_z = int(nel_z+1)
    elif element_type==2: #quadratic element
       nnod_x =  int((nel_x*2)+1)
       nnod_y =  int((nel_y*2)+1)
       nnod_z =  int((nel_z*2)+1)

    node = gen_rectangular_node(xdim, ydim, zdim, nnod_x, nnod_y, nnod_z)  # getting nodes

    if element_type == 1:  # linear element
        elems = cube_mesh_connectivity(nnod_x, nnod_y, nnod_z)  # getting elem connectivity
    elif element_type == 2:  # quadratic element
        elems = cube_mesh_connectivity_quadratic(nel_x, nel_y, nel_z, nnod_x, nnod_y,
                                                 nnod_z)  # getting element connectivity

    element_array = range(1, len(elems)+1)
    node_array = range(1, len(node)+1)
    if element_type==2:
        surfacenodes = identify_surface_node_quad(nel_x,nel_y,nel_z)
    else:
        print("This element type has no implemented surface node definition")
        surfacenodes = 0

    return{'nodes':node,'elems':elems,'element_array':element_array,
           'node_array':node_array,'surface_nodes':surfacenodes}


def gen_placental_mesh(nel_x,nel_y,nel_z,volume,thickness,ellipticity,element_type):

    """ Generates ellipsoid placental mesh to solve 3D problems

    Inputs:
       - nel: number of element in x,y,z axis , the more nel, the rounder the mesh
       - volume: volume of placental ellipsoid
       - thickness: placental thickness (z-dimension)
       - ellipticity: ratio of y to x axis dimensions
       
    Returns:
       - placental_node_coor: nodes location of mesh
       - placental_el_con: element connectivity of mesh (tetrahedral element)
       - node_array: array of nodes
       - element_array: array of elements
    """
    #creating cube between -1 and 1 with n number of element 
    #cubelength=2
    if element_type==1: #linear element
       nnod_x= int(nel_x+1)
       nnod_y = int(nel_y+1)
       nnod_z = int(nel_z+1)
    elif element_type==2: #quadratic element
       nnod_x =  int((nel_x*2)+1)
       nnod_y =  int((nel_y*2)+1)
       nnod_z =  int((nel_z*2)+1)

    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    cube_node = gen_rectangular_node(2*x_radius, 2*y_radius, 2*z_radius, nnod_x, nnod_y, nnod_z)
    if element_type == 1:  # linear element
        cube_elems = cube_mesh_connectivity(nnod_x, nnod_y, nnod_z)  # getting elem connectivity
    elif element_type == 2:  # quadratic element
        cube_elems = cube_mesh_connectivity_quadratic(nel_x, nel_y, nel_z, nnod_x, nnod_y,
                                                 nnod_z)  # getting element connectivity

    ellipsoid_coor = np.zeros((len(cube_node),3))

    for ii in range(0, len(cube_node)):
        ellipsoid_coor[ii, 0] = cube_node[ii, 0] * np.sqrt(1.0 - cube_node[ii, 1] ** 2 / (2.0*y_radius **2) -
                                cube_node[ii, 2] ** 2 / (2.0 * z_radius**2) + cube_node[ii, 1] ** 2 *
                                cube_node[ii,2] ** 2 / (3.0* y_radius**2 * z_radius**2))  # for  x_coor
        ellipsoid_coor[ii, 1] = cube_node[ii, 1] * np.sqrt(1.0 - cube_node[ii, 0] ** 2 / (2.0*x_radius **2) -
                                cube_node[ii, 2] ** 2 / (2.0 * z_radius**2) + cube_node[ii, 0] ** 2 * cube_node[ii, 2] ** 2
                                / (3.0* x_radius**2 * z_radius**2))  # for  y_coor
        ellipsoid_coor[ii, 2] = cube_node[ii, 2] * np.sqrt(1.0 - cube_node[ii, 1] ** 2 / (2.0*y_radius **2) -
                                cube_node[ii, 0] ** 2 / (2.0 * x_radius**2) + cube_node[ii, 1] ** 2 * cube_node[ii, 0] ** 2
                                / (3.0* y_radius**2 * x_radius**2))  # for  z_coor


    element_array = range(1, len(cube_elems)+1)
    node_array = range(1, len(ellipsoid_coor)+1)
    if element_type==2:
        surfacenodes = identify_surface_node_quad(nel_x,nel_y,nel_z)
    else:
        print("This element type has no implemented surface node definition")
        surfacenodes = 0
   
    return{'placental_node_coor':ellipsoid_coor,'placental_el_con':cube_elems,'element_array':element_array,
           'node_array':node_array,'surface_nodes':surfacenodes}


def cube_mesh_connectivity(nnod_x,nnod_y,nnod_z):
    """Generates element connectivity in cube mesh
      
       Inputs:
         - nnod_x:number of node in x axis
         - nnod_y:number of node in y axis
         - nnod_z:number of node in z axis

       Outputs:
         - elems: array of element connectivity

    """
    num_elems = (nnod_x - 1) * (nnod_y - 1) * (nnod_z - 1)
    elems = np.zeros((num_elems, 9),dtype=int)  # this stores first element number and then the nodes of each mesh element
    element_number = 0
    ne = 0
    # loop through elements
    for k in range(1, nnod_z):
        for j in range(1, nnod_y):
            for i in range(1, nnod_x):
                elems[ne][0] = ne  # store element number
                elems[ne][1] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1) #lowest coordinates
                elems[ne][2] = elems[ne][1] + 1 #add one in x
                elems[ne][3] = elems[ne][1] + nnod_x #go through x and find first in y
                elems[ne][4] = elems[ne][3] + 1 #add one in y
                elems[ne][5] = elems[ne][1] + nnod_x * nnod_y #same as 1 -4 but at higher z -coord
                elems[ne][6] = elems[ne][2] + nnod_x * nnod_y
                elems[ne][7] = elems[ne][3] + nnod_x * nnod_y
                elems[ne][8] = elems[ne][4] + nnod_x * nnod_y
                ne = ne + 1

    return elems

def cube_mesh_connectivity_quadratic(nel_x,nel_y,nel_z,nnod_x,nnod_y,nnod_z):
    """Generates element connectivity in quadratic cube mesh
      
       Inputs:
         - nnod_x:number of node in x axis
         - nnod_y:number of node in y axis
         - nnod_z:number of node in z axis

       Outputs:
         - elems: array of element connectivity in quadratic

    """

    num_elems = nel_x*nel_y*nel_z
    elems = np.zeros((num_elems, 28),dtype=int)  
    
    element_number = 0
    ne = 0
    #Got the element  
    for k in range(1,nnod_z,2):
        for j in range(1, nnod_y,2):
            for i in range(1,nnod_x, 2):
                #1st layer
                elems[ne][0] = ne  
                elems[ne][1] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1) #1st node
                elems[ne][2] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1) +1#right subsequent node
                elems[ne][3] = (i - 1) + (nnod_x) * (j - 1) + nnod_x * nnod_y * (k - 1) +2#right subsequent node
                elems[ne][4] = elems[ne][1] + nnod_x #1st node in another y layer
                elems[ne][5] = elems[ne][1] + nnod_x +1#right subsequent node
                elems[ne][6] = elems[ne][1] + nnod_x +2#right subsequent node
                elems[ne][7] = elems[ne][1] + 2*(nnod_x )#1st node in another y layer
                elems[ne][8] = elems[ne][1] + 2*(nnod_x )+1#right subsequent node
                elems[ne][9] = elems[ne][1] + 2*(nnod_x )+2#right subsequent node

                #2nd layer
                elems[ne][10] = elems[ne][1] + nnod_x * nnod_y #same in one z layer
                elems[ne][11] = elems[ne][2] + nnod_x * nnod_y
                elems[ne][12] = elems[ne][3] + nnod_x * nnod_y
                elems[ne][13] = elems[ne][4] + nnod_x * nnod_y
                elems[ne][14] = elems[ne][5] + nnod_x * nnod_y
                elems[ne][15] = elems[ne][6] + nnod_x * nnod_y             
                elems[ne][16] = elems[ne][7] + nnod_x * nnod_y
                elems[ne][17] = elems[ne][8] + nnod_x * nnod_y
                elems[ne][18] = elems[ne][9] + nnod_x * nnod_y

                #thrid layer

                elems[ne][19] = elems[ne][1] + nnod_x * nnod_y*2 #same in another z layer
                elems[ne][20] = elems[ne][2] + nnod_x * nnod_y*2
                elems[ne][21] = elems[ne][3] + nnod_x * nnod_y*2
                elems[ne][22] = elems[ne][4] + nnod_x * nnod_y*2
                elems[ne][23] = elems[ne][5] + nnod_x * nnod_y*2
                elems[ne][24] = elems[ne][6] + nnod_x * nnod_y*2                 
                elems[ne][25] = elems[ne][7] + nnod_x * nnod_y*2
                elems[ne][26] = elems[ne][8] + nnod_x * nnod_y*2
                elems[ne][27] = elems[ne][9] + nnod_x * nnod_y*2


                ne = ne + 1
               
    return elems



def identify_surface_node_quad(nel_x,nel_y,nel_z):
    """Generates collection of nodes that are on the surface of in quadratic placental mesh
      
       Inputs:
         - nel_x:number of elem in x axis
         - nel_y:number of elem in y axis
         - nel_z:number of elem in z axis

       Outputs:
         - surfacenode: collection of nodes on the surface of placental mesh

    """
    nnod_x =  int((nel_x*2)+1)#number of nodes in x axis
    nnod_y =  int((nel_y*2)+1)#number of nodes in y axis
    nnod_z =  int((nel_z*2)+1)#number of nodes in z axis
    #For left and right surface
    sIEN=np.zeros((9,nel_y*nel_z),dtype=int)#to store surface indiviaul element nodes (sIEN)
    e=0
    for k in range( 1,nnod_x*nnod_y*(nnod_z-1),(nnod_x*nnod_y)*2):#go up
            for j in range(  1,nnod_x*(nnod_y-1),2*nnod_x):#go left         
          
                sIEN[0,e] = j+(k-1) #1st node
                sIEN[1,e] = sIEN[0,e]+(nnod_x)*(nnod_y)#2nd node
                sIEN[2,e] = sIEN[1,e]+(nnod_x)*(nnod_y)#3rd node
                sIEN[3,e] = sIEN[0,e]+nnod_x#4th node
                sIEN[4,e] = sIEN[1,e]+nnod_x#5th node
                sIEN[5,e] = sIEN[2,e]+nnod_x#6th node
                sIEN[6,e] = sIEN[3,e]+nnod_x#7th node
                sIEN[7,e]=sIEN[4,e]+nnod_x#8th node
                sIEN[8,e]=sIEN[5,e]+nnod_x#9th node
                e = e+1            
   
            left=np.unique(sIEN)#collection of nodes of left surface
            right=np.unique(sIEN.T+(nnod_x-1))#collection of nodes on right surface

    #For front and back surface
    sIEN=np.zeros((9,nel_x*nel_z),dtype=int)
    e=0
    for k in range (1,nnod_x*nnod_y*(nnod_z-2),(nnod_x*nnod_y)*2):#go up
            for i in range( 1,nnod_x-1,2):#go right            
                sIEN[0,e] = i+(k-1)           
                sIEN[1,e] = sIEN[0,e]+1
                sIEN[2,e] = sIEN[0,e]+2                    
                sIEN[3,e]=sIEN[0,e]+(nnod_x*nnod_y)
                sIEN[4,e] = sIEN[3,e]+1
                sIEN[5,e] = sIEN[3,e]+2
                sIEN[6,e] = sIEN[3,e]+(nnod_x*nnod_y)
                sIEN[7,e] = sIEN[6,e]+1    
                sIEN[8,e] = sIEN[6,e]+2            
                e = e+1      
    
            front=np.unique(sIEN)#collection of nodes on front surface
            back=np.unique(sIEN.T+(nnod_x*(nnod_y-1)))#collection of nodes on back surface

    #For top and bottom surface
    sIEN=np.zeros((9,nel_x*nel_y),dtype=int)
    e=0
    for j in range( 1,nnod_x*(nnod_y-1),nnod_x*2):#go up
            for i in range (  1,nnod_x-1,2):#go back                                   
                sIEN[0,e] = i+(j-1)
                sIEN[1,e] = sIEN[0,e]+1
                sIEN[2,e] = sIEN[0,e]+2
                sIEN[3,e] = sIEN[0,e]+nnod_x          
                sIEN[4,e] = sIEN[3,e]+1
                sIEN[5,e] = sIEN[3,e]+2
                sIEN[6,e] = sIEN[3,e]+nnod_x          
                sIEN[7,e] = sIEN[6,e]+1
                sIEN[8,e] = sIEN[6,e]+2
                e = e+1

            bottom=np.unique(sIEN)    #collection of nodes on bottom surface
            top=np.unique(sIEN.T+(nnod_x*nnod_y)*(nnod_z-1))#collection of nodes on top surface

    surfacenode=np.hstack((front,back,left,right,bottom,top))
    surfacenode=np.unique(surfacenode)#collection of surface nodes from all surface
    return surfacenode

def identify_node_from_coord(nodes,filename):
    # reading in the node location
    xyz = open(filename, 'r')
    xyz_coor = xyz.readlines()  # readlines
    startLines = range(0, len(xyz_coor))

    for i in range(len(xyz_coor)):
        xyz_coor[i] = xyz_coor[i].split()
    xyzList = []
    for i in startLines:
        targetpoint = []
        targetpoint.append(float(xyz_coor[i][0]))  # x coor
        targetpoint.append((float(xyz_coor[i][1])))  # y coor
        targetpoint.append((float(xyz_coor[i][1])))  # y coor
        xyzList.append(targetpoint)
    xyz.close()

    node_list = np.zeros(len(xyzList))

    mindist = 100000
    for i in range(0,len(xyzList)):
        print xyzList[i]
        for j in range(0,len(nodes)):
            print(xyzList[i][0],nodes[j][0])

    return i


def identify_vessel_node(ellipsoid_coor,surfacenode,stem_file):
    """Generates array of spiral artery nodes and decidual vein nodes. Spiral artery nodes are mapped with stem villi.
      
       Inputs:
         - ellipsoid_coor:coordinate of nodes of placental mesh
         - surfacenode:array of surface nodes
         - stem_file:txt file that described stem villi locations

       Outputs:
         - spiral_array: array of spiral artery nodes 
         - decidual_array: array of decidual artery nodes
         - vesselnode: array of both spiral and decidual nodes
         - surfnode_ex_vessel: array of surface node excluding vessel nodes
    """    
    xyList=np.zeros((len(surfacenode),2))
    count=0
    for i in range (0,len(surfacenode)):#taking only x and y coordinates
     if ellipsoid_coor[surfacenode[i]-1,2]>0:#take if uppersurface node cos we are looking for the vessel node location from upper surface nodes only
      xyList[count,0]=ellipsoid_coor[surfacenode[i]-1,0]
      xyList[count,1]=ellipsoid_coor[surfacenode[i]-1,1]
      count=count+1

    xyList=xyList[0:count,:]

    vesselnode_temp = np.vstack({tuple(row) for row in xyList})#Remove duplicates#mignt not need this one

    #reading in the stem vessel to map the spiral artery location
    stem_xy = open(stem_file, 'r')
    stem_coor = stem_xy.readlines()#readlines
    startLines = range(0,len(stem_coor))

    for i in range(len(stem_coor)):
         stem_coor[i] = stem_coor[i].split()
    stem_xyList = []
    for i in startLines:
         node = []
         node.append(float(stem_coor[i][0]))#x coor of stem villi
         node.append((float(stem_coor[i][1]))) #y coor of stem villi     
         stem_xyList.append(node)
    stem_xy.close()

    vessel_mapped_stem=stem_xyList#this is the x,y location where we want to put spiral artery
    spiral_array=np.zeros((len(vessel_mapped_stem)),dtype=int)#store the node nuber of spiral artery
    decidual_array=np.zeros((len(vessel_mapped_stem)),dtype=int)#store the node number of decidual vein

    check=ellipsoid_coor[:,0:2]  
    np.random.seed(0) 
    for i in range(0,len(vessel_mapped_stem)): #for each blood vessel,Cycle through to find closest nodes
     
         distance = []
         for nodeX in vesselnode_temp:
           distance.append(np.sqrt((vessel_mapped_stem[i][0] - nodeX[0])**2 + (vessel_mapped_stem[i][1] - nodeX[1])**2))#distance from the nodes
         
         A=sorted(distance)[0]#taking the nearest node
         V=np.random.choice(sorted(distance)[1:])#choosing random , but it won't repeat with artery since it is started from 1, [0] is artery
          
         arterynode= vesselnode_temp[np.where(distance==A)]#this is artery node coor        
         veinnode=vesselnode_temp[np.where(distance==V)]#this is vein node coor

         xyNodes_A = np.where(np.all(check == arterynode[0], axis = 1))[0]#collection of nodes number that have same x,y of artery nodes 
         xyNodes_A = [x+1 for x in xyNodes_A]#adding 1 to get the right node number
     
         xyNodes_V = np.where(np.all(check == veinnode[0], axis = 1))[0]#collection of nodes number that have x,y of vein nodes 
         xyNodes_V = [x+1 for x in xyNodes_V]#adding 1 to get the right node number

         spiral_array[i] = xyNodes_A[len(xyNodes_A) - 1]#just taking the last node (cos we want from top surface)
         decidual_array[i] = xyNodes_V[len(xyNodes_V) - 1]#just taking the last node (cos we want from top surface)
         #remove the taken artery and vein nodes so that won't select again
         vesselnode_temp=np.delete(vesselnode_temp,np.where(np.all(arterynode[0] ==vesselnode_temp, axis=1)),0)#remove taken artery node             
         vesselnode_temp=np.delete(vesselnode_temp,np.where(np.all(veinnode [0]==vesselnode_temp, axis=1)),0)#remove taken vein node

    vesselnode=np.hstack((spiral_array,decidual_array))#array of vessel nodes (both artery and vein)
    index=[]
    for i in range (0,len(vesselnode)):
      idx= np.where(surfacenode==vesselnode[i])#index where surface node = vessel node to remove
      index.append(idx)

    surfnode_ex_vessel = np.delete(surfacenode, index)#array of surface node excluding vessel node

    return {'spiral_array':spiral_array,'decidual_array':decidual_array,'vesselnode':vesselnode,'surfnode_ex_vessel':surfnode_ex_vessel}
