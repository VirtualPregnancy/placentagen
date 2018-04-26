#!/usr/bin/env python
import numpy as np
from . import pg_utilities
import time
import sys
import math
import numpy.matlib

def calc_terminal_branch(node_loc, elems):
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

def terminals_in_sampling_grid_fast(rectangular_mesh, terminal_list, node_loc):
    # This function counts the number of terminals in a sampling grid element, will only work with
    # rectangular mesh created as in generate_shapes.gen_rectangular_mesh
    # inputs are:
    # Rectangular mesh - the sampling grid
    # terminal_list - a list of terminal,
    # node_loc - location of nodes
    num_terminals = terminal_list['total_terminals']
    terminals_in_grid = np.zeros(len(rectangular_mesh['elems']), dtype=int)
    terminal_elems = np.zeros(num_terminals, dtype=int)
    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']
    startx = np.min(nodes[:,0])
    xside = nodes[elems[0][8]][0]-nodes[elems[0][1]][0]
    endx=np.max(nodes[:,0])
    nelem_x=(endx-startx)/xside
    starty = np.min(nodes[:,1])
    yside = nodes[elems[0][8]][1]-nodes[elems[0][1]][1]
    endy=np.max(nodes[:,1])
    nelem_y = (endy - starty) / yside
    startz = np.min(nodes[:,2])
    zside = nodes[elems[0][8]][2]-nodes[elems[0][1]][2]
    endz=np.max(nodes[:,2])
    nelem_z = (endz - startz) / zside

    for nt in range(0, num_terminals):
        coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
        xelem_num = np.floor((coord_terminal[0] - startx) / xside)
        yelem_num = np.floor((coord_terminal[1] - starty) / yside)
        zelem_num = np.floor((coord_terminal[2] - startz) / zside)
        nelem = int(xelem_num + (yelem_num)*nelem_x + (zelem_num)*(nelem_x*nelem_y))
        terminals_in_grid[nelem]=terminals_in_grid[nelem]+1
        terminal_elems[nt] = nelem #record what element the terminal is in
    return{'terminals_in_grid': terminals_in_grid, 'terminal_elems': terminal_elems}


def terminals_in_sampling_grid(rectangular_mesh, placenta_list, terminal_list, node_loc):
    # This function counts the number of terminals in a sampling grid element
    # inputs are:
    # Rectangular mesh - the sampling grid
    # terminal_list - a list of terminal,
    # node_loc - location of nodes
    num_sample_elems = len(placenta_list)
    num_terminals = terminal_list['total_terminals']
    terminals_in_grid = np.zeros(len(rectangular_mesh['elems']), dtype=int)
    terminal_mapped = np.zeros(num_terminals, dtype=int)
    terminal_elems = np.zeros(num_terminals, dtype=int)

    for ne_i in range(0, num_sample_elems):
        # First node has min x,y,z and last node has max x,y,z
        ne = placenta_list[ne_i]
        if placenta_list[ne_i] > 0: #There is some placenta in this element (assuming none in el 0)
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
    return{'terminals_in_grid': terminals_in_grid, 'terminal_elems': terminal_elems}


def ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, num_test_points):
    # This subroutine calculates the placental volume associated with each element in a samplling grid
    # inputs are:
    # rectangular_mesh = the sampling grid nodes and elements
    # volume = placental volume
    # thickness = placental thickness
    # ellipiticity = placental ellipticity
    # num_test_points = resolution of integration quadrature
    total_elems = rectangular_mesh['total_elems']
    elems = rectangular_mesh['elems']
    nodes = rectangular_mesh['nodes']
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_radius = radii['z_radius']
    x_radius = radii['x_radius']
    y_radius = radii['y_radius']

    #Initialise the array that defines the volume of placenta in each grid element
    pl_vol_in_grid = np.zeros(total_elems)
    non_empty_loc = np.zeros(total_elems, dtype = int)
    non_empty_count = 0

    for ne in range(0, len(elems)):  # looping through elements
        count_in_range = 0
        nod_in_range = np.zeros(8, dtype=int)
        #define range of x, y , and z in the element
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
            pl_vol_in_grid[ne] = (endx - startx) * (endy - starty) * (endz - startz)  # the placental vol in that samp_grid_el is same as vol of samp_grid_el
        elif count_in_range == 0:  # if all 8 nodes are outside the ellpsiod
            # since this samp_grid_el is completely outside, the placental vol is zero
            pl_vol_in_grid[ne] = 0
        else:  # if some nodes in and some nodes out, the samp_grid_el is at the edge of ellipsoid
        #Use trapezoidal quadrature to caculate the volume under the surface of the ellipsoid in each element
            non_empty_loc[non_empty_count] = ne
            non_empty_count = non_empty_count + 1
            # need to map to positive quadrant
            repeat = False
            if (startz < 0 and endz <=0):
                #need to project to positive z axis
                startz = abs(nodes[elems[ne][8]][2])
                endz = abs(nodes[elems[ne][1]][2])
            elif(startz < 0 and endz > 0):
                #Need to split into components above and below the axis and sum the two
                startz = 0
                endz = abs(nodes[elems[ne][1]][2])
                startz_2=0
                endz_2 =nodes[elems[ne][8]][2]
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
                pl_vol_in_grid[ne] =pl_vol_in_grid[ne] + (Value1 - startz_2 * (endx - startx) * (
                                endy - starty))

    print('Number of Non-empty cells: ' + str(non_empty_count))
    print('Total number of cells: ' + str(total_elems))
    non_empty_loc = np.resize(non_empty_loc, non_empty_count)

    return {'pl_vol_in_grid': pl_vol_in_grid, 'non_empty_rects': non_empty_loc}



def cal_br_vol_samp_grid(rectangular_mesh,eldata,nodedata,volume,thickness,ellipticity,p_vol):
    '''
    This subroutine is to:
    1. calculate total volume of branches in each samp_grid_el (to use when calculate vol_frac/porosity)
    2. calculate total daimeter variable of branches in each samp_grid_el(to use when calculate wt_diam)
    3. count number of branches in each samp_grid_el
    4. calculate the volume of individual branch and total vol of all branches in the whole tree
         
    '''
    total_elems=rectangular_mesh['total_elems']#total number of samp_gr_element
    branch_node=nodedata['nodes']#node coordinate of branches whole tree
    branch_el=eldata['elems']#element connectivity of branches whole tree
    
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)#calculate radii of ellipsoid
    z_rad = radii['z_radius']
    x_rad = radii['x_radius']
    y_rad= radii['y_radius']
    
    pi=math.pi
    br_num_in_samp_gr=np.zeros((total_elems,1),dtype=int)# number of branch in each samp_grid element (useful to see distribution)
    vol_each_br=np.zeros((len(branch_el),1))#  the vol of each whole branch
    total_vol_samp_gr=np.zeros((total_elems, 2))#1st col of stores total vol of braches in each samp_grid el, 2nd col stores diameter variable of branches in each samp_grid el
    b_count=0
    
    for ne in range (0,len(branch_el)):#looping for all branchs in tree
       N1=branch_node[branch_el[ne][1]][1:4]#coor of start node of a branch element
       N2=branch_node[branch_el[ne][2]][1:4]#coor of end node of a branch element

       #check the branch is located inside the ellipsoid
       check_in_range_N1 = pg_utilities.check_in_ellipsoid(N1[0], N1[1],N1[2], x_rad, y_rad, z_rad)
       check_on_range_N1 = pg_utilities.check_on_ellipsoid(N1[0], N1[1],N1[2], x_rad, y_rad, z_rad)
       check_in_range_N2 = pg_utilities.check_in_ellipsoid(N2[0], N2[1],N2[2], x_rad, y_rad, z_rad)
       check_on_range_N2 = pg_utilities.check_on_ellipsoid(N2[0], N2[1],N2[2], x_rad, y_rad, z_rad)
       
       if (check_in_range_N1 ==False and check_on_range_N1 ==False) or (check_in_range_N2 ==False and check_in_range_N2 ==False):
          sys.exit('branch number: ' + str(ne)+' is located outside the ellispoid (whole or partial). Check ellipsoid vol/coordinates of br')
       
       r=0.1#######artifical value at the moment, radius of individual branch
       interval=0.01
       points=8
       unit_Vx=[1, 0, 0]# Defining Unit vector along the X-direction
       angle_N1N2 = np.arccos(np.dot(unit_Vx,np.subtract(N2,N1))/(np.linalg.norm(unit_Vx)*np.linalg.norm(np.subtract(N2,N1))) )*180/pi#branching angle from unit vector
       axis_rot = np.cross(unit_Vx,np.subtract(N2,N1) )# the axis of rotation (single rotation) to rotate the branch in % X-direction to the required arbitrary direction through cross product
       length_br=np.linalg.norm(np.subtract(N2,N1))#length of individual branch 
       branch_vol = pi*r**2*length_br#volume of individual branch 
       dp_along_length=np.linspace(interval,length_br,10)#linspace points along the length of branch
          
       #generate evenly spaced points in the y and z vector cross section of branch
       yp=np.linspace(-r,r,points)#linspace along of cross-section of br
       zp=np.linspace(-r,r,points)#linspace along of cross-section of br
       [yd,zd]=np.meshgrid(yp,zp)#generate meshgrid 
       
       counter=0
       dp_cross_section=np.zeros((2,len(yd)*len(yd[0])))
       
       for a in range(0,len(yd[0])):
           for b in range(0,len(yd[0])):
              if r-(math.sqrt(yd[a,b]**2+zd[a,b]**2)) >= 0: # include points only if points are inside or on boundary of cylinderical cross-section of br
           
                 dp_cross_section[0,counter]=yd[a,b]
                 dp_cross_section[1,counter]=zd[a,b]
                 counter=counter+1

       dp_cross_section=dp_cross_section[:,0:counter]#yz vector for one layer
       dp_cross_section=np.matlib.repmat(dp_cross_section[:],1,len(dp_along_length)) #repeating yz vector for all layers along length
       dp_whole_br=np.matlib.repmat(dp_along_length[0], 1, counter)#x vector for one layer
       
       for c in range (1,len(dp_along_length)):#x vector for all layer
          dp_whole_br=np.hstack([dp_whole_br, np.matlib.repmat(dp_along_length[c],1,counter)]) 
              
       dp_whole_br=np.vstack((dp_whole_br,dp_cross_section)) #combine x and yz vector of datapoints
       
       # Rotate branch cylinder to correct position as specified by its two end node locations
       if angle_N1N2!=0 or angle_N1N2!=180:#if angle is not 0 or 180
         axis_rot=axis_rot[:]/np.linalg.norm(axis_rot)
         R=np.eye(3)

         for ii in range(0,3):
            v=R[:,ii]
            R[:,ii] = v*math.cos(math.radians(angle_N1N2)) + np.cross(axis_rot,v)*math.sin(math.radians(angle_N1N2)) + sum((axis_rot*v))*(1-math.cos(math.radians(angle_N1N2)))*axis_rot#Rodrigues' formula
 
       dp_whole_br_old=np.dot(R,dp_whole_br)
       N1=np.array([N1]).T
       dp_whole_br_new=dp_whole_br_old+N1#datapoints inside branch are also corrected accordingly, element by element binary operation 
       N2=np.array([N2])
     
       if angle_N1N2==0: #if angle is 0, datapoints are corrected accordigly
          
          dp_whole_br_new = dp_whole_br
   
          dp_whole_br_new[0,:] = -1*dp_whole_br_new[0,:]+N2[0,0]
          dp_whole_br_new[1,:] = dp_whole_br_new[1,:]+N2[0,1]
          dp_whole_br_new[2,:] = dp_whole_br_new[2,:]+N2[0,2]
    
       if angle_N1N2==180: #if  angle is 180, datapoints are corrected accordinlgy
         
         dp_whole_br_new = dp_whole_br
         dp_whole_br_new[0,:] = dp_whole_br[0,:]+N2[0,0]
         dp_whole_br_new[1,:]= dp_whole_br_new[1,:]+N2[0,1]
         dp_whole_br_new[2,:] = dp_whole_br_new[2,:]+N2[0,2]
    

       datapoints = dp_whole_br_new.T# corrected coor of datapoints in one branch element
       
       vol_each_br[ne,0]=branch_vol
       b_count=b_count+1
       vol_samp_gr = np.zeros((total_elems,1))
       temp_br_num_in_samp_gr=np.zeros((total_elems,1),dtype=int)
       points_in_grid = np.zeros(total_elems,dtype = int)
       
       elems = rectangular_mesh['elems']
       nodes = rectangular_mesh['nodes']
       startx = np.min(nodes[:,0])
       xside = nodes[elems[0][8]][0]-nodes[elems[0][1]][0]
       endx=np.max(nodes[:,0])
       nelem_x=(endx-startx)/xside
       starty = np.min(nodes[:,1])
       yside = nodes[elems[0][8]][1]-nodes[elems[0][1]][1]
       endy=np.max(nodes[:,1])
       nelem_y = (endy - starty) / yside
       startz = np.min(nodes[:,2])
       zside = nodes[elems[0][8]][2]-nodes[elems[0][1]][2]
       endz=np.max(nodes[:,2])
       nelem_z = (endz - startz) / zside
       
       for num_points in range(0, len(datapoints)):
           coord_datapoints = datapoints[num_points]
           xelem_num = np.floor((coord_datapoints[0] - startx) / xside)
           yelem_num = np.floor((coord_datapoints[1] - starty) / yside)
           zelem_num = np.floor((coord_datapoints[2] - startz) / zside)
           nelem = int(xelem_num + (yelem_num)*nelem_x + (zelem_num)*(nelem_x*nelem_y))
           points_in_grid[nelem]=points_in_grid[nelem]+1
              
      
       points_in_grid=points_in_grid.T
    
       vol_samp_gr= np.true_divide(points_in_grid,len(datapoints))*branch_vol#distribute vol in different sampling grid
     
       br_search=np.where(vol_samp_gr!=0)#the samp_gr_el where br are allocated
       temp_br_num_in_samp_gr[br_search]=1
       br_num_in_samp_gr[:,0]=br_num_in_samp_gr[:,0] + temp_br_num_in_samp_gr[:,0]#adding up the number of branches in each loop
       
       total_vol_samp_gr[:,0] = total_vol_samp_gr[:,0] + vol_samp_gr# adding up volume  as one samp_grid_el may have more than one branch element
       total_vol_samp_gr[:,1] = total_vol_samp_gr[:,1] + vol_samp_gr*r*2;#adding up diam related variable
    

    pl_vol_in_grid=p_vol['pl_vol_in_grid']
    for i in range (0,len(total_vol_samp_gr)):#just countercheck
         if pl_vol_in_grid[i]==0 and total_vol_samp_gr[i,0]!=0 :#this should not happen
            sys.exit("some datapoints of branches are allocated outside ellipsoid")

    print('Total number of branch assessed, branch in ellipsoid =  ' + str(b_count))
    return{'total_vol_samp_gr': total_vol_samp_gr, 'br_num_in_samp_gr':br_num_in_samp_gr,'vol_each_br':vol_each_br,'total_br_vol':np.sum(vol_each_br)}
