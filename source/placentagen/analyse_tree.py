#!/usr/bin/env python
import numpy as np
from . import pg_utilities

def calc_terminal_branch(node_loc,elems):
    #This function generates a list of terminal nodes associated with a branching geometry
    #inputs are node locations and elements
    num_elems = len(elems)
    num_nodes = len(node_loc)
    elem_cnct = pg_utilities.element_connectivity_1D(node_loc, elems)

    terminal_branches = np.zeros(num_elems, dtype = int)
    terminal_nodes = np.zeros(num_nodes, dtype = int)

    num_term = 0
    for ne in range(0,num_elems):
        if elem_cnct['elem_down'][ne][0] == 0: #no downstream element
            terminal_branches[num_term] = ne
            terminal_nodes[num_term] = elems[ne][2] #node that leaves the terminal element
            num_term = num_term + 1

    terminal_branches = np.resize(terminal_branches,num_term)
    terminal_nodes = np.resize(terminal_nodes,num_term)

    print('Total number of terminals assessed, num_terminals =  ' + str(num_term))

    return {'terminal_elems': terminal_branches, 'terminal_nodes': terminal_nodes, 'total_terminals': num_term}


def terminals_in_sampling_grid(rectangular_mesh,terminal_list,node_loc):
    #This function counts the number of terminals in a sampling grid element
    #inputs are:
    #Rectangular mesh - the sampling grid
    #terminal_list - a list of terminals
    #node_loc - location of nodes
    num_sample_elems = rectangular_mesh['total_elems']
    num_terminals = terminal_list['total_terminals']
    terminals_in_grid = np.zeros(num_sample_elems,dtype = int)

    for ne in range(0,num_sample_elems):
        #First node has min x,y,z and last node has max x,y,z
        first_node = rectangular_mesh['elems'][ne][1]
        last_node = rectangular_mesh['elems'][ne][8]
        min_coords = rectangular_mesh['nodes'][first_node][0:3]
        max_coords = rectangular_mesh['nodes'][last_node][0:3]
        for nt in range(0,num_terminals):
            in_element = [False, False, False]
            coord_terminal = node_loc[terminal_list['terminal_nodes'][nt]][1:4]
            if coord_terminal[0] >= min_coords[0] and coord_terminal[0] < max_coords[0]:
                in_element[0] = True
                if coord_terminal[1] >= min_coords[1] and coord_terminal[1] < max_coords[1]:
                    in_element[1] = True
                    if coord_terminal[2] >= min_coords[2] and coord_terminal[2] < max_coords[2]:
                        in_element[2] = True
            if(np.all(in_element)):
                terminals_in_grid[ne] = terminals_in_grid[ne]+1

    return terminals_in_grid


def placental_vol(rectangular_mesh, volume, thickness, ellipticity,x_spacing,y_spacing,z_spacing):
    elems=rectangular_mesh['elems']
    nodes=rectangular_mesh['nodes']
    x=nodes[:,0]
    y=nodes[:,1]
    z=nodes[:,2]
    total_elems=rectangular_mesh['total_elems']

    
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)
    z_rad = radii['z_radius']
    x_rad = radii['x_radius']
    y_rad= radii['y_radius']

    pl_vol_in_grid=np.zeros((total_elems,3),dtype=object)#1st col stores the number of samp_grid_el, 2nd col stores 0/1/2 (depend on type of samp_grid_el), 3rd col stores the pl_vol in that sam_grid_el
    pl_vol_in_grid[:,0:1]=int(0)
    pl_vol_in_grid[:,2]=float(0)
    
    for ii in range(0,len(elems)):
    
       if ((x[elems[ii,1]])**2/x_rad**2+(y[elems[ii,1]])**2/y_rad**2+(z[elems[ii,1]])**2/z_rad**2)<=1 and ((x[elems[ii,2]])**2/x_rad**2+(y[elems[ii,2]])**2/y_rad**2+(z[elems[ii,2]])**2/z_rad**2)<=1 and ((x[elems[ii,3]])**2/x_rad**2+(y[elems[ii,3]])**2/y_rad**2+(z[elems[ii,3]])**2/z_rad**2)<=1 and ((x[elems[ii,4]])**2/x_rad**2+(y[elems[ii,4]])**2/y_rad**2+(z[elems[ii,4]])**2/z_rad**2)<=1 and ((x[elems[ii,5]])**2/x_rad**2+(y[elems[ii,5]])**2/y_rad**2+(z[elems[ii,5]])**2/z_rad**2)<=1 and ((x[elems[ii,6]])**2/x_rad**2+(y[elems[ii,6]])**2/y_rad**2+(z[elems[ii,6]])**2/z_rad**2)<=1 and ((x[elems[ii,7]])**2/x_rad**2+(y[elems[ii,7]])**2/y_rad**2+(z[elems[ii,7]])**2/z_rad**2)<=1 and ((x[elems[ii,8]])**2/x_rad**2+(y[elems[ii,8]])**2/y_rad**2+(z[elems[ii,8]])**2/z_rad**2)<=1 :#if all 8 nodes are inside the ellipsoid
         pl_vol_in_grid[ii,0]=ii
         pl_vol_in_grid[ii,1]=2 #if the sam_grid_element is completely inside the placental ellipsoid, it is defined as 2
         pl_vol_in_grid[ii,2]=x_spacing*y_spacing*z_spacing#the placental vol in that samp_grid_el is same as vol of samp_grid_el
         

       elif ((x[elems[ii,1]])**2/x_rad**2+(y[elems[ii,1]])**2/y_rad**2+(z[elems[ii,1]])**2/z_rad**2)>1 and ((x[elems[ii,2]])**2/x_rad**2+(y[elems[ii,2]])**2/y_rad**2+(z[elems[ii,2]])**2/z_rad**2)>1 and ((x[elems[ii,3]])**2/x_rad**2+(y[elems[ii,3]])**2/y_rad**2+(z[elems[ii,3]])**2/z_rad**2)>1 and ((x[elems[ii,4]])**2/x_rad**2+(y[elems[ii,4]])**2/y_rad**2+(z[elems[ii,4]])**2/z_rad**2)>1 and ((x[elems[ii,5]])**2/x_rad**2+(y[elems[ii,5]])**2/y_rad**2+(z[elems[ii,5]])**2/z_rad**2)>1 and ((x[elems[ii,6]])**2/x_rad**2+(y[elems[ii,6]])**2/y_rad**2+(z[elems[ii,6]])**2/z_rad**2)>1 and ((x[elems[ii,7]])**2/x_rad**2+(y[elems[ii,7]])**2/y_rad**2+(z[elems[ii,7]])**2/z_rad**2)>1 and ((x[elems[ii,8]])**2/x_rad**2+(y[elems[ii,8]])**2/y_rad**2+(z[elems[ii,8]])**2/z_rad**2)>1 :#if all 8 nodes are outside the ellpsiod

         pl_vol_in_grid[ii,0]=ii  
         pl_vol_in_grid[ii,1]=0#this type of samp_grid_el is defined as 0
         pl_vol_in_grid[ii,2]=0 #since this samp_grid_el is completely outside, the placental vol is zero (there will be no tree here and no need to worried about the vol)
         
       else: #if some nodes in and some nodes out, the samp_grid_el is at the edge of ellipsoid

                       
         pl_vol_in_grid[ii,1]=1 #this type of samp_grid_el is defined as 1. The placental vol cannot be fully cube, so has to calculate proportion of pl_vol
         #Filling the equally spaced datapoints in this samp_grid_el and looking for how many datapoints falls in the ellipsoid and can calculate volume 
         startx=x[elems[ii,1]]
         endx=x[elems[ii,2]]

         starty=y[elems[ii,1]]
         endy=y[elems[ii,3]]

         startz=z[elems[ii,1]]
         endz=z[elems[ii,5]]

         xVector = np.linspace(startx, endx, 30)
         yVector = np.linspace(starty, endy, 30)
         zVector = np.linspace(startz, endz, 30)
         
         nodes = np.vstack(np.meshgrid(xVector,yVector,zVector)).reshape(3,-1).T
         
         pointcount=0
         for jj in range (0, len(nodes)):
             if (nodes[jj,0])**2/x_rad**2+(nodes[jj,1])**2/y_rad**2+(nodes[jj,2])**2/z_rad**2<=1:#if the point fall inside the ellipsoid
                pointcount=pointcount+1
                                        
 
         pl_vol_in_grid[ii,2]=float(pointcount)/float(len(nodes))*(x_spacing*y_spacing*z_spacing)#calculate the proportion of placental vol in that samp_gr_el
         pl_vol_in_grid[ii,0]=ii
        

    return{'pl_vol_in_grid':pl_vol_in_grid,'total_pl_vol':np.sum(pl_vol_in_grid[:,2])}    
