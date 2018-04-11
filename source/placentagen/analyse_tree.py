#!/usr/bin/env python
import numpy as np
from . import pg_utilities
import sys
import math
import numpy.matlib 
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



def cal_vol_voxel(rectangular_mesh,eldata,nodedata,volume,thickness,ellipticity):
    '''
    This subroutine is to:
    1. calculate total volume of branches in each samp_grid_el (this is to calculate vol_fraction)
    2. calculate variables that will use to calculate weighted_daimeter of branches in each samp_grid_el in next subroutine
    3. count the number of branches in each samp_grid_el (to see the distribution)
    4. calculate the volume of individual branch in the whole tree
    5. calculate total volume of all branches in the whole tree (summation of 4 should be equal to 5)
     
    '''
    total_elems=rectangular_mesh['total_elems']#total number of samp_gr_element
    branch_node=nodedata['nodes']#node coordinate of branches whole tree
    
    branch_el=eldata['elems']#element connectivity of branches whole tree
    
    radii = pg_utilities.calculate_ellipse_radii(volume, thickness, ellipticity)#calculate radii of ellipsoid
    z_rad = radii['z_radius']
    x_rad = radii['x_radius']
    y_rad= radii['y_radius']
    
    pi=math.pi
    br_counting=np.zeros((total_elems,1),dtype=int)# this array stores the number of branch in each meshgrid element (to see distribution)
    vol_each_br=np.zeros((len(branch_el),1))# this array stores the vol of individual branch
    master_vol_voxel=np.zeros((total_elems, 2))#1st col of stores total vol of braches in each samp_grid el, 2nd col stores variable that will be used to calculate weighted_diameter later

    
    for ne in range (0,len(branch_el)):#looping for all branchs in tree
       N1=branch_node[branch_el[ne][1]][1:4]#coor of start node of a branch element
       N2=branch_node[branch_el[ne][2]][1:4]#coor of end node of a branch element
       
       r=0.05#######artifical value at the moment, radius of individual branch
       interval=0.01
       points=5
       unit_Vx=[1, 0, 0]# Defining Unit vector along the X-direction
       angle_X1X2 = np.arccos(np.dot(unit_Vx,np.subtract(N2,N1))/(np.linalg.norm(unit_Vx)*np.linalg.norm(np.subtract(N2,N1))) )*180/pi#branching angle from unit vector
       axis_rot = np.cross(unit_Vx,np.subtract(N2,N1) )# the axis of rotation (single rotation) to rotate the branch in % X-direction to the required arbitrary direction through cross product
       length_br=np.linalg.norm(np.subtract(N2,N1))#length of individual branch 
       branch_vol = pi*r**2*length_br#volume of individual branch 
       vol_each_br[ne,0]=branch_vol
       dp_along_length=np.linspace(interval,np.subtract(length_br,interval),10)#points along the length of branch

          
       #generate evenly spaced points in the y and z cross section of branch
       yp=np.linspace(-r,r,points)#linspace along cross-section of br
       zp=np.linspace(-r,r,points)#linspace along cross-section of br
       [yd,zd]=np.meshgrid(yp,zp)#generate meshgrid 

       counter=0
       dp_cross_section=np.zeros((2,len(yd)*len(yd[0])))
       
       for a in range(0,len(yd[0])):
           for b in range(0,len(yd[0])):
              if r-(math.sqrt(yd[a,b]**2+zd[a,b]**2)) >= 0: # include points only if points are inside or on boundary of cross-section (circle) of branch 
           
                 dp_cross_section[0,counter]=yd[a,b]
                 dp_cross_section[1,counter]=zd[a,b]
                 counter=counter+1

       dp_cross_section=dp_cross_section[:,0:counter]#datapoints inside or on boundary of cross-section
      
       dp_whole_br=np.matlib.repmat(dp_along_length[0], 1, counter)#create x vector
       
       for c in range (1,len(dp_along_length)):
          dp_whole_br=np.hstack([dp_whole_br, np.matlib.repmat(dp_along_length[c],1,counter)]) #x vector for multiple slices of branch cross-section
       
       dp_cross_section=np.matlib.repmat(dp_cross_section[:],1,len(dp_along_length)) #repeating yz coor of datapoint n times(n=dp_along_length)
       dp_whole_br=np.vstack((dp_whole_br,dp_cross_section)) #combine x and yz coor of datapoints
       
       # Rotate branch cylinder to correct position as specified by its two end node locations
       if angle_X1X2!=0 or angle_X1X2!=180:#if branching angle is not 0 or 180
         axis_rot=axis_rot[:]/np.linalg.norm(axis_rot)
         R=np.eye(3)

         for ii in range(0,3):
            v=R[:,ii]
            R[:,ii] = v*math.cos(math.radians(angle_X1X2)) + np.cross(axis_rot,v)*math.sin(math.radians(angle_X1X2)) + sum((axis_rot*v))*(1-math.cos(math.radians(angle_X1X2)))*axis_rot#Rodrigues' formula
 
       dp_whole_br_old=np.dot(R,dp_whole_br)
       N1=np.array([N1]).T
       dp_whole_br_new=dp_whole_br_old+N1#datapoints inside branch are corrected accordingly, element by element binary operation 
       N2=np.array([N2])
     
       if angle_X1X2==0: #if angle is 0, datapoints are corrected accordigly
          print('ZERO')
          dp_whole_br_new = dp_whole_br
   
          dp_whole_br_new[0,:] = -1*dp_whole_br_new[0,:]+N2[0,0]
          dp_whole_br_new[1,:] = dp_whole_br_new[1,:]+N2[0,1]
          dp_whole_br_new[2,:] = dp_whole_br_new[2,:]+N2[0,2]
    
       if angle_X1X2==180: #if angle is 180, datapoints are corrected accordinlgy
         print('ONE EIGHTY')
         dp_whole_br_new = dp_whole_br
         dp_whole_br_new[0,:] = dp_whole_br[0,:]+N2[0,0]
         dp_whole_br_new[1,:]= dp_whole_br_new[1,:]+N2[0,1]
         dp_whole_br_new[2,:] = dp_whole_br_new[2,:]+N2[0,2]
    

       datapoints = dp_whole_br_new.T# collection of coor of datapoints in each branch element
               
       #Rarely a very few datapoints may fall outside the ellipsoid. This can happen when the br is at the very peripheral part of ellipsoid and diameter is thick, so include only the datapoints inside
       Dx=[]
       Dy=[]
       Dz=[]
       for ndp in range (0,len(datapoints)):#checking if all datapoints fall within placental ellipsoid

          if (datapoints[ndp,0])**2/x_rad**2+(datapoints[ndp,1])**2/y_rad**2+(datapoints[ndp,2])**2/z_rad**2<=1:#if datapoints are inside or on the bourndary of ellipsoid

              Dx.append(datapoints[ndp][0])
              Dy.append(datapoints[ndp][1])
              Dz.append(datapoints[ndp][2])
       
     
       DP=np.vstack((Dx,Dy,Dz)).T
       
       if len(DP)==0:#countercheck
          print ('the whole branch is outside the ellipsoid.Something wrong',ne)
          sys.exit("CHECK ERROR")
       vol_voxel = np.zeros((total_elems,1))#this stores total vol in each samp_grid_el
       temp_br_counting=np.zeros((total_elems,1),dtype=int)
       points_in_grid = np.zeros(total_elems,dtype = int)
       #looking for how the datapoints(i.e.DP) are distributed in samp_gr_element (i.e. how individual branch are distributed in how many samp_grid_el)
       for samp_ne in range(0,total_elems):
        
        first_node = rectangular_mesh['elems'][samp_ne][1]#first node of sampl_grid_el
        last_node = rectangular_mesh['elems'][samp_ne][8]#last node of samp_grid_el
        min_coords = rectangular_mesh['nodes'][first_node][0:3]#coor of first node
        max_coords = rectangular_mesh['nodes'][last_node][0:3]#coor of last node
         
        for nt in range(0,len(DP)):
            in_element = [False, False, False]
            coord_DP = DP[nt]
            
            if coord_DP[0] >= min_coords[0] and coord_DP[0] < max_coords[0]:
                in_element[0] = True
                if coord_DP[1] >= min_coords[1] and coord_DP[1] < max_coords[1]:
                    in_element[1] = True
                    if coord_DP[2] >= min_coords[2] and coord_DP[2] < max_coords[2]:
                        in_element[2] = True
            
            
            if(np.all(in_element)):

                points_in_grid[samp_ne] = points_in_grid[samp_ne]+1

       
       if np.sum(points_in_grid) != len(DP):# to countercheck
          sys.exit("CHECK ERROR:not all datapoints are allocated to voxels")

       points_in_grid=points_in_grid.T
       
       vol_voxel= np.true_divide(points_in_grid,len(DP))*branch_vol

       if np.sum(vol_voxel)-vol_each_br[ne,0]>10**-6:#countercheck
          print ('branch number',ne)
          sys.exit("Error: vol not equal")
       
       br_search=np.where(vol_voxel!=0)#the samp_gr_el where parts of br are allocated
       
       temp_br_counting[br_search]=1
       
       br_counting[:,0]=br_counting[:,0] + temp_br_counting[:,0]#adding up the number of branches in each loop
       
       master_vol_voxel[:,0] = master_vol_voxel[:,0] + vol_voxel# adding up volume  as one samp_grid_el may have more than one branch element
       master_vol_voxel[:,1] = master_vol_voxel[:,1] + vol_voxel*r*2;#adding up variables
       
    if np.sum(master_vol_voxel[:,0])-np.sum(vol_each_br)>10**-6:#countercheck
       sys.exit("Error: total vol not equal")
    return{'master_vol_voxel': master_vol_voxel, 'br_counting':br_counting,'vol_each_br':vol_each_br,'total_br_vol':np.sum(vol_each_br),'total_vol_check':np.sum(master_vol_voxel[:,0])}










