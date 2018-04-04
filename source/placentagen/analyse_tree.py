#!/usr/bin/env python
import numpy as np

def terminal_branch(Br_El):#Br_El is el_array output from import_exeelem_tree subroutine

    T_Br=np.ones((len(Br_El),1))#if it is 1, it is terminal br, if it is zero it is not, initialise all branches as terminal

    for i in range(0, len(Br_El)):
        for j in range(0,len(Br_El)):     
        
            if Br_El[j][1] == Br_El[i][2]:# if the first node of any elment is equal to 2nd node another element, 
               T_Br[i] = 0 #the element is not terminal (zero means 'not terminal')
             
       
    Terminal_El_num=np.array(np.where(T_Br[:,0]==1))#looking for the element number which are terminal branch 
    
    return {'terminal_el': Terminal_El_num, 'total_terminal_el': len(Terminal_El_num[0])}
def pos_axis(NODE,x_rad,y_rad,z_rad):#NODE is the output array from nodedata['node_array'], x_rad,y_rad,z_rad are radius of ellipsoid#this function is to shift the whole tree into positive axis as the meshgrid in positive axis works easier and faster
    
    for i in range (0,len(NODE)):
      
         NODE[i][1] = NODE[i][1]+x_rad#to shift the x axis into positve axis
         NODE[i][2] = NODE[i][2]+y_rad#to shift the y axis into positive axis
         NODE[i][3] = NODE[i][3]+z_rad#to shift the z axis into positive axis
    return{'new_node_array':NODE}

def terminal_block(Terminal_El_num,nel_x,nel_y,nel_z,Br_El,node,x_min,y_min,z_min,x_width,y_width,z_width):
    term_block = np.zeros((nel_x*nel_y*nel_z,1))#this stores the 1 and 0 value. If it is 1, that meshgrid element contain terminal branch of tree
    for j in range(0,len(Terminal_El_num[0])):

       N=Br_El[Terminal_El_num[0][j]]# element connectivity of each terminal branch
  
       Endpoints=node[int(N[2])-1] # coordinates of 2nd node in element connectivity
  
       #searching where that node is located in the mesh grid cube element
       num_x = np.floor((Endpoints[1]-x_min)/x_width)+1;
       num_y = np.floor((Endpoints[2]-y_min)/y_width)+1;
       num_z = np.floor((Endpoints[3]-z_min)/z_width)+1;
   
       if num_x > nel_x:
          num_x = nel_x

       if num_y >= nel_y:
          num_y = nel_y
   
       if num_z >= nel_z:
          num_z = nel_z

       T_e = ((num_z-1)*nel_x*nel_y + (num_y-1)*nel_x + num_x)-1#the number of mesh grid element where that end node of term_br lie
       term_block[T_e,0] = 1;# if 1, that mesh grid cube contains terminal branch
 
    return {'term_block':term_block, 'total_term_block': int(np.sum(term_block))}


