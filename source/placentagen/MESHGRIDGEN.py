#!/usr/bin/env python

import numpy as np

def meshgridgen(x_min,x_max,y_min,y_max,z_min,z_max,nel_x,nel_y,nel_z,x_width,y_width,z_width):
    x = np.linspace(x_min,x_max, x_max/x_width+1)#linspace for x axis
    y = np.linspace(y_min, y_max, y_max/y_width+1)#linspace for y axis
    z = np.linspace(z_min, z_max, z_max/z_width+1)#linspace for z axis
    nodes = np.vstack(np.meshgrid(y,z,x)).reshape(3,-1)#generate mesh

    y=np.array(nodes[0,:])#y coordinates
    z=np.array(nodes[1,:])#z coordinates
    x=np.array(nodes[2,:])#x coordinates

    #Generating the node number of each cube element, 8 nodes for each 3D cube element
    nodeOfelement=np.zeros((8,nel_x*nel_y*nel_z))
    E1=0

    for K1 in range (1,nel_z+1):
       for J1 in range (1,nel_y+1):
          for I1 in range(1,nel_x+1):
           
            nodeOfelement[0,E1] = I1+(nel_x+1)*(J1-1)+(nel_x+1)*(nel_y+1)*(K1-1);
            nodeOfelement[1,E1] = nodeOfelement[0,E1]+1;
            nodeOfelement[2,E1] = nodeOfelement[0,E1]+nel_x+1;
            nodeOfelement[3,E1] = nodeOfelement[2,E1]+1;
            nodeOfelement[4,E1] = nodeOfelement[0,E1]+(nel_x+1)*(nel_y+1);
            nodeOfelement[5,E1] = nodeOfelement[1,E1]+(nel_x+1)*(nel_y+1);
            nodeOfelement[6,E1] = nodeOfelement[2,E1]+(nel_x+1)*(nel_y+1);
            nodeOfelement[7,E1] = nodeOfelement[3,E1]+(nel_x+1)*(nel_y+1);
            
            E1 = E1+1;

    nodeOfelement=nodeOfelement.T
    return {'nodeOfelement': nodeOfelement, 'x_coor': x,'y_coor':y,'z_coor':z,'total_mesh_el':nel_x*nel_y*nel_z,'total_mesh_node':(nel_x+1)*(nel_y+1)*(nel_z+1)}


