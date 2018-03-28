import numpy as np

aA=float(2.5)# x radius of ellipsoid
y_xratio=float(2)#y_x ratio of ellipsoid
bB=y_xratio*aA # y radius of ellipsoid
cC=float(3.54)#z radius of ellipsoid

Xlen=np.ceil(aA*2)# X lenght of meshgrid
Ylen=np.ceil(bB*2)# Y lenght of meshgrid
Zlen=np.ceil(cC*2)# Z lenght of meshgrid
x_min=0
y_min=0
z_min=0

x_max=Xlen
y_max=Ylen
z_max=Zlen
nel_x=int(Xlen/2) # number of element in x axis, 2 x 2 x 2 size is optimal
nel_y=int(Ylen/2) # number of element in y axis
nel_z=int(Zlen/2) # number of element in z axis

x_width=(x_max-x_min)/nel_x # the x width of each mesh grid cube
y_width=(y_max-y_min)/nel_y # the y width of each mesh grid cube
z_width=(z_max-z_min)/nel_z # the z width of each mesh grid cube


