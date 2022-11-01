import os
import placentagen as pg
#from grow_tree import *
from skeleton_to_tree import *
from analyse_tree import *
path = '/hpc/vsri355/Modelling/Modelling-files/test-inputs' #points to the folder containing the ex nodes and elem files
os.chdir(path)
file_name = "trialtree"  # file names of nodes and elems and takes in the node and elem files to be used
#for removing multiple connections in the read in tree
chorion_and_stem = {}  #initializing the tree geometry
chorion_and_stem['nodes'] = pg.import_exnode_tree(file_name+'.exnode')['nodes'][:,1:4]
chorion_and_stem['elems'] = pg.import_exelem_tree(file_name+'.exelem')['elems']
#populate element connectivity
elem_connectivity = pg.element_connectivity_1D(chorion_and_stem['nodes'],chorion_and_stem['elems'])
#removing multiple connections
chorion_and_stem,elem_connectivity=remove_multiple_elements(chorion_and_stem,elem_connectivity,type=None)
output = 1  # set whether output exfiles are made
if output == 1:
    # cmgui files to be exported
    path = '/hpc/vsri355/Modelling/Modelling-files/CMGUI_files-output/test'
    os.chdir(path) # filepath to save exported CMGUIfiles
pg.export_ex_coords(chorion_and_stem['nodes'], 'vessels', file_name, 'exnode') # exports nodes
pg.export_exelem_1d(chorion_and_stem['elems'], 'vessels', file_name) # exports elems


