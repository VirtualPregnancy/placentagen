
from unittest import TestCase

import numpy as np
import unittest
import placentagen
import os
TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exnode')
TESTDATA_FILENAME1 = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exelem')

class Test_Terminal_Br(TestCase):
        
    def test_terminal_br(self):
        eldata   = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br  = placentagen.calc_terminal_branch(noddata['nodes'],eldata['elems'])
        print(term_br['terminal_elems'])
        self.assertTrue(term_br['terminal_elems'][0] == 1)
        
    def test_terminal_br_total(self):
        eldata   = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br  = placentagen.calc_terminal_branch(noddata['nodes'],eldata['elems'])
        self.assertTrue(term_br['total_terminals'] == 2)



        
class test_term_br_location(TestCase):
        
    def test_term_br_loc(self):
        
        eldata=placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        nodedata=placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br  = placentagen.calc_terminal_branch(nodedata['nodes'],eldata['elems'])
        volume = 1  # mm^3
        thickness = 1  # mm
        ellipticity = 1.00  # no units
        spacing = 1.0  # mm
        rectangular_mesh = placentagen.gen_rectangular_mesh(volume, thickness, ellipticity, spacing, spacing, spacing)
        terminals_in_grid = placentagen.terminals_in_sampling_grid(rectangular_mesh, term_br, nodedata['nodes'])
        test_array=terminals_in_grid == [[0, 1, 0, 1]]
        self.assertTrue(test_array.all)
        

class Test_placental_vol(TestCase):
        
    def test_pl_vol(self):
        placenta_vol=1
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        rectangular_mesh=placentagen.gen_rectangular_mesh(placenta_vol,thickness, 1,1,1,1)
        p_vol=placentagen.placental_vol(rectangular_mesh, placenta_vol, thickness, 1,1,1,1)
               
        self.assertTrue(placenta_vol-0.02<= p_vol['total_pl_vol'] <= placenta_vol+0.02)



class Test_cal_vol_voxel(TestCase):
        
    def test_vol_br_in_voxel(self):
        placenta_vol=5
        thickness = 2
        rectangular_mesh=placentagen.gen_rectangular_mesh(placenta_vol,thickness, 1,1,1,1)
        eldata=placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        nodedata=placentagen.import_exnode_tree(TESTDATA_FILENAME)
        vol_voxel=placentagen.cal_vol_voxel(rectangular_mesh,eldata,nodedata,placenta_vol,thickness,1)

        self.assertTrue(vol_voxel['br_counting'][1,0] , 1)
        self.assertTrue(np.isclose(vol_voxel['total_vol_check']-vol_voxel['total_br_vol'],0))
        self.assertTrue(np.isclose(vol_voxel['total_br_vol'],0.0150341981624))
        self.assertTrue(np.isclose(vol_voxel['vol_each_br'][2,0],0.0055536036727))
        self.assertTrue(np.isclose(vol_voxel['master_vol_voxel'][1,0],2.13600141e-04))
        self.assertTrue(np.isclose(vol_voxel['master_vol_voxel'][1,1],2.13600141e-05))


       

if __name__ == '__main__':
   unittest.main()
