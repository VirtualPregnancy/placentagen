
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

      
class test_pl_vol_in_grid(TestCase):
        
    def test_pl_vol_margin(self):
        thickness =  (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0  # mm
        ellipticity = 1.00  # no units
        spacing = 1.0  # mm
        volume=1
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = [[0., 0., 0.], [ thickness/2.0, 0., 0.],[0., thickness/2.0, 0.],[ thickness/2.0, thickness/2.0, 0.],[0., 0., thickness/2.0], [ thickness/2.0, 0., thickness/2.0],[0., thickness/2.0,thickness/2.0],[ thickness/2.0, thickness/2.0, thickness/2.0]]
        rectangular_mesh['elems'] = [[ 0,  0,  1,  2,  3,  4, 5, 6, 7]]
        rectangular_mesh['total_nodes'] =8
        rectangular_mesh['total_elems'] = 1
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 25)
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][0], 0.12485807941))
        self.assertTrue(abs(pl_vol['pl_vol_in_grid']-1./8.)/(1./8)<1e-2)#looking for less than 1% error in expected volume of 1/8

    def test_pl_vol_complete_inside(self):
        thickness =  2  # mm
        ellipticity = 1.6  # no units
        spacing = 0.5  # mm
        volume=5
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = [[-1., -1.5, -1.],[-0.5 ,-1.5, -1.],[-1., -1., -1.] ,[-0.5,-1., -1.],[-1.,-1.5,-0.5],[-0.5,-1.5,-0.5],[-1.,-1.,-0.5] ,[-0.5,-1.,-0.5]]
        rectangular_mesh['elems'] = [[0, 0, 1, 2, 3, 4, 5, 6, 7]]
        rectangular_mesh['total_nodes'] =8
        rectangular_mesh['total_elems'] = 1
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 25)
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][0], 0.0))
      

    def test_pl_vol_complete_outside(self):
        thickness =  2  # mm
        ellipticity = 1.6  # no units
        spacing = 0.5  # mm
        volume=5
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = [[-0.5,-0.5,-0.5],[ 0., -0.5,-0.5],[-0.5, 0.,-0.5],[ 0., 0. ,-0.5],[-0.5, -0.5 ,0. ],[ 0., -0.5 ,0.],[-0.5, 0.,0.],[0.,0.,0.]]
        rectangular_mesh['elems'] = [[0,  0,  1,  2,  3,  4, 5, 6, 7]]
        rectangular_mesh['total_nodes'] =8
        rectangular_mesh['total_elems'] = 1
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 0.125)
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][0], 0.125))
        self.assertTrue(abs(pl_vol['pl_vol_in_grid'][0]- (spacing*spacing*spacing))<1e-2)#not sure this one is need as it is supposed to be exactly the same

    
if __name__ == '__main__':
   unittest.main()
