
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
        
class test_pl_vol_in_grid(TestCase):
        
    def test_pl_vol_margin(self):
        
        volume = 1  # mm^3
        thickness =  (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0  # mm
        ellipticity = 1.00  # no units
        spacing = 1.0  # mm
        rectangular_mesh = placentagen.gen_rectangular_mesh(volume, thickness, ellipticity, spacing, spacing, spacing)
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 25)
        
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][0], 0.124841925586))


    def test_pl_vol_complete_outside(self):
        
        volume = 5  # mm^3
        thickness =  2 # mm
        ellipticity = 1.6  # no units
        spacing = 0.5  # mm
        rectangular_mesh = placentagen.gen_rectangular_mesh(volume, thickness, ellipticity, spacing, spacing, spacing)
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 25)
        
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][0], 0.0))

    def test_pl_vol_complete_inside(self):
        
        volume = 5  # mm^3
        thickness =  2 # mm
        ellipticity = 1.6  # no units
        spacing = 0.5  # mm
        rectangular_mesh = placentagen.gen_rectangular_mesh(volume, thickness, ellipticity, spacing, spacing, spacing)
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 25)
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][33], 0.125))

if __name__ == '__main__':
   unittest.main()
