
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
        term_br  = placentagen.terminal_branch(eldata['el_array'])
        self.assertTrue(np.isclose(term_br['terminal_el'][0,1], 2))
        
    def test_terminal_br_total(self):
        eldata   = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        term_br_total  = placentagen.terminal_branch(eldata['el_array'])
        self.assertTrue(np.isclose(term_br_total['total_terminal_el'], 2))


class Test_pos_axis(TestCase):
        
    def test_axis(self):
        NODE = placentagen. import_exnode_tree(TESTDATA_FILENAME)
        pos_coor  = placentagen.pos_axis(NODE['node_array'],2.5,5,3.54)
        
        self.assertTrue(np.isclose(pos_coor['new_node_array'][3][3], 0.2999999999999998))
        
class Test_term_br_location(TestCase):
        
    def test_term_br_loc(self):
        
        eldata=import_exelem_tree(TESTDATA_FILENAME1)
        nodedata=import_exnode_tree(TESTDATA_FILENAME)
        term_br=terminal_branch(eldata['el_array'])
        pos_node= pos_axis(nodedata['node_array'],2.5,5,3.54)
        term_block=terminal_block(term_br['terminal_el'],nel_x,nel_y,nel_z,eldata['el_array'],pos_node['new_node_array'],x_min,y_min,z_min,x_width,y_width,z_width)
        self.assertTrue(np.isclose(term_block['total_term_block'], 2))
        

if __name__ == '__main__':
   unittest.main()
