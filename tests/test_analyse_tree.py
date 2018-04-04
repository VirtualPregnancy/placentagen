
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
        
        eldata=placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        nodedata=placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br=placentagen.terminal_branch(eldata['el_array'])
        pos_node= placentagen.pos_axis(nodedata['node_array'],2.5,5,3.54)
        print term_br['terminal_el']
        print eldata['el_array']
        print pos_node['new_node_array']
        term_block=placentagen.terminal_block(term_br['terminal_el'],2,5,4,eldata['el_array'],pos_node['new_node_array'],0,0,0,2.5,2,2)
        self.assertTrue(np.isclose(term_block['total_term_block'], 2))

if __name__ == '__main__':
   unittest.main()
