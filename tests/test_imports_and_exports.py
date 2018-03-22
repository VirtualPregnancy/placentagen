import os
from unittest import TestCase
import unittest
import numpy as np

import placentagen

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exnode')
TESTDATA_FILENAME1 = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exelem')


class Test_import_exnode_exelem(TestCase):

    def test_num_nodes(self):
        nodedata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        self.assertTrue(nodedata['total_nodes'] is 4)

    def test_node_array_setup(self):
        nodedata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        node_array = nodedata['node_array']
        self.assertTrue(np.isclose(node_array[2][2], -0.5000000000000000E+01))



    def test_el_num(self):
       el_num= placentagen.import_exelem_tree(TESTDATA_FILENAME1)
       self.assertTrue(el_num['total_el'] is 3)


    def test_el_array_setup(self):
        eldata = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        el_array = eldata['el_array']
        self.assertTrue(np.isclose(el_array[0][2], 2))



if __name__ == '__main__':
   
    unittest.main()
