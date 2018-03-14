from unittest import TestCase
import os
import placentagen
import numpy as np

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'TestData/Small.exnode')

class Test_import_exnode_tree(TestCase):
    def test_num_nodes(self):
        nodedata=placentagen.import_exnode_tree(TESTDATA_FILENAME)
        self.assertTrue(nodedata['total_nodes'] is 4)
    def test_node_array_setup(self):
        nodedata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        node_array = nodedata['node_array']
        self.assertTrue(np.isclose(node_array[2][2],-0.5000000000000000E+01))

