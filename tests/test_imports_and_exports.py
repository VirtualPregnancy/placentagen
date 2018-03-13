from unittest import TestCase
import os
import placentagen

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'TestData/Small.exnode')

class Test_import_exnode_tree(TestCase):
    def test_num_nodes(self):
        nodedata=placentagen.import_exnode_tree(TESTDATA_FILENAME)
        self.assertTrue(nodedata['total_nodes'] is 4)