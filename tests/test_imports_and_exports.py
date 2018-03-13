from unittest import TestCase

import placentagen

class Test_import_exnode_tree(TestCase):
    def test_num_nodes(self):
        nodedata=placentagen.import_exnode_tree('TestData/Small.exnode')
        self.assertTrue(nodedata['total_nodes'] is 4)