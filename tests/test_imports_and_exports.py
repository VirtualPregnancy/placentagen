from unittest import TestCase
import os
import placentagen
import numpy as np

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'TestData/Small.exnode')


class Test_import_exnode_tree(TestCase):

    def test_num_nodes(self):
        def list_files(startpath):
            for root, dirs, files in os.walk(startpath):
                level = root.replace(startpath, '').count(os.sep)
                indent = ' ' * 4 * (level)
                print('{}{}/'.format(indent, os.path.basename(root)))
                subindent = ' ' * 4 * (level + 1)
                for f in files:
                    print('{}{}'.format(subindent, f))

        list_files(os.path.dirname(__file__))
        list_files(os.path.realpath(os.path.dirname(__file__)))
        nodedata=placentagen.import_exnode_tree(TESTDATA_FILENAME)
        self.assertTrue(nodedata['total_nodes'] is 4)

    def test_node_array_setup(self):
        nodedata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        node_array = nodedata['node_array']
        self.assertTrue(np.isclose(node_array[2][2],-0.5000000000000000E+01))

