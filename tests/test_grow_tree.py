from unittest import TestCase

import numpy as np

import placentagen


class Test_create_trees(TestCase):
    def test_umilical_node(self):
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        datapoints = placentagen.umbilical_seed_geometry(1.0, thickness, 1.0, 0.0, 0.0, 0.1, 20.0)
        self.assertTrue(datapoints['umb_nodes'][5][3], 0.61833222)

    def test_umilical_node(self):
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        datapoints = placentagen.umbilical_seed_geometry(1.0, thickness, 1.0, 0.0, 0.0, 0.1, 20.0)
        self.assertTrue(datapoints['umb_elems'][2][2], 4)
