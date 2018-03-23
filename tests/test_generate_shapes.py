from unittest import TestCase
import unittest
import numpy as np

import placentagen


class Test_generate_data(TestCase):
    def test_data_in_ellipsoid(self):
        thickness = (3.0 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        datapoints = placentagen.equispaced_data_in_ellipsoid(1, 1.0, thickness, 1.0)
        array_test = np.isclose(datapoints, [0.0, 0.0, 0.0])
        self.assertTrue(array_test.all)


class Test_create_trees(TestCase):
    def test_umilical_node(self):
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        datapoints = placentagen.umbilical_seed_geometry(1.0, thickness, 1.0, 0.0, 0.0, 0.1)
        self.assertTrue(datapoints['umb_nodes'][5][3],0.61833222)

    def test_umilical_node(self):
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        datapoints = placentagen.umbilical_seed_geometry(1.0, thickness, 1.0, 0.0, 0.0, 0.1)
        self.assertTrue(datapoints['umb_elems'][2][2], 4)