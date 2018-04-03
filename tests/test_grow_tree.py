from unittest import TestCase

import numpy as np

import placentagen


class Test_create_trees(TestCase):
    def test_umilical_node(self):
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        data_input = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        datapoints = placentagen.umbilical_seed_geometry(1.0, thickness, 1.0, 0.0, 0.0, 0.1, 20.0, data_input)
        self.assertTrue(datapoints['umb_nodes'][5][3], 0.61833222)

    def test_umilical_node(self):
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        data_input = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        datapoints = placentagen.umbilical_seed_geometry(1.0, thickness, 1.0, 0.0, 0.0, 0.1, 20.0, data_input)
        self.assertTrue(datapoints['umb_elems'][2][2], 4)


class Test_grow_trees(TestCase):
    def test_grow_node(self):
        seed_geom = {}
        seed_geom['umb_nodes'] = [[0, 0, 1, 0], [1, 0, .1, 0], [2, -0.1, 0.1, 0], [3, 0.1, 0.1, 0]]
        seed_geom['umb_elems'] = [[0, 0, 1], [1, 1, 2], [2, 1, 3]]
        seed_geom['elem_up'] = [[0, 0, 0], [1, 0, 0], [1, 0, 0]]
        seed_geom['elem_down'] = [[2, 1, 2], [0, 0, 0], [0, 0, 0]]
        data = [[-0.1, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.0, 0.3], [0.1, 0.1, 0.1], [-0.2, -0.2, -0.2]]
        chorion_geom = placentagen.grow_chorionic_surface(90 * np.pi / 180, 45 * np.pi / 180, 0.5, 0.1, 1,
                                                 1, 1, 1, data, seed_geom, 'surface')
        self.assertTrue(chorion_geom['nodes'][6][3], 0.48527182)
