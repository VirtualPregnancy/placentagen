from unittest import TestCase

import numpy as np
import unittest
import placentagen
import os

class Test_create_trees(TestCase):
    def test_umilical_node(self):
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        data_input = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        datapoints = placentagen.umbilical_seed_geometry(1.0, thickness, 1.0, 0.0, 0.0, 0.1, 20.0, data_input)
        self.assertTrue(datapoints['nodes'][5][3], 0.61833222)

    def test_umilical_node(self):
        thickness = (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        data_input = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        datapoints = placentagen.umbilical_seed_geometry(1.0, thickness, 1.0, 0.0, 0.0, 0.1, 20.0, data_input)
        self.assertTrue(datapoints['elems'][2][2], 4)


class Test_grow_trees(TestCase):
    def test_grow_node(self):
        seed_geom = {}
        seed_geom['nodes'] = [[0, 0, 1, 0], [1, 0, .1, 0], [2, -0.1, 0.1, 0], [3, 0.1, 0.1, 0]]
        seed_geom['elems'] = [[0, 0, 1], [1, 1, 2], [2, 1, 3]]
        seed_geom['elem_up'] = [[0, 0, 0], [1, 0, 0], [1, 0, 0]]
        seed_geom['elem_down'] = [[2, 1, 2], [0, 0, 0], [0, 0, 0]]
        data = [[-0.1, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.0, 0.3], [0.1, 0.1, 0.1], [-0.2, -0.2, -0.2]]
        chorion_geom = placentagen.grow_chorionic_surface(90 * np.pi / 180, 45 * np.pi / 180, 0.5, 0.1, 1,
                                                          1, 1, 1, data, seed_geom, 'surface')
        self.assertTrue(chorion_geom['nodes'][6][3], 0.48527182)

    def test_grow_node_alt(self):
        seed_geom = {}
        seed_geom['nodes'] = [[0, 0, 1, 0], [1, 0, .1, 0], [2, -0.1, 0.1, 0], [3, 0.1, 0.1, 0]]
        seed_geom['elems'] = [[0, 0, 1], [1, 1, 2], [2, 1, 3]]
        seed_geom['elem_up'] = [[0, 0, 0], [1, 0, 0], [1, 0, 0]]
        seed_geom['elem_down'] = [[2, 1, 2], [0, 0, 0], [0, 0, 0]]
        data = placentagen.uniform_data_on_ellipsoid(5,1,1,1,0)
        geom = placentagen.grow_large_tree(90 * np.pi / 180, 45 * np.pi / 180, 0.5, 0.1, 1,
                                                              1, 1, 1, data, seed_geom)
        self.assertTrue(geom['nodes'][6][3], 0.22301891)


class Test_add_villi(TestCase):
    def test_add_villi(self):
        from_elem = 0
        initial_geom = {}
        initial_geom['nodes'] = [[0, 0.0, 0.0, 0.0], [1, 0.0, 0.5, 0.0], [2, 0.0, 1.0, 0.0]]
        initial_geom['elems'] = [[0, 0, 1], [1, 1, 2]]
        initial_geom['elem_up'] = [[0, 0, 0], [1, 0, 0]]
        initial_geom['elem_down'] = [[1, 1, 0], [0, 0, 0]]
        chorion_and_stem = placentagen.add_stem_villi(initial_geom, from_elem, 0.2)
        self.assertTrue(chorion_and_stem['nodes'][3][3], -0.2)
