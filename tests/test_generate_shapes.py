from unittest import TestCase

import placentagen
import unittest
import numpy as np


class Test_generate_data(TestCase):
    def test_data_in_ellipsoid(self):
        thickness = (3.0 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        datapoints = placentagen.equispaced_data_in_ellipsoid(1, 1.0, thickness, 1.0)
        array_test = np.isclose(datapoints, [0.0, 0.0, 0.0])
        self.assertTrue(array_test.all)

    def test_data_on_ellipsoid(self):
        thickness = (3.0 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0
        datapoints = placentagen.uniform_data_on_ellipsoid(3, 1.0, thickness, 1.0, 0)
        array_test = np.isclose(datapoints[1][:], [0.57526684, -0.14461422, 0.18163017])
        self.assertTrue(array_test.all)

class Test_gen_rectangular_mesh(TestCase):

    def test_rect_el_num(self):
        mesh_el = placentagen.gen_rectangular_mesh(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.assertTrue(mesh_el['total_elems'] == 4)

    def test_rect_el_val(self):
        mesh_el = placentagen.gen_rectangular_mesh(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.assertTrue(mesh_el['elems'][0][5] == 9)

    def test_rect_node_num(self):
        mesh_el = placentagen.gen_rectangular_mesh(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.assertTrue(mesh_el['total_nodes'] == 18)

    def test_rect_node_val(self):
        mesh_el = placentagen.gen_rectangular_mesh(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.assertTrue(np.isclose(mesh_el['nodes'][14][2],0.5))


class Test_darcy_mesh(TestCase):

      def test_darcy_node(self):
          volume=5
          thickness=2.1
          ellipticity=1
          mesh_node = placentagen.gen_mesh_darcy(volume, thickness, ellipticity,28)
          self.assertTrue(np.isclose(mesh_node['nodes'][0,0],0))
          self.assertTrue(np.isclose(mesh_node['nodes'][0,1],-1.06621809))
          self.assertTrue(np.isclose(mesh_node['nodes'][0,2],0))
      def test_darcy_el(self):
          volume=5
          thickness=2.1
          ellipticity=1
          mesh_node = placentagen.gen_mesh_darcy(volume, thickness, ellipticity,28)
          self.assertTrue(np.isclose(mesh_node['elems'][0,0],4))
          self.assertTrue(np.isclose(mesh_node['elems'][0,1],5))
          self.assertTrue(np.isclose(mesh_node['elems'][0,2],7))
          self.assertTrue(np.isclose(mesh_node['elems'][0,3],2))
 
      def test_darcy_el_node_array(self):
          volume=5
          thickness=2.1
          ellipticity=1
          mesh_node = placentagen.gen_mesh_darcy(volume, thickness, ellipticity,28)
          self.assertTrue(mesh_node['element_array'][0] == 1)
          self.assertTrue(mesh_node['node_array'][0] == 1)

if __name__ == '__main__':
    unittest.main()
