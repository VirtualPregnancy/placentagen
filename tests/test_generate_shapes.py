from unittest import TestCase

import numpy as np
import unittest
import placentagen


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
        
    def test_meshgrid_el(self):
        mesh_el  = placentagen.gen_rectangular_mesh(0,5,0,10,0,8,2,5,4,2.5,2,2)
        self.assertTrue(np.isclose(len(mesh_el['nodeOfelement']), 40))
        self.assertTrue(np.isclose(mesh_el['total_mesh_el'],40))
        
    def test_meshgrid_node(self):
        mesh_el  = placentagen.gen_rectangular_mesh(0,5,0,10,0,8,2,5,4,2.5,2,2)
        self.assertTrue(np.isclose(len(mesh_el['z_coor']),90))
        self.assertTrue(np.isclose(mesh_el['total_mesh_node'],90))

if __name__ == '__main__':
   unittest.main()
