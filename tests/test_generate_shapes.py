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
          volume=50
          thickness=5
          ellipticity=1
          rectangular_mesh={}
          rectangular_mesh['nodes']=np.array([[-3.,-3.,-3.],[-1.,-3.,-3.],[1.,-3.,-3.],[3.,-3.,-3.],[-3.,-1.,-3.],[-1.,-1.,-3.],[1.,-1.,-3.],[3.,-1.,-3.],[-3.,1.,-3.],[-1.,1.,-3.],[3.,1.,-3.],[-3.,3.,-3.],[-1.,3.,-3.],[1.,3.,-3.],[3.,3.,-3.],[-3.,-3.,-1.],[-1.,-3.,-1.],[1.,-3.,-1.],[3.,-3.,-1.],[-3.,-1.,-1.],[-1.,-1.,-1.],[1.,-1.,-1.],[3.,-1.,-1.],[-3.,1.,-1.],[-1.,1.,-1.],[1.,1.,-1.],[3.,1.,-1.],[-3.,3.,-1.], [-1.,  3., -1.], [ 1. , 3., -1.], [ 3. , 3., -1.], [-3., -3. , 1.], [-1. ,-3., 1.],[ 1. ,-3. , 1.], [ 3., -3. , 1.], [-3. ,-1.,  1.], [-1. ,-1.,  1.], [ 1., -1.,  1.], [ 3. ,-1. , 1.], [-3. , 1., 1.], [-1. , 1.,1.],[1.,1.,1.],[3.,1.,1.],[-3.,3.,1.],[-1.,3.,1.],[1.,3.,1.],[3.,3.,1.],[-3.,-3.,3.],[-1.,-3.,3.],[1.,-3.,3.],[3.,-3.,3.],[-3.,-1.,3.],[-1.,-1.,3.],[1.,-1.,3.],[3.,-1.,3.],[-3.,1.,3.],[-1.,1.,3.],[ 1.,1., 3.],[3.,1.,3.],[-3.,3.,3.],[-1.,3.,3.],[1.,3.,3.],[3.,3.,3.]]
)
          mesh_node = placentagen.gen_mesh_darcy(rectangular_mesh,volume, thickness, ellipticity)
          self.assertTrue(mesh_node['nodes'][0,0] == -1)
          self.assertTrue(mesh_node['nodes'][0,1] == -1)
          self.assertTrue(np.isclose(mesh_node['nodes'][0,2],-1.90578229))

      def test_darcy_el(self):
          volume=50
          thickness=5
          ellipticity=1
          rectangular_mesh={}
          rectangular_mesh['nodes']=np.array([[-3.,-3.,-3.],[-1.,-3.,-3.],[1.,-3.,-3.],[3.,-3.,-3.],[-3.,-1.,-3.],[-1.,-1.,-3.],[1.,-1.,-3.],[3.,-1.,-3.],[-3.,1.,-3.],[-1.,1.,-3.],[3.,1.,-3.],[-3.,3.,-3.],[-1.,3.,-3.],[1.,3.,-3.],[3.,3.,-3.],[-3.,-3.,-1.],[-1.,-3.,-1.],[1.,-3.,-1.],[3.,-3.,-1.],[-3.,-1.,-1.],[-1.,-1.,-1.],[1.,-1.,-1.],[3.,-1.,-1.],[-3.,1.,-1.],[-1.,1.,-1.],[1.,1.,-1.],[3.,1.,-1.],[-3.,3.,-1.], [-1.,  3., -1.], [ 1. , 3., -1.], [ 3. , 3., -1.], [-3., -3. , 1.], [-1. ,-3., 1.],[ 1. ,-3. , 1.], [ 3., -3. , 1.], [-3. ,-1.,  1.], [-1. ,-1.,  1.], [ 1., -1.,  1.], [ 3. ,-1. , 1.], [-3. , 1., 1.], [-1. , 1.,1.],[1.,1.,1.],[3.,1.,1.],[-3.,3.,1.],[-1.,3.,1.],[1.,3.,1.],[3.,3.,1.],[-3.,-3.,3.],[-1.,-3.,3.],[1.,-3.,3.],[3.,-3.,3.],[-3.,-1.,3.],[-1.,-1.,3.],[1.,-1.,3.],[3.,-1.,3.],[-3.,1.,3.],[-1.,1.,3.],[ 1.,1., 3.],[3.,1.,3.],[-3.,3.,3.],[-1.,3.,3.],[1.,3.,3.],[3.,3.,3.]]
)
          mesh_node = placentagen.gen_mesh_darcy(rectangular_mesh,volume, thickness, ellipticity)
          self.assertTrue(mesh_node['elems'][0,0] == 7 )
          self.assertTrue(mesh_node['elems'][0,1] == 3 )
          self.assertTrue(mesh_node['elems'][0,2] == 2 )
          self.assertTrue(mesh_node['elems'][0,3] == 1 )

      def test_darcy_el_node_array(self):
          volume=50
          thickness=5
          ellipticity=1
          rectangular_mesh={}
          rectangular_mesh['nodes']=np.array([[-3.,-3.,-3.],[-1.,-3.,-3.],[1.,-3.,-3.],[3.,-3.,-3.],[-3.,-1.,-3.],[-1.,-1.,-3.],[1.,-1.,-3.],[3.,-1.,-3.],[-3.,1.,-3.],[-1.,1.,-3.],[3.,1.,-3.],[-3.,3.,-3.],[-1.,3.,-3.],[1.,3.,-3.],[3.,3.,-3.],[-3.,-3.,-1.],[-1.,-3.,-1.],[1.,-3.,-1.],[3.,-3.,-1.],[-3.,-1.,-1.],[-1.,-1.,-1.],[1.,-1.,-1.],[3.,-1.,-1.],[-3.,1.,-1.],[-1.,1.,-1.],[1.,1.,-1.],[3.,1.,-1.],[-3.,3.,-1.], [-1.,  3., -1.], [ 1. , 3., -1.], [ 3. , 3., -1.], [-3., -3. , 1.], [-1. ,-3., 1.],[ 1. ,-3. , 1.], [ 3., -3. , 1.], [-3. ,-1.,  1.], [-1. ,-1.,  1.], [ 1., -1.,  1.], [ 3. ,-1. , 1.], [-3. , 1., 1.], [-1. , 1.,1.],[1.,1.,1.],[3.,1.,1.],[-3.,3.,1.],[-1.,3.,1.],[1.,3.,1.],[3.,3.,1.],[-3.,-3.,3.],[-1.,-3.,3.],[1.,-3.,3.],[3.,-3.,3.],[-3.,-1.,3.],[-1.,-1.,3.],[1.,-1.,3.],[3.,-1.,3.],[-3.,1.,3.],[-1.,1.,3.],[ 1.,1., 3.],[3.,1.,3.],[-3.,3.,3.],[-1.,3.,3.],[1.,3.,3.],[3.,3.,3.]]
)
          mesh_node = placentagen.gen_mesh_darcy(rectangular_mesh,volume, thickness, ellipticity)
          self.assertTrue(mesh_node['element_array'][0] == 1)
          self.assertTrue(mesh_node['node_array'][7] == 8)

if __name__ == '__main__':
    unittest.main()
