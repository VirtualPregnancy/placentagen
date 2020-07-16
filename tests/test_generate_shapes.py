from unittest import TestCase

import placentagen
import unittest
import numpy as np
import os

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

class Test_gen_rect_cover_ellipsoid(TestCase):

    def test_rect_el_num(self):
        mesh_el = placentagen.gen_rect_cover_ellipsoid(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.assertTrue(mesh_el['total_elems'] == 4)

    def test_rect_el_val(self):
        mesh_el = placentagen.gen_rect_cover_ellipsoid(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.assertTrue(mesh_el['elems'][0][5] == 9)

    def test_rect_node_num(self):
        mesh_el = placentagen.gen_rect_cover_ellipsoid(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.assertTrue(mesh_el['total_nodes'] == 18)

    def test_rect_node_val(self):
        mesh_el = placentagen.gen_rect_cover_ellipsoid(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        self.assertTrue(np.isclose(mesh_el['nodes'][14][2],0.5))

class Test_tet_mesh(TestCase):

      def test_darcy_node(self):
          volume=5
          thickness=2.1
          ellipticity=1
          mesh_node = placentagen.gen_ellip_mesh_tet(volume, thickness, ellipticity,28)
          self.assertTrue(np.isclose(mesh_node['nodes'][0,0],0))
          self.assertTrue(np.isclose(mesh_node['nodes'][0,1],0))
          self.assertTrue(np.isclose(mesh_node['nodes'][0,2],-1.05))
      def test_darcy_el(self):
          volume=5
          thickness=2.1
          ellipticity=1
          mesh_node = placentagen.gen_ellip_mesh_tet(volume, thickness, ellipticity,28)
          self.assertTrue(np.isclose(mesh_node['elems'][0,0],4))
          self.assertTrue(np.isclose(mesh_node['elems'][0,1],7))
          self.assertTrue(np.isclose(mesh_node['elems'][0,2],6))
          self.assertTrue(np.isclose(mesh_node['elems'][0,3],3))
 
      def test_darcy_el_node_array(self):
          volume=5
          thickness=2.1
          ellipticity=1
          mesh_node = placentagen.gen_ellip_mesh_tet(volume, thickness, ellipticity,28)
          self.assertTrue(mesh_node['element_array'][0] == 1)
          self.assertTrue(mesh_node['node_array'][0] == 1)

class Test_rect_node(TestCase):
      def test_rectangular_node(self):
          rect_node = placentagen.gen_rectangular_node(1, 1, 1, 2, 2, 2)
          self.assertTrue(rect_node[0,0],-0.5)
          self.assertTrue(rect_node[0,1],-0.5)
          self.assertTrue(rect_node[0,2],-0.5)

class Test_cube_mesh_con(TestCase):
      def test_mesh_con(self):
          mesh_connectivity = placentagen.cube_mesh_connectivity(2,2,2)
          self.assertTrue((mesh_connectivity[0] == [0,0,1,2,3,4,5,6,7]).all())        
          
class Test_pl_mesh_linear(TestCase):
      def test_placental_node(self):
          pl_mesh = placentagen.gen_3d_ellipsoid(1,1,1,5,1,1,1)
          self.assertTrue(np.isclose(pl_mesh['placental_node_coor'][0,0],-0.892062058076))
          self.assertTrue(np.isclose(pl_mesh['placental_node_coor'][0,1],-0.892062058076))
          self.assertTrue(np.isclose(pl_mesh['placental_node_coor'][0,2],-0.28867513))

      def test_placental_el(self):
          pl_mesh = placentagen.gen_3d_ellipsoid(1,1,1,5,1,1,1)
          self.assertTrue((pl_mesh['placental_el_con'][0] == [0,0,1,2,3,4,5,6,7]).all())

      def test_el_node_array(self):
          pl_mesh = placentagen.gen_3d_ellipsoid(1,1,1,5,1,1,1)

          self.assertTrue(pl_mesh['element_array'][0] == 1)
          self.assertTrue(pl_mesh['node_array'][0] == 1)
          self.assertTrue(pl_mesh['node_array'][7] == 8)

class Test_pl_mesh_qua(TestCase):
      def test_placental_node_q(self):
          pl_mesh = placentagen.gen_3d_ellipsoid(1,1,1,5,1,1,2)
          self.assertTrue(np.isclose(pl_mesh['placental_node_coor'][0,0],-0.892062058076))
          self.assertTrue(np.isclose(pl_mesh['placental_node_coor'][0,1],-0.892062058076))
          self.assertTrue(np.isclose(pl_mesh['placental_node_coor'][0,2],-0.288675134595))
      def test_placental_el_q(self):
          pl_mesh = placentagen.gen_3d_ellipsoid(1,1,1,5,1,1,2)
          self.assertTrue((pl_mesh['placental_el_con'][0] == [0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]).all())


class Test_surface_node(TestCase):
      def test_surf_node(self):
          surfacenode = placentagen.identify_surface_node_quad(1,1,1)
          self.assertTrue((surfacenode == [1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27]).all())

class Test_half_ellipse(TestCase):
      def test_half_ellipse(self):
          pl_mesh = placentagen.gen_half_ellipsoid_structured(1,(np.pi*4/3),2,1,0.75,0.85,1,False)
          self.assertTrue(np.isclose(pl_mesh['nodes'][0][1],-0.5303300858899107))
          self.assertTrue(np.isclose(pl_mesh['nodes'][0][2], 0.5303300858899107))
          self.assertTrue(np.isclose(pl_mesh['nodes'][0][3], 0.0))

#TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'Testdata/stem_xy.txt')
#class Test_vessel_node(TestCase):
#      def test_vs_node(self):
          #volume=5
          #thickness=2.1
          #ellipticity=1
          #ellipsoid_coor=np.array([[-58.1785453,-58.1785453,6.37972047],[-55.26961804,-59.57977309, 6.53337577],
          ## [-52.36069077,-60.87933747,-6.6758829 ]])
          #surfacenode=np.array([1,2,3])
          #v_node = placentagen.identify_vessel_node(ellipsoid_coor,surfacenode,TESTDATA_FILENAME,volume,thickness,
          ## ellipticity)
          #self.assertTrue(v_node['spiral_array'] == 2)
          #self.assertTrue(v_node['decidual_array'] == 1)
          #self.assertTrue((v_node['vesselnode'] == [2,1]).all())
          #self.assertTrue(v_node['surfnode_ex_vessel'] == 3)
if __name__ == '__main__':
    unittest.main()
