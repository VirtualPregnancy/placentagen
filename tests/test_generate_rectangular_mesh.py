import os
from unittest import TestCase
import unittest
import numpy as np
import placentagen as pg
from parameter import *

class Test_meshgrid_gen(TestCase):
        
    def test_meshgrid_el(self):
        mesh_el  = pg.generate_rectangular_mesh(x_min,x_max,y_min,y_max,z_min,z_max,nel_x,nel_y,nel_z,x_width,y_width,z_width)
        self.assertTrue(np.isclose(len(mesh_el['nodeOfelement']), 40))
        self.assertTrue(np.isclose(mesh_el['total_mesh_el'],40))
        
    def test_meshgrid_node(self):
        mesh_el  = pg.generate_rectangular_mesh(x_min,x_max,y_min,y_max,z_min,z_max,nel_x,nel_y,nel_z,x_width,y_width,z_width)
        self.assertTrue(np.isclose(len(mesh_el['z_coor']),90))
        self.assertTrue(np.isclose(mesh_el['total_mesh_node'],90))

if __name__ == '__main__':
   unittest.main()
