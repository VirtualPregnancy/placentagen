
from unittest import TestCase

import numpy as np
import unittest
import placentagen
import os
TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exnode')
TESTDATA_FILENAME1 = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exelem')

class Test_Terminal_Br(TestCase):
        
    def test_terminal_br(self):
        eldata   = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br  = placentagen.calc_terminal_branch(noddata['nodes'],eldata['elems'])
        print(term_br['terminal_elems'])
        self.assertTrue(term_br['terminal_elems'][0] == 1)
        
    def test_terminal_br_total(self):
        eldata   = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br  = placentagen.calc_terminal_branch(noddata['nodes'],eldata['elems'])
        self.assertTrue(term_br['total_terminals'] == 2)

      
class test_pl_vol_in_grid(TestCase):
        
    def test_pl_vol_margin(self):
        thickness =  (3.0 * 1 / (4.0 * np.pi)) ** (1.0 / 3.0) * 2.0  # mm
        ellipticity = 1.00  # no units
        spacing = 1.0  # mm
        volume=1
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = [[0., 0., 0.], [ thickness/2.0, 0., 0.],[0., thickness/2.0, 0.],[ thickness/2.0, thickness/2.0, 0.],[0., 0., thickness/2.0], [ thickness/2.0, 0., thickness/2.0],[0., thickness/2.0,thickness/2.0],[ thickness/2.0, thickness/2.0, thickness/2.0]]
        rectangular_mesh['elems'] = [[ 0,  0,  1,  2,  3,  4, 5, 6, 7]]
        rectangular_mesh['total_nodes'] =8
        rectangular_mesh['total_elems'] = 1
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 25)
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][0], 0.12485807941))
        self.assertTrue(abs(pl_vol['pl_vol_in_grid']-1./8.)/(1./8)<1e-2)#looking for less than 1% error in expected volume of 1/8

    def test_pl_vol_complete_inside(self):
        thickness =  2  # mm
        ellipticity = 1.6  # no units
        spacing = 0.5  # mm
        volume=5
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = [[-1., -1.5, -1.],[-0.5 ,-1.5, -1.],[-1., -1., -1.] ,[-0.5,-1., -1.],[-1.,-1.5,-0.5],[-0.5,-1.5,-0.5],[-1.,-1.,-0.5] ,[-0.5,-1.,-0.5]]
        rectangular_mesh['elems'] = [[0, 0, 1, 2, 3, 4, 5, 6, 7]]
        rectangular_mesh['total_nodes'] =8
        rectangular_mesh['total_elems'] = 1
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 25)
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][0], 0.0))
      

    def test_pl_vol_complete_outside(self):
        thickness =  2  # mm
        ellipticity = 1.6  # no units
        spacing = 0.5  # mm
        volume=5
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = [[-0.5,-0.5,-0.5],[ 0., -0.5,-0.5],[-0.5, 0.,-0.5],[ 0., 0. ,-0.5],[-0.5, -0.5 ,0. ],[ 0., -0.5 ,0.],[-0.5, 0.,0.],[0.,0.,0.]]
        rectangular_mesh['elems'] = [[0,  0,  1,  2,  3,  4, 5, 6, 7]]
        rectangular_mesh['total_nodes'] =8
        rectangular_mesh['total_elems'] = 1
        pl_vol=placentagen.ellipse_volume_to_grid(rectangular_mesh, volume, thickness, ellipticity, 0.125)
        self.assertTrue(np.isclose(pl_vol['pl_vol_in_grid'][0], spacing*spacing*spacing))
   


class test_br_vol_in_grid(TestCase):
        
    def test_br_vol_sampling_grid(self):
        thickness =  2  # mm
        ellipticity = 1.00  # no units
        volume=5
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = np.array([[-0.5,-0.5, -1. ],[ 0.5 ,-0.5 ,-1. ],[-0.5 , 0.5 ,-1. ],[ 0.5  ,0.5, -1. ],[-0.5, -0.5 , 0. ],[ 0.5 ,-0.5 , 0. ],[-0.5  ,0.5  ,0. ],[ 0.5,  0.5 , 0. ]])
        rectangular_mesh['elems'] = [[ 0,  0,  1,  2,  3,  4, 5, 6, 7]]
        rectangular_mesh['total_elems'] = 1
        p_vol={}
        p_vol['pl_vol_in_grid']=[ 0.926015574057]
        eldata={}
        eldata['elems']=[[0 ,0, 1]]
        nodedata={}
        nodedata['nodes']=[[ 0.,0.,0., -1., 2.,0.,0.],[ 1.,0.,0.,-0.5 ,2.,0.,0.]]
        br_vol_in_grid=placentagen.cal_br_vol_samp_grid(rectangular_mesh,eldata,nodedata,5,2,1,p_vol)
        self.assertTrue(np.isclose(br_vol_in_grid['total_vol_samp_gr'][0,0],  0.01570796))
        self.assertTrue(np.isclose(br_vol_in_grid['total_vol_samp_gr'][0,1],  0.00314159))


           
    def test_each_and_total_br_vol(self):
        thickness =  2  # mm
        ellipticity = 1.00  # no units
        volume=5
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = np.array([[-0.5,-0.5, -1. ],[ 0.5 ,-0.5 ,-1. ],[-0.5 , 0.5 ,-1. ],[ 0.5  ,0.5, -1. ],[-0.5, -0.5 , 0. ],[ 0.5 ,-0.5 , 0. ],[-0.5  ,0.5  ,0. ],[ 0.5,  0.5 , 0. ]])
        rectangular_mesh['elems'] = [[ 0,  0,  1,  2,  3,  4, 5, 6, 7]]
        rectangular_mesh['total_elems'] = 1
        p_vol={}
        p_vol['pl_vol_in_grid']=[ 0.926015574057]
        eldata={}
        eldata['elems']=[[0 ,0, 1]]
        nodedata={}
        nodedata['nodes']=[[ 0.,0.,0., -1., 2.,0.,0.],[ 1.,0.,0.,-0.5 ,2.,0.,0.]]
        br_vol_in_grid=placentagen.cal_br_vol_samp_grid(rectangular_mesh,eldata,nodedata,5,2,1,p_vol)
        self.assertTrue(np.isclose(br_vol_in_grid['total_br_vol'],  0.01570796))
        self.assertTrue(abs(br_vol_in_grid['total_br_vol']-0.0157)<1e-2)#error tolerance for the value when calculate manually using pi is 3.14
        self.assertTrue(np.isclose(br_vol_in_grid['vol_each_br'][0],  0.01570796))
        self.assertTrue(abs(br_vol_in_grid['vol_each_br'][0]-0.0157)<1e-2)#error tolerance for the value when calculate manually using pi is 3.14

    def test_br_num_samp_gr(self):
        thickness =  2  # mm
        ellipticity = 1.00  # no units
        volume=5
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = np.array([[-0.5,-0.5, -1. ],[ 0.5 ,-0.5 ,-1. ],[-0.5 , 0.5 ,-1. ],[ 0.5  ,0.5, -1. ],[-0.5, -0.5 , 0. ],[ 0.5 ,-0.5 , 0. ],[-0.5  ,0.5  ,0. ],[ 0.5,  0.5 , 0. ]])
        rectangular_mesh['elems'] = [[ 0, 0, 1, 2, 3, 4,5 ,6,7]]
        rectangular_mesh['total_elems'] = 1
        p_vol={}
        p_vol['pl_vol_in_grid']=[0.926015574057]
        eldata={}
        eldata['elems']=[[0 ,0, 1]]
        nodedata={}
        nodedata['nodes']=[[ 0.,0.,0., -1., 2.,0.,0.],[ 1.,0.,0.,-0.5 ,2.,0.,0.]]
        br_vol_in_grid=placentagen.cal_br_vol_samp_grid(rectangular_mesh,eldata,nodedata,5,2,1,p_vol)
        self.assertTrue(np.isclose(br_vol_in_grid['br_num_in_samp_gr'][0],1))
        
if __name__ == '__main__':
   unittest.main()





    
    

