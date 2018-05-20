
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
        thickness =  2.1  # mm
        ellipticity = 1.00  # no units
        volume=5
       
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = np.array([[-0.5, -0.5, -1.5],[ 0.5, -0.5,-1.5],[-0.5,  0.5 ,-1.5],[ 0.5 , 0.5, -1.5],[-0.5 ,-0.5, -0.5],[ 0.5 ,-0.5 ,-0.5],[-0.5 , 0.5 ,-0.5],[ 0.5 , 0.5 ,-0.5],[-0.5, -0.5 , 0.5],[ 0.5, -0.5 , 0.5],[-0.5  ,0.5 , 0.5],[ 0.5 , 0.5  ,0.5]])
        rectangular_mesh['elems'] = [[ 0,  0,  1,  2,  3,  4, 5, 6, 7],[1,4,5,6,7,8,9,10,11]]
        rectangular_mesh['total_elems'] = 2
        branch_elems={}
        branch_elems['elems']=[[0 ,0, 1]]
        branch_nodes={}
        branch_nodes['nodes']=np.array([[ 0.,0.,0., -1., 2.,0.,0.],[ 1.,0.,0.,-0.5 ,2.,0.,0.]])
        branch_radius=[0.1]
        br_vol_in_grid=placentagen.cal_br_vol_samp_grid(rectangular_mesh,  branch_nodes['nodes'], branch_elems['elems'],branch_radius, volume, thickness,ellipticity, 0)
        self.assertTrue(np.isclose(br_vol_in_grid['br_vol_in_grid'][0],   0.01396263))
        self.assertTrue(np.isclose(br_vol_in_grid['br_vol_in_grid'][1],    0.00174533))

    
    def test_br_diameter_sampling_grid(self):
        thickness =  2.1  # mm
        ellipticity = 1.00  # no units
        volume=5
       
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = np.array([[-0.5, -0.5, -1.5],[ 0.5, -0.5,-1.5],[-0.5,  0.5 ,-1.5],[ 0.5 , 0.5, -1.5],[-0.5 ,-0.5, -0.5],[ 0.5 ,-0.5 ,-0.5],[-0.5 , 0.5 ,-0.5],[ 0.5 , 0.5 ,-0.5],[-0.5, -0.5 , 0.5],[ 0.5, -0.5 , 0.5],[-0.5  ,0.5 , 0.5],[ 0.5 , 0.5  ,0.5]])
        rectangular_mesh['elems'] = [[ 0,  0,  1,  2,  3,  4, 5, 6, 7],[1,4,5,6,7,8,9,10,11]]
        rectangular_mesh['total_elems'] = 2
        branch_elems={}
        branch_elems['elems']=[[0 ,0, 1]]
        branch_nodes={}
        branch_nodes['nodes']=np.array([[ 0.,0.,0., -1., 2.,0.,0.],[ 1.,0.,0.,-0.5 ,2.,0.,0.]])
        branch_radius=[0.1]
        br_vol_in_grid=placentagen.cal_br_vol_samp_grid(rectangular_mesh,  branch_nodes['nodes'], branch_elems['elems'],branch_radius, volume, thickness,ellipticity, 0)
        
        self.assertTrue(np.isclose(br_vol_in_grid['br_diameter_in_grid'][0],  0.00279253))
        self.assertTrue(np.isclose(br_vol_in_grid['br_diameter_in_grid'][1],  0.00034907))




class Test_terminals_in_sampling_grid_fast(TestCase):
        
    def test_terminals_in_grid_present(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br={}
        term_br['terminal_nodes']=[3]
        term_br['total_terminals']=1
        rectangular_mesh = {}
        rectangular_mesh['nodes'] =np.array( [[ 0.,  0.,  0.],[ 1.,  0. , 0.],[ 0.,  1. , 0.],[ 1. , 1. , 0.],[ 0.,  0. , 1.],[ 1.,  0. , 1.],[ 0. , 1. , 1.],[ 1. , 1. , 1.]])
        rectangular_mesh['elems']=[[0, 0, 1, 2, 3, 4, 5, 6, 7]]
        term_grid =placentagen.terminals_in_sampling_grid_fast(rectangular_mesh, term_br, noddata['nodes'])
        self.assertTrue(term_grid['terminals_in_grid'][0] == 1)
        
     
    def test_terminal_elems_present(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br={}
        term_br['terminal_nodes']=[3]
        term_br['total_terminals']=1
        rectangular_mesh = {}
        rectangular_mesh['nodes'] =np.array( [[ 0.,  0.,  0.],[ 1.,  0. , 0.],[ 0.,  1. , 0.],[ 1. , 1. , 0.],[ 0.,  0. , 1.],[ 1.,  0. , 1.],[ 0. , 1. , 1.],[ 1. , 1. , 1.]])
        rectangular_mesh['elems']=[[0, 0, 1, 2, 3, 4, 5, 6, 7]]
        term_grid =placentagen.terminals_in_sampling_grid_fast(rectangular_mesh, term_br, noddata['nodes'])
        self.assertTrue(term_grid['terminal_elems'][0] == 0)#this zero does not mean branch are not located. it means samp_grid el 0


class Test_terminals_in_sampling_grid_general(TestCase):

    def test_terminals_in_grid_general_present(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br = {}
        term_br['terminal_nodes'] = [3]
        term_br['total_terminals'] = 1
        placenta_list = [7]
        rectangular_mesh = {}
        rectangular_mesh['elems'] = np.zeros((8, 9), dtype=int)
        rectangular_mesh['nodes'] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.], [0., 0., 1.], [1., 0., 1.],
                                     [0., 1., 1.], [1., 1., 1.]]
        rectangular_mesh['elems'][7] = [0, 0, 1, 2, 3, 4, 5, 6, 7]
        term_grid = placentagen.terminals_in_sampling_grid(rectangular_mesh, placenta_list, term_br, noddata['nodes'])
        self.assertTrue(term_grid['terminals_in_grid'][7] == 1)

    def test_terminals_elem_general_present(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br = {}
        term_br['terminal_nodes'] = [3]
        term_br['total_terminals'] = 1
        placenta_list = [7]
        rectangular_mesh = {}
        rectangular_mesh['elems'] = np.zeros((8, 9), dtype=int)
        rectangular_mesh['nodes'] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.], [0., 0., 1.], [1., 0., 1.],
                                     [0., 1., 1.], [1., 1., 1.]]
        rectangular_mesh['elems'][7] = [0, 0, 1, 2, 3, 4, 5, 6, 7]
        term_grid = placentagen.terminals_in_sampling_grid(rectangular_mesh, placenta_list, term_br, noddata['nodes'])
        self.assertTrue(term_grid['terminal_elems'][0] == 7)

    def test_terminals_in_grid_general_absent(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        eldata = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        term_br = placentagen.calc_terminal_branch(noddata['nodes'], eldata['elems'])
        placenta_list = [1]
        rectangular_mesh = {}
        rectangular_mesh['elems'] = np.zeros((8, 9), dtype=int)
        rectangular_mesh['nodes'] = [[0., -1., -1.], [1., -1., -1.], [0., 0., -1.], [1., 0., -1.], [0., -1., 0.],
                                     [1., -1., 0.], [0., 0., 0.], [1., 0., 0.]]
        rectangular_mesh['elems'][1] = [0, 0, 1, 2, 3, 4, 5, 6, 7]
        term_grid = placentagen.terminals_in_sampling_grid(rectangular_mesh, placenta_list, term_br, noddata['nodes'])
        self.assertTrue(
            np.sum(term_grid['terminals_in_grid']) == 0)  # all must be zero as could not locate any terminal br

    def test_terminals_elem_general_absent(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        eldata = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        term_br = placentagen.calc_terminal_branch(noddata['nodes'], eldata['elems'])
        placenta_list = [1]
        rectangular_mesh = {}
        rectangular_mesh['elems'] = np.zeros((8, 9), dtype=int)
        rectangular_mesh['nodes'] = [[0., -1., -1.], [1., -1., -1.], [0., 0., -1.], [1., 0., -1.], [0., -1., 0.],
                                     [1., -1., 0.], [0., 0., 0.], [1., 0., 0.]]
        rectangular_mesh['elems'][1] = [0, 0, 1, 2, 3, 4, 5, 6, 7]
        term_grid = placentagen.terminals_in_sampling_grid(rectangular_mesh, placenta_list, term_br, noddata['nodes'])
        self.assertTrue(
            np.sum(term_grid['terminal_elems']) == 0)  # all must be zero as could not locate any terminal br



class Test_terminals_villous_volume(TestCase):
        
    def test_terminals_vill_vol(self):

        num_int_gens = 3
        num_convolutes = 10
        len_int = 1.5 #mm
        rad_int = 0.03 #mm
        len_convolute = 3.0 #mm
        rad_convolute = 0.025 #mm
        term_vill_vol=placentagen.terminal_villous_volume(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute)
        self.assertTrue(np.isclose(term_vill_vol,1.77657064561))
      

class Test_tissue_volume_gr(TestCase):
        
    def test_tissue_vol(self):

        tissue_vol=placentagen.tissue_vol_in_samp_gr(1.77, 0.014,2)
        
        self.assertTrue(np.isclose(tissue_vol,3.554))




class Test_terminals_villous_diameter(TestCase):
        
    def test_terminals_vill_diameter(self):

        num_int_gens = 3
        num_convolutes = 10
        len_int = 1.5 #mm
        rad_int = 0.03 #mm
        len_convolute = 3.0 #mm
        rad_convolute = 0.025 #mm
        term_vill_diameter=placentagen.terminal_villous_diameter(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute)
        
        self.assertTrue(np.isclose(term_vill_diameter,0.090100877305))


class Test_weighted_diameter(TestCase):
        
    def test_wt_diameter(self):
        br_diameter_in_grid=np.array([0.0023])
        tissue_vol=np.array([3.554])
        wt_D=placentagen.weighted_diameter_in_samp_gr(0.09,br_diameter_in_grid,2,tissue_vol)
        self.assertTrue(np.isclose(wt_D,0.05129432))
        
      
if __name__ == '__main__':
   unittest.main()





    
    

