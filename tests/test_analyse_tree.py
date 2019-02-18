
from unittest import TestCase

import numpy as np
import unittest
import placentagen
import os
TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exnode')
TESTDATA_FILENAME1 = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exelem')
TESTDDATA_FILENAME2 = os.path.join(os.path.dirname(__file__), 'Testdata/stem_xy.txt')

class test_arrange_by_br(TestCase):

    def arrange_by_branch_no_doubles(self):
        geom = {}
        geom['nodes'] = np.array(
            [[0., 0., 0., -1., 2., 0., 0.], [1., 0., 0., -0.5, 2., 0., 0.], [1., 0., -0.5, -0.5, 2., 0., 0.],
             [1., 0., 0.5, -0.5, 2., 0., 0.]])
        geom['elems'] = np.array([[0, 0, 1], [1, 1, 2], [2, 1, 3]], dtype=int)
        geom['radii'] = [0.1, 0.1, 0.1]
        geom['length'] = [0.5, 0.5, 0.5]
        geom['euclidean length'] = geom['length']
        elem_up = np.zeros((3, 3), dtype=int)
        elem_up[0, 0] = 0
        elem_up[0, 1] = 0
        elem_up[1, 0] = 1
        elem_up[1, 1] = 0
        elem_up[2, 0] = 1
        elem_up[2, 1] = 0

        order = [2, 1, 1]
        generation = [1, 2, 2]

        arranged = placentagen.arrange_by_branches(geom, elem_up, order, generation)
        self.assertTrue(len(arranged['branches']) == 3)

    def arrange_by_branch_doubles(self):
        geom = {}
        geom['nodes'] = np.array(
            [[0., 0., 0., -1., 2., 0., 0.], [0., 0., 0., -0.75, 2., 0., 0.],[1., 0., 0., -0.5, 2., 0., 0.], [1., 0., -0.5, -0.5, 2., 0., 0.],
             [1., 0., 0.5, -0.5, 2., 0., 0.]])
        geom['elems'] = np.array([[0, 0, 1], [1,1,2],[0, 2, 3], [2, 3, 4]], dtype=int)
        geom['radii'] = [0.1,0.1, 0.1, 0.1]
        geom['length'] = [0.25,0.25, 0.5, 0.5]
        geom['euclidean length'] = geom['length']
        elem_up = np.zeros((3, 3), dtype=int)
        elem_up[0, 0] = 0
        elem_up[0, 1] = 0
        elem_up[1, 0] = 1
        elem_up[1, 1] = 0
        elem_up[2, 0] = 1
        elem_up[2, 1] = 1
        elem_up[3, 0] = 1
        elem_up[3, 1] = 1

        order = [2, 2, 1, 1]
        generation = [1, 1, 2, 2]

        arranged = placentagen.arrange_by_branches(geom, elem_up, order, generation)
        self.assertTrue(len(arranged['branches']) == 4)

class test_arrange_strahler_order(TestCase):

    def test_simple_arrange_strahler_order(self):
        geom = {}
        geom['nodes'] = np.array(
            [[0., 0., 0., -1., 2., 0., 0.], [1., 0., 0., -0.5, 2., 0., 0.], [1., 0., -0.5, -0.5, 2., 0., 0.],
             [1., 0., 0.5, -0.5, 2., 0., 0.]])
        geom['elems'] = np.array([[0, 0, 1], [1, 1, 2], [2, 1, 3]], dtype=int)
        geom['radii'] = [0.1, 0.1, 0.1]
        geom['length'] = [0.5, 0.5, 0.5]
        geom['euclidean length'] = geom['length']

        arranged = placentagen.arrange_by_strahler_order(geom, 1, [0, 0, 0])
        self.assertTrue((arranged['elems'][0] == np.array([0,0,1])).all)

    def test_peturbed_arrange_strahler_order(self):
        geom = {}
        geom['nodes'] = np.array(
             [[0., 0., 0., -1., 2., 0., 0.], [1., 0., 0., -0.5, 2., 0., 0.], [1., 0., -0.5, -0.5, 2., 0., 0.],
             [1., 0., 0.5, -0.5, 2., 0., 0.]])
        geom['elems'] = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 3]], dtype=int)
        geom['radii'] = [0.1, 0.1, 0.1]
        geom['length'] = [0.5, 0.5, 0.5]
        geom['euclidean length'] = geom['length']

        arranged = placentagen.arrange_by_strahler_order(geom, 1, [0, 0, 0])
        print(np.array([0, 0, 1]))
        self.assertTrue((arranged['elems'][0] == np.array([0, 0, 1])).all)

    def test_definlet_arrange_strahler_order(self):
        geom = {}
        geom['nodes'] = np.array(
            [[0., 0., 0., -1., 2., 0., 0.], [1., 0., 0., -0.5, 2., 0., 0.], [1., 0., -0.5, -0.5, 2., 0., 0.],
            [1., 0., 0.5, -0.5, 2., 0., 0.]])
        geom['elems'] = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 3]], dtype=int)
        geom['radii'] = [0.1, 0.1, 0.1]
        geom['length'] = [0.5, 0.5, 0.5]
        geom['euclidean length'] = geom['length']

        arranged = placentagen.arrange_by_strahler_order(geom, 0, [0., 0., -1.0])
        self.assertTrue((arranged['elems'][0] == np.array([0, 0, 1])).all)

class test_evaluate_orders(TestCase):

    def test_evaluate_orders(self):
        geom = {}
        geom['nodes'] = np.array(
            [[0., 0., 0., -1., 2., 0., 0.], [1., 0., 0., -0.5, 2., 0., 0.], [1., 0., -0.5, -0.5, 2., 0., 0.],
             [1., 0., 0.5, -0.5, 2., 0., 0.]])
        geom['elems'] = np.array([[0, 0, 1], [1, 1, 2], [2, 1, 3]], dtype=int)
        orders = placentagen.evaluate_orders(geom['nodes'],  geom['elems'])

        self.assertTrue(orders['strahler'][0] == 2 and orders['generation'][0] == 1 )

class test_summary_statistics(TestCase):

    def test_generation_summary_statistics(self):
        geom = {}
        geom['nodes'] = np.array(
            [[0., 0., 0., -1., 2., 0., 0.], [1., 0., 0., -0.5, 2., 0., 0.], [1., 0., -0.3, -0.5, 2., 0., 0.],
             [1., 0., 0.5, -0.5, 2., 0., 0.]])
        geom['elems'] = np.array([[0, 1, 2], [0, 0, 1], [0, 1, 3]], dtype=int)
        geom['radii'] = [0.1, 0.1, 0.05]
        geom['length'] = [0.5, 0.5, 0.5]
        geom['euclidean length'] = geom['length']
        geom['branch angles'] = [0.1, 0.5, 0.4]
        geom['diam_ratio'] = [0.1, 0.3, 0.4]
        geom['length_ratio'] = [0.1, 1.0, 0.6]
        elem_down = np.zeros((3, 3), dtype=int)
        elem_down[0, 0] = 2
        elem_down[0, 1] = 1
        elem_down[0, 2] = 2
        elem_down[1, 0] = 0
        elem_down[2, 0] = 0

        orders = {}
        orders['strahler'] = np.array([2, 1, 1], dtype=int)
        orders['generation'] = np.array([1, 2, 2], dtype=int)

        major_minor_results = placentagen.major_minor(geom, elem_down)

        arranged = placentagen.generation_summary_statistics(geom, orders, major_minor_results)
        self.assertTrue(np.isclose(np.array(arranged[2][0:19]), np.array(
            [-1., 3., 0.5, 0., 0.08333333, 0.02357023, 0.5, 0., 6.66666667, 2.3570226, 1., 0., 0.33333333, 0.16996732,
             0.4, 0., 0.5, 0., 0.56666667])).all)

class test_terminal_br(TestCase):
        
    def test_terminal_br(self):
        eldata   = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        term_br  = placentagen.calc_terminal_branch(noddata['nodes'],eldata['elems'])
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
        self.assertTrue(np.sum(term_grid['terminals_in_grid']) == 0)  # all must be zero as could not locate any terminal br

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
        self.assertTrue(np.sum(term_grid['terminal_elems']) == 0)  # all must be zero as could not locate any terminal br

#class Test_terminals_villous_volume(TestCase):
#
#    def test_terminals_vill_vol(self):#
#
#        num_int_gens = 3
#        num_convolutes = 10
#        len_int = 1.5 #mm
#        rad_int = 0.03 #mm
#        len_convolute = 3.0 #mm
#        rad_convolute = 0.025 #mm
#        term_vill_vol=placentagen.terminal_villous_volume(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute)
#        self.assertTrue(np.isclose(term_vill_vol,1.77657064561))

class Test_tissue_volume_gr(TestCase):
        
    def test_tissue_vol(self):
        tissue_vol=placentagen.tissue_vol_in_samp_gr(0.444, 0.008)   
        self.assertTrue(np.isclose(tissue_vol,0.452))

#class Test_terminals_villous_diameter(TestCase):
#
#    def test_terminals_vill_diameter(self):
#
#        num_int_gens = 3
 #       num_convolutes = 10
#        len_int = 1.5 #mm
##        rad_int = 0.03 #mm
 #       len_convolute = 3.0 #mm
 #       rad_convolute = 0.025 #mm
 #       term_vill_diameter=placentagen.terminal_villous_diameter(num_int_gens,num_convolutes,len_int,rad_int,len_convolute,rad_convolute)
 #       self.assertTrue(np.isclose(term_vill_diameter,0.090100877305))

class Test_conductivity_samp_gr(TestCase):
        
    def test_conductivity(self):

        vol_frac= [0.72401065]
        weighted_diameter=[0.17988357]
        non_empties=[0]
        conductivity=placentagen.conductivity_samp_gr(vol_frac,weighted_diameter,non_empties)
        self.assertTrue(np.isclose(conductivity, 7.20937313e-06))

#class Test_vol_frac_samp_gr(TestCase):
 #
 #   def test_volume_fraction(self):
 #       tissue_vol=[0.453]
 #       placental_volume={}
 #       placental_volume['non_empty_rects']=[0]
 #      placental_volume['pl_vol_in_grid']=[0.625]
 #       vol_frac=placentagen.vol_frac_in_samp_gr(tissue_vol,placental_volume)
 #       self.assertTrue(np.isclose(vol_frac, 0.7248))

class Test_term_vol_grid(TestCase):
        
    def test_terminal_vol_grid(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = np.array([[-2., -2. ,-2.],[ 0. ,-2. ,-2.],[ 2. ,-2. ,-2.],[-2. , 0., -2.],[ 0. , 0. ,-2.],[ 2. , 0. ,-2.],[-2. ,-2. , 0.],[ 0. ,-2. , 0.],[ 2. ,-2. , 0.],[-2. , 0. ,0.],[ 0. , 0.,  0.],[ 2.,  0. , 0.],[-2. ,-2. , 2.],[ 0. ,-2. , 2.],[ 2., -2.,  2.],[-2. , 0. , 2.],[ 0.,  0. , 2.],[ 2. , 0.,  2.]])
        rectangular_mesh['elems'] = [[ 0,0,1,3,4,6,7,9,10],[ 1,  1,2,4,5,7,8,10,11],[2,6,7,9,10,12,13,15,16],[4,7,8,10,11,13,14,16,17]]
        rectangular_mesh['total_elems'] = 4
        terminal_list={}
        terminal_list['total_terminals']=1
        terminal_list['terminal_nodes']=[2]
        volume=5
        thickness=2.1
        ellipticity=1
        term_total_vol=0.04#artificial value to match with a smaller ellipsoid
        term_tissue_vol=1.77657064561
        term_tissue_diam=0.090100877305
        term_vol=placentagen.terminal_volume_to_grid(rectangular_mesh, terminal_list, noddata['nodes'],volume, thickness,ellipticity,term_total_vol,term_tissue_vol, term_tissue_diam)
        self.assertTrue(np.isclose(term_vol['term_vol_in_grid'][0],0.44414266))

    def test_terminal_diam_grid(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        rectangular_mesh = {}
        rectangular_mesh['nodes'] = np.array([[-2., -2. ,-2.],[ 0. ,-2. ,-2.],[ 2. ,-2. ,-2.],[-2. , 0., -2.],[ 0. , 0. ,-2.],[ 2. , 0. ,-2.],[-2. ,-2. , 0.],[ 0. ,-2. , 0.],[ 2. ,-2. , 0.],[-2. , 0. ,0.],[ 0. , 0.,  0.],[ 2.,  0. , 0.],[-2. ,-2. , 2.],[ 0. ,-2. , 2.],[ 2., -2.,  2.],[-2. , 0. , 2.],[ 0.,  0. , 2.],[ 2. , 0.,  2.]])
        rectangular_mesh['elems'] = [[ 0,0,1,3,4,6,7,9,10],[ 1,  1,2,4,5,7,8,10,11],[2,6,7,9,10,12,13,15,16],[4,7,8,10,11,13,14,16,17]]
        rectangular_mesh['total_elems'] = 4
        terminal_list={}
        terminal_list['total_terminals']=1
        terminal_list['terminal_nodes']=[2]
        volume=5
        thickness=2.1
        ellipticity=1
        term_total_vol=0.04#artificial value to match with a smaller ellipsoid
        term_tissue_vol=1.77657064561
        term_tissue_diam=0.090100877305
        term_vol=placentagen.terminal_volume_to_grid(rectangular_mesh, terminal_list, noddata['nodes'],volume, thickness,ellipticity,term_total_vol,term_tissue_vol, term_tissue_diam)
        self.assertTrue(np.isclose(term_vol['term_diameter_in_grid'][0],0.08003529))

class Test_weighted_diameter(TestCase):
        
    def test_wt_diameter(self):
        term_diameter_in_grid= 0.08003529
        br_diameter_in_grid= 1.37118148e-03
        tissue_vol=0.45255089
        wt_D=placentagen.weighted_diameter_in_samp_gr(term_diameter_in_grid,br_diameter_in_grid,tissue_vol)
        self.assertTrue(np.isclose(wt_D,0.17988357))

class test_radius_br(TestCase):
        
    def test_radius_by_order(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        eldata = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        system='strahler'
        inlet_elem=0
        inlet_radius=0.1
        radius_ratio=1.53
        radius=placentagen.define_radius_by_order(noddata['nodes'], eldata['elems'], system, inlet_elem, inlet_radius, radius_ratio)
        self.assertTrue(np.isclose(radius[1],0.0653594771242))

    def test_radius_by_order_stem(self):
        noddata = placentagen.import_exnode_tree(TESTDATA_FILENAME)
        eldata = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        system='strahler'
        inlet_elem=0
        inlet_radius=0.1
        radius_ratio=1.53
        radius=placentagen.define_radius_by_order_stem(noddata['nodes'], eldata['elems'], system, TESTDDATA_FILENAME2, inlet_radius, radius_ratio)
        print(radius)
        self.assertTrue(np.isclose(radius[1],0.1))


class Test_porosity(TestCase):
      def test_porosity(self):
          vol_frac=np.array([0.3])
          porosity=placentagen.porosity(vol_frac)
          self.assertTrue(np.isclose(porosity,0.7))

class Test_node_in_sampling_grid(TestCase):
        
      def test_node_in_grid(self):
         node_loc=np.array([[0., 0.,0.,-1.05]])
         rectangular_mesh = {}
         rectangular_mesh['nodes'] =np.array( [[0.,0.,-2.],[2.,0.,-2.],[0.,2.,-2.],[2.,2.,-2.],[0.,0.,0.],[2.,0.,0.],[0.,2.,0.],[2.,2., 0.]])
         rectangular_mesh['elems']=[[0, 0, 1, 2, 3, 4, 5, 6, 7]]
         node_grid =placentagen.node_in_sampling_grid(rectangular_mesh, node_loc)
         self.assertTrue(node_grid[0][1] == 0)

#class Test_mapping_node(TestCase):
#      def test_mapping(self):
#          comp_node_elems=np.array([3])
#          non_empty_rects=np.array([2,3])
#          conductivity=np.array([0.4,0.5])
#          porosity=np.array([0.3,0.7])
#          mapping=placentagen.mapping_mesh_sampl_gr(comp_node_elems, non_empty_rects,conductivity,porosity,False,'test.txt')
#          self.assertTrue(np.isclose(mapping[0,0], 1))
#          self.assertTrue(np.isclose(mapping[0,1], 0.5)) 
#          self.assertTrue(np.isclose(mapping[0,2],0.7))

if __name__ == '__main__':
   unittest.main()   

