
from unittest import TestCase

import numpy as np
import unittest
import placentagen
import os
TESTDATA_FILENAME1 = os.path.join(os.path.dirname(__file__), 'Testdata/Small.exelem')

class Test_Terminal_Br(TestCase):
        
    def test_terminal_br(self):
        eldata   = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        term_br  = placentagen.terminal_branch(eldata['el_array'])
        self.assertTrue(np.isclose(term_br['terminal_el'][0,1], 2))
        
    def test_terminal_br_total(self):
        eldata   = placentagen.import_exelem_tree(TESTDATA_FILENAME1)
        term_br_total  = placentagen.terminal_branch(eldata['el_array'])
        self.assertTrue(np.isclose(term_br_total['total_terminal_el'], 2))


if __name__ == '__main__':
   unittest.main()
