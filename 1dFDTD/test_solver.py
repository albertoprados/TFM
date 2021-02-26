import unittest
import numpy as np
from mesh import Mesh
from solver import FDTD, Utilities, Source

class TestSolver(unittest.TestCase):

    def test_FFT(self):
        
        pulso=Source('gauss',40,12,20)

        et1k1, et1k2= FDTD(Mesh(200,0.001,4,0.04,110,140), \
                      pulso,5e-9).FDTDLoop(40,160)
        et2k1, et2k2= FDTD(Mesh(200,0.001,1,0,110,140), \
                      pulso,5e-9).FDTDLoop(40,160)

        result=Utilities().FFT(et1k1,et2k1, et1k2, et2k2)[0]    

        for i in range(len(result)):
            self.assertLess(result[i], 1)       







if __name__ == '__main__':
    unittest.main()