import unittest
import numpy as np
from mesh import Mesh
from solver import FDTD, Utilities, Source

class TestSolver(unittest.TestCase):

    def test_FFT(self):
        malla1=Mesh(200,0.001,4,0.04,110,140)
        malla2=Mesh(200,0.001,1,0,110,140)
        pulso=Source('gauss',40,12,20)

        et1k1= FDTD(malla1,pulso,5e-9).FDTDLoop(40,160)[0]
        e2tk1= FDTD(malla2,pulso,5e-9).FDTDLoop(40,160)[0]
        et1k2= FDTD(malla1,pulso,5e-9).FDTDLoop(40,160)[1]
        e2tk2= FDTD(malla2,pulso,5e-9).FDTDLoop(40,160)[1]

        result=Utilities().FFT(et1k1,e2tk1, et1k2, e2tk2)[0]    

        for i in range(len(result)):
            self.assertLess(result[i], 2.6)       







if __name__ == '__main__':
    unittest.main()