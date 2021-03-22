import unittest
import numpy as np
from utilities import Utilities, Source, Variables_FDTD
from mesh import Mesh
from solver import FDTD

class TestSolver(unittest.TestCase):

    def test_FFT(self):
        
        set_1_var=Variables_FDTD(200,0.001,5e-9,4,0,110,140)
        set_2_var=Variables_FDTD(200,0.001,5e-9,1,0,110,140)

        pulso=Source('gauss',40,12,20)
        
        malla1=Mesh(set_1_var)
        malla2=Mesh(set_2_var)

        #Mesh with a specific permitivity, no conductivity
        et1k1_test1, et1k2_test1= FDTD(malla1,pulso,set_1_var).FDTDLoop()
        et2k1, et2k2= FDTD(malla2,pulso,set_2_var).FDTDLoop()
        #Mesh with no permitivity, no conductivity
        et1k1_test2, et1k2_test2= FDTD(malla2,pulso,set_2_var).FDTDLoop()    


        result1,result2= Utilities().FFT(et1k1_test1,et2k1, et1k2_test1, et2k2)  
        result3,result4= Utilities().FFT(et1k1_test2,et2k1, et1k2_test2, et2k2)  
        

        for i in range(200):
            #Mesh without conductivity should return R*R+T*T=1
            self.assertLess(result1[i], 1)       
            self.assertAlmostEqual(result1[i]*result1[i]+result2[i]*result2[i],1,delta=0.1)
            #Mesh without material should return R=0, T=1
            self.assertEqual(result3[i],0)
            self.assertEqual(result4[i],1)
            






if __name__ == '__main__':
    unittest.main()