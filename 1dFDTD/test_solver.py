import unittest
import numpy as np
from utilities import Utilities, Source
from mesh import Mesh, Materials
from solver import FDTD

class TestSolver(unittest.TestCase):

    def test_FFT(self):
        
        #Permitivity,Conductivity,Start_Point,End_Point
        material_1=[4,0,110,140]
        material_2=[3,0,80,110]
        set_m=[material_1,material_2]
        #Std_Permitivity,Std_Conductivity,Corr_Eps_E,Corr_Sigma_E,Corr_Eps_H,Corr_Sigma_H
        s_material_1=[1,0.002,1,1,1,1]
        s_material_2=[0.75,0.003,1,1,1,1]
        set_s_m=[s_material_1,s_material_2]

        parameters=Materials(set_m,set_s_m)

        air=[[1,0,0,200]]
        s_air=[[0,0,0,0,0,0]]
        parameters_air=Materials(air,s_air)

        malla1=Mesh(200,0.001,parameters)
        malla2=Mesh(200,0.001,parameters_air)

        pulso=Source('gauss',40,12,20)
        
        #Mesh with a specific permitivity, no conductivity
        et1k1_test1= FDTD(malla1,pulso).FDTDLoop(5e-9)[0]
        et1k2_test1= FDTD(malla1,pulso).FDTDLoop(5e-9)[1]
        et2k1= FDTD(malla2,pulso).FDTDLoop(5e-9)[0]
        et2k2= FDTD(malla2,pulso).FDTDLoop(5e-9)[1]
        #Mesh with no permitivity, no conductivity
        et1k1_test2=FDTD(malla2,pulso).FDTDLoop(5e-9)[0]
        et1k2_test2= FDTD(malla2,pulso).FDTDLoop(5e-9)[1]    


        result1,result2= Utilities().FFT(et1k1_test1,et2k1, et1k2_test1, et2k2)  
        result3,result4= Utilities().FFT(et1k1_test2,et2k1, et1k2_test2, et2k2)  
        

        for i in range(100):
            #Mesh without conductivity should return R*R+T*T=1
            self.assertLess(result1[i], 1)       
            self.assertAlmostEqual(result1[i]*result1[i]+result2[i]*result2[i],1,delta=0.1)
            #Mesh without material should return R=0, T=1
            self.assertEqual(result3[i],0)
            self.assertEqual(result4[i],1)
            
    def test_coef_SFDTD(self):
        
        #Permitivity,Conductivity,Start_Point,End_Point
        material_1=[4,0.04,110,140]
        set_m=[material_1]
        s_material_1=[1,0.002,1,1,1,1]
        set_s_m=[s_material_1]

        parameters=Materials(set_m,set_s_m)

        air=[[1,0,0,200]]
        s_air=[[0,0,0,0,0,0]]
        parameters_air=Materials(air,s_air)

        malla1=Mesh(200,0.001,parameters)
        malla2=Mesh(200,0.001,parameters_air)
  
        coef_1_m=malla1.c1_StDe()
        
        

        for i in range(110,140):
            #Coef_1 should return a value between 0 & 1 in material
            self.assertLess(coef_1_m[i], 1)       
            





if __name__ == '__main__':
    unittest.main()