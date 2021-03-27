import numpy as np
import scipy.constants as sp
from math import pi, sin, exp


class Mesh:
    def __init__(self, ncells, ddx, parameters):    
        self.ncells=ncells
        self.ddx=ddx
        self.parameters=parameters

    def dt(self):
        return self.ddx/(2*sp.c)    
    
    def FFTpoints(self):
        return self.ncells - 120, self.ncells - 40


    def materials(self):
        try:
            self.parameters.shape[1]
            num_materials=len(self.parameters)
        except IndexError:
            num_materials=1

        eaf = np.empty(num_materials)
        ca = np.ones(self.ncells+1)
        cb = np.ones(self.ncells+1) * 0.5
       
        if num_materials==1:
            eaf = self.dt() * self.parameters[1] \
                /(2 * sp.epsilon_0 * self.parameters[0]) 
            ca[int(self.parameters[2]) : int(self.parameters[3])] = \
                (1 - eaf ) / (1 + eaf )
            cb[int(self.parameters[2]) : int(self.parameters[3])] = \
                0.5 / (self.parameters[0] * (1 + eaf ))        
        else:
            for i in range(num_materials):     
                eaf[i] = self.dt() * self.parameters[i][1] \
                    /(2 * sp.epsilon_0 * self.parameters[i][0]) 
                ca[int(self.parameters[i][2]) : int(self.parameters[i][3])] = \
                    (1 - eaf[i] ) / (1 + eaf[i] )
                cb[int(self.parameters[i][2]) : int(self.parameters[i][3])] = \
                    0.5 / (self.parameters[i][0] * (1 + eaf[i] ))
       
        
        return  ca, cb


class Materials:

    def material(self,epsilon_r,sigma,start_m,end_m):
        parameters=list(locals().values())

        #Avoid getting self parameter 
        material=np.empty(len(parameters)-1)

        for i in range(1,len(parameters)):
            material[i-1] = parameters[i]

        return material

    def material_matrix(self, materials):
        num_materials=len(materials)

        material_matrix=np.empty((num_materials,len(materials[0])))

        for i in range(num_materials):
            material_matrix[i][:]=materials[i]
        
        return material_matrix

    def s_material(self,std_eps,std_sigma,c_eps_E,c_sigma_E,c_eps_H,c_sigma_H):    
        parameters=list(locals().values())

        #Avoid getting self parameter 
        stochastic_material=np.empty(len(parameters)-1)

        for i in range(1,len(parameters)):
            stochastic_material[i-1] = parameters[i]

        return stochastic_material

    def s_material_matrix(self, s_materials):
        num_materials=len(s_materials)

        stochastic_material_matrix=np.empty((num_materials,len(s_materials[0])))

        for i in range(num_materials):
            stochastic_material_matrix[i][:]=s_materials[i]
        
        return stochastic_material_matrix    