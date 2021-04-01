import numpy as np
import scipy.constants as sp
import math


class Mesh:
    def __init__(self, ncells, ddx, parameters, s_parameters):    
        self.ncells=ncells
        self.ddx=ddx
        self.parameters=parameters
        self.s_parameters=s_parameters
        self.num_materials=len(parameters)
        
    def dt(self):
        return self.ddx/(2*sp.c)    
    
    def FFTpoints(self):
        return self.ncells - 120, self.ncells - 40


    def materials(self):
        

        eaf = np.empty(self.num_materials)
        ca = np.ones(self.ncells+1)
        cb = np.ones(self.ncells+1) * 0.5
       
        for i in range(self.num_materials):     
            eaf[i] = self.dt() * self.parameters[i][1] \
                    /(2 * sp.epsilon_0 * self.parameters[i][0]) 
            ca[int(self.parameters[i][2]) : int(self.parameters[i][3])] = \
                    (1 - eaf[i] ) / (1 + eaf[i] )
            cb[int(self.parameters[i][2]) : int(self.parameters[i][3])] = \
                    0.5 / (self.parameters[i][0] * (1 + eaf[i] ))
       
        
        return  ca, cb



    def coef_aux(self):
        c_aux = np.empty(self.num_materials)

        for i in range(self.num_materials):     
            c_aux[i] = 2 * sp.epsilon_0 * self.parameters[i][0] + \
                self.dt() * self.parameters[i][1]   
        
        return c_aux

    def c1_StDe(self):
        c1_coef = np.empty(self.num_materials)
        c1_malla = np.ones(self.ncells+1)


        for i in range(self.num_materials):
            c1_coef[i] = ( (2 * sp.epsilon_0 * self.parameters[i][0] - \
                self.dt() * self.parameters[i][1]) \
                / self.coef_aux()[i] ) \
                * math.sqrt(sp.mu_0 / sp.epsilon_0) 

            c1_malla[int(self.parameters[i][2]) : int(self.parameters[i][3])]= c1_coef[i]    
                             
        return c1_malla
    
    def c2_StDe(self):
        c2_coef = np.empty(self.num_materials)
        c2_malla = np.ones(self.ncells+1) * 0.5 * math.sqrt(sp.mu_0 / sp.epsilon_0)

        for i in range(self.num_materials): 
            c2_coef[i]= 1.0 / (sp.c * self.coef_aux()[i]) 

            c2_malla[int(self.parameters[i][2]) : int(self.parameters[i][3])]= c2_coef[i]     
                      
        return c2_malla

    def c3_StDe(self):
        c3_coef = np.empty(self.num_materials)
        c3_malla =  np.zeros(self.ncells+1)

        for i in range(self.num_materials):
            c3_coef[i] = 4 * self.dt()*\
            (self.parameters[i][1]*self.s_parameters[i][2]*self.s_parameters[i][0]-\
            self.parameters[i][0]*self.s_parameters[i][3]*self.s_parameters[i][1])\
            / (sp.c * np.power(self.coef_aux()[i],2))    

            c3_malla[int(self.parameters[i][2]) : int(self.parameters[i][3])]= c3_coef[i]      

        return c3_malla

    def c4_StDe(self):
        c4_coef = np.empty(self.num_materials)
        c4_malla = np.zeros(self.ncells+1)

        for i in range(self.num_materials):
            c4_coef[i] = (1.0/(sp.c*self.coef_aux()[i]))\
                *((2*sp.epsilon_0*self.s_parameters[i][0]*\
                self.s_parameters[i][4] + self.dt()*self.s_parameters[i][1] * \
                self.s_parameters[i][5])/self.coef_aux()[i])        

            c4_malla[int(self.parameters[i][2]) : int(self.parameters[i][3])]= c4_coef[i]

        return c4_malla



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