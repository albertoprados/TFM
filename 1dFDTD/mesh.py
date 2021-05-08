import numpy as np
import scipy.constants as sp
import math


class Mesh:
    def __init__(self, ncells, ddx, par):    
        self.ncells=ncells
        self.ddx=ddx
        self.par=par
        
    def dt(self):
        return self.ddx/(2*sp.c)    
    
    def FFTpoints(self):
        return 250, 700


    def materials(self):
        eaf = np.empty(self.par.num_materials)
        #Coef. for update E equation
        ca = np.ones(self.ncells+1)
        c_aux=(self.dt()/(self.ddx*sp.epsilon_0))*math.sqrt(sp.epsilon_0/sp.mu_0)
        cb = np.ones(self.ncells+1) * c_aux

        #Coef. for update H equation
        cc=(self.dt()/(sp.mu_0*self.ddx))*math.sqrt(sp.mu_0/sp.epsilon_0)

        for i in range(self.par.num_materials):     
            eaf[i] = self.dt() * self.par.sigma()[i] \
                    /(2 * sp.epsilon_0 * self.par.epsilon_r()[i]) 
            ca[self.par.start_m()[i] : self.par.end_m()[i]] = \
                    (1 - eaf[i] ) / (1 + eaf[i] )
            cb[self.par.start_m()[i] : self.par.end_m()[i]] = \
                    c_aux / (self.par.epsilon_r()[i] * (1 + eaf[i] ))
       
        return  ca, cb, cc

    


class Materials:
    def __init__(self, materiales):
        self.materiales=materiales
        self.num_materials=len(materiales)

    def epsilon_r(self):
        eps=np.empty(self.num_materials)

        for i in range(self.num_materials):
            eps[i]=self.materiales[i][0]

        return eps

    def sigma(self):
        sigma=np.empty(self.num_materials)
         
        for i in range(self.num_materials):
            sigma[i]=self.materiales[i][1]

        return sigma    

    def start_m(self):
        start_m=np.empty(self.num_materials,dtype=int)

        for i in range(self.num_materials):
            start_m[i]=self.materiales[i][2]

        return start_m    

    def end_m(self):
        end_m=np.empty(self.num_materials,dtype=int)

        for i in range(self.num_materials):
            end_m[i]=self.materiales[i][3]

        return end_m   


class S_Materials:         
    def __init__(self, s_materiales):
        self.s_materiales=s_materiales
        self.num_materials=len(s_materiales)

    def std_eps_r(self):
        std_eps_r=np.empty(self.num_materials)

        for i in range(self.num_materials):
            std_eps_r[i]=self.s_materiales[i][0]

        return std_eps_r   

    def std_sigma(self):
        std_sigma=np.empty(self.num_materials)

        for i in range(self.num_materials):
            std_sigma[i]=self.s_materiales[i][1]

        return std_sigma 
          
    def c_eps_E(self):
        c_eps_E=np.empty(self.num_materials)

        for i in range(self.num_materials):
            c_eps_E[i]=self.s_materiales[i][2]

        return c_eps_E
  
    def c_sigma_E(self):
        c_sigma_E=np.empty(self.num_materials)

        for i in range(self.num_materials):
            c_sigma_E[i]=self.s_materiales[i][3]

        return c_sigma_E

    def c_eps_H(self):
        c_eps_H=np.empty(self.num_materials)

        for i in range(self.num_materials):
            c_eps_H[i]=self.s_materiales[i][4]

        return c_eps_H
        
    def c_sigma_H(self):
        c_sigma_H=np.empty(self.num_materials)

        for i in range(self.num_materials):
            c_sigma_H[i]=self.s_materiales[i][5]

        return c_sigma_H



