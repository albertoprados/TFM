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
        return self.ncells - 260, self.ncells - 20


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

    

    def coef_aux(self):
        c_aux = np.empty(self.par.num_materials)

        for i in range(self.par.num_materials):     
            c_aux[i] = 2 * sp.epsilon_0 * self.par.epsilon_r()[i] + \
                self.dt() * self.par.sigma()[i]   
        
        return c_aux


    def c1_StDe(self):
        c1_coef = np.empty(self.par.num_materials)
        c1_malla = np.ones(self.ncells+1)

        for i in range(self.par.num_materials):
            c1_coef[i] = ( (2 * sp.epsilon_0 * self.par.epsilon_r()[i] - \
                self.dt() * self.par.sigma()[i]) \
                / self.coef_aux()[i] ) 
                

            c1_malla[self.par.start_m()[i] : self.par.end_m()[i]]= c1_coef[i]    
                             
        return c1_malla
    

    def c2_StDe(self):
        c2_coef = np.empty(self.par.num_materials)
        c2_malla = np.ones(self.ncells+1) * (self.dt()/(self.ddx*sp.epsilon_0)) \
                * math.sqrt(sp.epsilon_0/sp.mu_0)

        for i in range(self.par.num_materials): 
            c2_coef[i]= (2.0 * self.dt() / (self.ddx * self.coef_aux()[i])) \
                        * math.sqrt(sp.epsilon_0/sp.mu_0)

            c2_malla[self.par.start_m()[i] : self.par.end_m()[i]]= c2_coef[i]     
                      
        return c2_malla


    def c3_StDe(self):
        c3_coef = np.empty(self.par.num_materials)
        c3_malla =  np.zeros(self.ncells+1)

        for i in range(self.par.num_materials):
            c3_coef[i] = 4 * self.dt() * sp.epsilon_0 * \
            (self.par.sigma()[i]*self.par.c_eps_E()[i]*self.par.std_eps_r()[i]-\
            self.par.epsilon_r()[i]*self.par.c_sigma_E()[i]*self.par.std_sigma()[i])\
            / (np.power(self.coef_aux()[i],2))    

            c3_malla[self.par.start_m()[i] : self.par.end_m()[i]]= c3_coef[i]      

        return c3_malla


    def c4_StDe(self):
        c4_coef = np.empty(self.par.num_materials)
        c4_malla = np.zeros(self.ncells+1)

        for i in range(self.par.num_materials):
            c4_coef[i] = ((2.0 * self.dt()/(self.ddx*self.coef_aux()[i]))\
                *((2*sp.epsilon_0*self.par.std_eps_r()[i]*\
                self.par.c_eps_H()[i] + self.dt()*self.par.std_sigma()[i] * \
                self.par.c_sigma_H()[i])/self.coef_aux()[i])) * \
                math.sqrt(sp.epsilon_0/sp.mu_0)           

            c4_malla[self.par.start_m()[i] : self.par.end_m()[i]]= c4_coef[i]

        return c4_malla





class Materials:
    def __init__(self, materiales, s_materiales):
        self.materiales=materiales
        self.s_materiales=s_materiales
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



