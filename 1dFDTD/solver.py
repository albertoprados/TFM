import numpy as np
import copy
import math
import scipy.constants as sp


class FDTD:
    def __init__(self, mesh, pulse):
        self.mesh=mesh
        self.pulse=pulse

    def boundarymur(self, ex, boundary_low, boundary_high): 
        ex[0] = boundary_low.pop(0)
        boundary_low.append(ex[1])      

        ex[self.mesh.ncells] = boundary_high.pop(0)
        boundary_high.append(ex[self.mesh.ncells-1])

    """
    def Include_SFDTD_Analysis(self,std_h,std_e,e,h):
        self.sfdtd.StandardDeviation_E(std_h,std_e,e,h)
        self.sfdtd.StandardDeviation_H(std_h,std_e)
    """
    

    def FDTDLoop(self,time):
        
        dt=self.mesh.dt()
        nsteps= int(time / dt)

        #Def
        ex=np.zeros(self.mesh.ncells+1)
        hy=np.zeros(self.mesh.ncells+1)
        
        ex_save_k1=np.empty(nsteps+1)
        ex_save_k2=np.empty(nsteps+1)
        
        #std_h=np.zeros(self.var.ncells+1)
        #std_e=np.zeros(self.var.ncells+1)

        ex_save_film=np.empty((nsteps+1,self.mesh.ncells+1))
        #std_e_save_film=np.empty((nsteps+1,self.var.ncells+1))

        ca=self.mesh.material()[0]
        cb=self.mesh.material()[1]

        boundary_low = [0, 0]
        boundary_high = [0, 0]
       
        for time_step in range(1, nsteps + 1):

            ex[1:-1] = ca[1:-1] * ex[1:-1] + cb[1:-1] * (hy[:-2] - hy[1:-1])
            
            #Guardo los valores a representar
            ex_save_film[time_step][:]=ex[:]
            
            #Guardo los valores para calcular la transformada
            ex_save_k1[time_step]=ex[self.mesh.FFTpoints()[0]]
            ex_save_k2[time_step]=ex[self.mesh.FFTpoints()[1]]
           
            ex[self.pulse.k_ini] +=  0.5* self.pulse.pulse(time_step) 
            
            self.boundarymur(ex,boundary_low,boundary_high)  
            
            
            hy[:-1] = hy[:-1] + 0.5 * (ex[:-1] - ex[1:])   

            
            #Stochastic_FDTD(self.var,self.svar).StandardDeviation_E(ex,hy,std_h,std_e)
            #Stochastic_FDTD(self.var,self.svar).StandardDeviation_H(std_h,std_e)
            #std_e_save_film[time_step][:]=std_e[:]

            t= time_step+1/2
            hy[self.pulse.k_ini] += 0.25 * self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.25 * self.pulse.pulse(t)   

       
        return ex_save_k1, ex_save_k2, ex_save_film #std_e, std_h, std_e_save_film 


"""

class Stochastic_FDTD:
    def __init__(self, set_par_1, set_par_2):
        self.var1=var1
        self.var2=var2
       

        
    def StandardDeviation_H(self,std_h,std_e):
        std_h[:-1] = std_h[:-1] - 0.5 * (std_e[:-1] - std_e[1:])

    
    def StandardDeviation_E(self,parameters,e,h,std_h,std_e):   
        std_e[1:-1]=self.c1_StDe(parameters)[1:-1] * std_e[1:-1]+ \
                self.c2_StDe(parameters)[1:-1] * (std_h[1:-1] - std_h[:-2]) + \
                self.c3_StDe(parameters)[1:-1] * e[1:-1] + \
                self.c4_StDe(parameters)[1:-1] * (h[:-2] - h[1:-1])               


    def coef_aux(self,parameters):
        c_aux= 2 * sp.epsilon_0 * parameters.epsilon_r + \
            parameters.dt() * parameters.sigma
        
        return c_aux

    def c1_StDe(self,parameters):
        c1=( (2 * sp.epsilon_0 * parameters.epsilon_r - \
            parameters.dt() * parameters.sigma) \
            / self.coef_aux(parameters) ) \
            * math.sqrt(sp.mu_0 / sp.epsilon_0)  
    
        return c1
    
    def c2_StDe(self,parameters):
        c2= 1.0 / (sp.c * self.coef_aux(parameters))
    
        return c2

    def c3_StDe(self,parameters):
        c3=4*parameters.dt()*\
           (parameters.sigma*parameters.c_eps_E*parameters.std_eps-\
            parameters.epsilon_r*parameters.c_sigma_E*parameters.std_sigma)\
            / (sp.c * np.power(self.coef_aux(parameters),2))

        return c3

    def c4_StDe(self,parameters):
        c4= (1.0/(sp.c*self.coef_aux(parameters)))\
            *((2*sp.epsilon_0*parameters.std_eps*\
            parameters.c_eps_H + parameters.dt()*parameters.std_sigma * \
            parameters.c_sigma_H)/self.coef_aux(parameters))

        return c4

               

"""







