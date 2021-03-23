import numpy as np
import copy
import math
import scipy.constants as sp


class FDTD:
    def __init__(self, mesh, pulse, var, svar):
        self.mesh=mesh
        self.pulse=pulse
        self.var=var
        self.svar=svar

    def boundarymur(self, ex, boundary_low, boundary_high): 
        ex[0] = boundary_low.pop(0)
        boundary_low.append(ex[1])      

        ex[self.var.ncells] = boundary_high.pop(0)
        boundary_high.append(ex[self.var.ncells-1])

    """
    def Include_SFDTD_Analysis(self,std_h,std_e,e,h):
        self.sfdtd.StandardDeviation_E(std_h,std_e,e,h)
        self.sfdtd.StandardDeviation_H(std_h,std_e)
    """
    

    def FDTDLoop(self):
        
        dt=self.var.dt()
        nsteps= int(self.var.time / dt)

        #Def
        ex=np.zeros(self.var.ncells+1)
        hy=np.zeros(self.var.ncells+1)

        std_e=np.zeros(self.var.ncells+1)
        std_h=np.zeros(self.var.ncells+1)
        
        ex_save_k1=np.empty(nsteps+1)
        ex_save_k2=np.empty(nsteps+1)

        #ex_save_film=np.empty((nsteps+1,self.var.ncells+1))
        
        ca=self.mesh.material()[0][1:-1]
        cb=self.mesh.material()[1][1:-1]

        boundary_low = [0, 0]
        boundary_high = [0, 0]
       
        for time_step in range(1, nsteps + 1):

            ex[1:-1] = ca * ex[1:-1] + cb * (hy[:-2] - hy[1:-1])
            
            #Guardo los valores a representar
            #ex_save_film[time_step][:]=ex[:]
            
            #Guardo los valores para calcular la transformada
            ex_save_k1[time_step]=ex[self.var.FFTpoints()[0]]
            ex_save_k2[time_step]=ex[self.var.FFTpoints()[1]]
           
            ex[self.pulse.k_ini] +=  0.5* self.pulse.pulse(time_step) 
            
            self.boundarymur(ex,boundary_low,boundary_high)  
            
            
            hy[:-1] = hy[:-1] + 0.5 * (ex[:-1] - ex[1:])   

            Stochastic_FDTD(self.var,self.svar).StandardDeviation_E(std_h,std_e,ex,hy)
            Stochastic_FDTD(self.var,self.svar).StandardDeviation_H(std_h,std_e)

            t= time_step+1/2
            hy[self.pulse.k_ini] += 0.25 * self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.25 * self.pulse.pulse(t)   

            


        return ex_save_k1, ex_save_k2, std_h, std_e #ex_save_film




class Stochastic_FDTD:
    def __init__(self, var1, var2):
        self.var1=var1
        self.var2=var2

    def StandardDeviation_H(self,std_h,std_e):
        std_h[:-1] = std_h[:-1] - 0.5 * (std_e[:-1] - std_e[1:])

    
    def StandardDeviation_E(self,std_h,std_e,e,h):   
        std_e[1:-1]=self.c1_StDe() * std_e[1:-1] + \
                    self.c2_StDe() * (std_h[1:-1] - std_h[:-2]) + \
                    self.c3_StDe() * e[1:-1] + \
                    self.c4_StDe() * (h[:-2] - h[1:-1])     


    def coef_aux(self):
        return 2 * sp.epsilon_0 * self.var1.epsilon_r + \
               self.var1.dt() * self.var1.sigma

    def c1_StDe(self):
        c1=( (2 * sp.epsilon_0 * self.var1.epsilon_r - \
            self.var1.dt() * self.var1.sigma) \
            / self.coef_aux()) * math.sqrt(sp.mu_0 / sp.epsilon_0)  

        return c1
    
    def c2_StDe(self):
        return 1.0 / (sp.c * self.coef_aux())

    def c3_StDe(self):
        c3=4*self.var1.dt()*\
           (self.var1.sigma*self.var2.c_eps_E*self.var2.std_eps-\
            self.var1.epsilon_r*self.var2.c_sigma_E*self.var2.std_sigma)\
            / (sp.c * math.pow(self.coef_aux(),2))

        return c3

    def c4_StDe(self):
        c4=(1.0/(sp.c*self.coef_aux()))*((2*sp.epsilon_0*self.var2.std_eps*\
            self.var2.c_eps_H + self.var1.dt()*self.var2.std_sigma * \
            self.var2.c_sigma_H)/self.coef_aux())

        return c4

               









