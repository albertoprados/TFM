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

    
    def Include_SFDTD_Analysis(self,std_h,std_e,e,h,bound_low_s,bound_high_s):
        Stochastic_FDTD(self.mesh).StandardDeviation_E(e,h,std_h,std_e)
        Stochastic_FDTD(self.mesh).BoundaryCondition(std_e,bound_low_s,bound_high_s)
        
        Stochastic_FDTD(self.mesh).StandardDeviation_H(std_h,std_e)
        Stochastic_FDTD(self.mesh).BoundaryCondition2(std_h)
    
    

    def FDTDLoop(self,time):
        
        dt=self.mesh.dt()
        nsteps= int(time / dt)

        #Def
        ex=np.zeros(self.mesh.ncells+1)
        hy=np.zeros(self.mesh.ncells+1)
        
        ex_save_k1=np.empty(nsteps+1)
        ex_save_k2=np.empty(nsteps+1)
        
        std_h=np.zeros(self.mesh.ncells+1)
        std_e=np.zeros(self.mesh.ncells+1)

        #Saving values film
        ex_save_film=np.empty((nsteps+1,self.mesh.ncells+1))
        std_e_save_film=np.empty((nsteps+1,self.mesh.ncells+1))

        ca=self.mesh.materials()[0]
        cb=self.mesh.materials()[1]

        boundary_low = [0, 0]
        boundary_high = [0, 0]
        bound_low_s = [0, 0]
        bound_high_s = [0, 0]
        
       
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


            self.Include_SFDTD_Analysis(std_h,std_e,ex,hy,bound_low_s,bound_high_s)
            std_e_save_film[time_step][:]=std_e[:]
          
            t= time_step+1/2
            hy[self.pulse.k_ini] += 0.25 * self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.25 * self.pulse.pulse(t)   

       
        return ex_save_k1, ex_save_k2, ex_save_film, std_e_save_film 




class Stochastic_FDTD:
    def __init__(self, malla):
        self.malla=malla
        
               
    def StandardDeviation_H(self,std_h,std_e):
        std_h[:-1] = std_h[:-1] - 0.5 * (std_e[:-1] - std_e[1:])

    
    def StandardDeviation_E(self,e,h,std_h,std_e):   
        std_e[1:-1]=self.malla.c1_StDe()[1:-1] * std_e[1:-1]+ \
                self.malla.c2_StDe()[1:-1] * (std_h[1:-1] - std_h[:-2]) + \
                self.malla.c3_StDe()[1:-1] * e[1:-1] + \
                self.malla.c4_StDe()[1:-1] * (h[:-2] - h[1:-1])               

    def BoundaryCondition(self,std_e,bound_low_s,bound_high_s):
        std_e[0] = bound_low_s.pop(0)
        bound_low_s.append(std_e[1])      

        std_e[self.malla.ncells] = bound_high_s.pop(0)
        bound_high_s.append(std_e[self.malla.ncells-1])        
    
    def  BoundaryCondition2(self, std_h):
        std_h[79]=std_h[80]
        std_h[109]=std_h[110]
        std_h[139]=std_h[140]
        








