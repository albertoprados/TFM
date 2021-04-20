import numpy as np
import copy
import math
import scipy.constants as sp


class FDTD:
    def __init__(self, mesh, pulse):
        self.mesh=mesh
        self.pulse=pulse


    def boundarymur(self, ex, ex_old):
        ncells, dt, ddx= self.mesh.ncells, self.mesh.dt(), self.mesh.ddx

        c_bound=(sp.c*dt-ddx)/(sp.c*dt+ddx)

        ex[0]=ex_old[1] + c_bound * (ex[1]-ex_old[0])
        ex[ncells]=ex_old[ncells-1] + c_bound * (ex[ncells-1]-ex_old[ncells])


    def Include_SFDTD_Analysis(self,std_h,std_e,e,h,std_e_old):
        Stochastic_FDTD(self.mesh).StandardDeviation_E(std_h,std_e,e,h)
        Stochastic_FDTD(self.mesh).BoundaryCondition(std_e,std_e_old)

        Stochastic_FDTD(self.mesh).StandardDeviation_H(std_h,std_e)
        
        
        
    def FDTDLoop(self,time):
        dt=self.mesh.dt()
        nsteps= int(time / dt)

        ex=np.zeros(self.mesh.ncells+1)
        hy=np.zeros(self.mesh.ncells)
        ex_old=np.zeros(self.mesh.ncells+1)

        #Fourier transform
        ex_save_k1=np.empty(nsteps+1)
        ex_save_k2=np.empty(nsteps+1)
        
        #Standard Deviation
        std_e=np.zeros(self.mesh.ncells+1)
        std_h=np.zeros(self.mesh.ncells)
        std_e_old=std_e=np.zeros(self.mesh.ncells+1)

        #Saving values for film
        ex_save_film=np.empty((nsteps+1,self.mesh.ncells+1))
        std_e_save_film=np.empty((nsteps+1,self.mesh.ncells+1))

        ca=self.mesh.materials()[0]
        cb=self.mesh.materials()[1]
        cc=self.mesh.materials()[2]

       
        for time_step in range(1, nsteps + 1):
            ex_old=copy.deepcopy(ex)
            
            ex[1:-1] = ca[1:-1] * ex[1:-1] + cb[1:-1] * (hy[:-1] - hy[1:])
            
            #Guardo los valores a representar
            ex_save_film[time_step][:]=ex[:]
            
            #Guardo los valores para calcular la transformada
            ex_save_k1[time_step]=ex[self.mesh.FFTpoints()[0]]
            ex_save_k2[time_step]=ex[self.mesh.FFTpoints()[1]]
           
            ex[self.pulse.k_ini] += 0.5 * self.pulse.pulse(time_step) 
            
            self.boundarymur(ex,ex_old)  
            
            hy[:] = hy[:] + cc * (ex[:-1] - ex[1:])   


            std_e_old=copy.deepcopy(std_e)          
            self.Include_SFDTD_Analysis(std_h,std_e,ex_old,hy,std_e_old)
            std_e_save_film[time_step][:]=std_e[:]
            
            
            t= time_step + 1/2
            hy[self.pulse.k_ini] += 0.25 * self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.25 * self.pulse.pulse(t)   
            
       
        return ex_save_k1, ex_save_k2, ex_save_film, std_e_save_film * std_e_save_film




class Stochastic_FDTD:
    def __init__(self, mesh):
        self.mesh=mesh
        
               
    def StandardDeviation_H(self,std_h,std_e):
        std_h[:] = std_h[:] - 0.5 * (std_e[:-1] - std_e[1:])

    
    def StandardDeviation_E(self,std_h,std_e,e,h):   
        std_e[1:-1]=self.mesh.c1_StDe()[1:-1] * std_e[1:-1]+ \
                self.mesh.c2_StDe()[1:-1] * (std_h[1:] - std_h[:-1]) + \
                self.mesh.c3_StDe()[1:-1] * e[1:-1] + \
                self.mesh.c4_StDe()[1:-1] * (h[:-1] - h[1:])               


    def BoundaryCondition(self,std_e,std_e_old):
        ncells, dt, ddx= self.mesh.ncells, self.mesh.dt(), self.mesh.ddx

        c_bound=(sp.c*dt-ddx)/(sp.c*dt+ddx)

        std_e[0]=std_e_old[1] + c_bound * (std_e[1]-std_e_old[0])
        std_e[ncells]=std_e_old[ncells-1] + c_bound * \
             (std_e[ncells-1]-std_e_old[ncells])    
    
   








