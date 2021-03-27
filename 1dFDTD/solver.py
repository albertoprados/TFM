import numpy as np
import copy
import math
import scipy.constants as sp


class FDTD:
    def __init__(self, mesh, pulse, set_m, set_s_m):
        self.mesh=mesh
        self.pulse=pulse
        self.set_m=set_m
        self.set_s_m=set_s_m

    def boundarymur(self, ex, boundary_low, boundary_high): 
        ex[0] = boundary_low.pop(0)
        boundary_low.append(ex[1])      

        ex[self.mesh.ncells] = boundary_high.pop(0)
        boundary_high.append(ex[self.mesh.ncells-1])

    
    def Include_SFDTD_Analysis(self,std_h,std_e,e,h):
        Stochastic_FDTD(self.mesh, self.set_m, self.set_s_m).StandardDeviation_E(e,h,std_h,std_e)
        Stochastic_FDTD(self.mesh, self.set_m, self.set_s_m).StandardDeviation_H(std_h,std_e)
    
    

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

            
            self.Include_SFDTD_Analysis(std_h,std_e,ex,hy)
            std_e_save_film[time_step][:]=std_e[:]

            t= time_step+1/2
            hy[self.pulse.k_ini] += 0.25 * self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.25 * self.pulse.pulse(t)   

       
        return ex_save_k1, ex_save_k2, ex_save_film, std_e_save_film 




class Stochastic_FDTD:
    def __init__(self, malla, set_material, set_s_material):
        self.malla=malla
        self.set_material=set_material
        self.set_s_material=set_s_material
        try:
            set_material.shape[1]
            self.num_materials=len(set_material)
        except IndexError:
            self.num_materials=1
       

        
    def StandardDeviation_H(self,std_h,std_e):
        std_h[:-1] = std_h[:-1] - 0.5 * (std_e[:-1] - std_e[1:])

    
    def StandardDeviation_E(self,e,h,std_h,std_e):   
        std_e[1:-1]=self.c1_StDe()[1:-1] * std_e[1:-1]+ \
                self.c2_StDe()[1:-1] * (std_h[1:-1] - std_h[:-2]) + \
                self.c3_StDe()[1:-1] * e[1:-1] + \
                self.c4_StDe()[1:-1] * (h[:-2] - h[1:-1])               


    def coef_aux(self):
        c_aux = np.empty(self.num_materials)

        if self.num_materials==1:
            c_aux= 2 * sp.epsilon_0 * self.set_material[0] + \
                self.malla.dt() * self.set_material[1]
        else:
            for i in range(self.num_materials):     
                c_aux[i] = 2 * sp.epsilon_0 * self.set_material[i][0] + \
                    self.malla.dt() * self.set_material[i][1]   
        
        return c_aux

    def c1_StDe(self):
        c1_coef = np.empty(self.num_materials)
        c1_malla = np.ones(self.malla.ncells+1)

        if self.num_materials==1:
            c1_coef = ( (2 * sp.epsilon_0 * self.set_material[0] - \
                self.malla.dt() * self.set_material[1]) \
                / self.coef_aux() ) \
                * math.sqrt(sp.mu_0 / sp.epsilon_0) 

            c1_malla[int(self.set_material[2]) : int(self.set_material[3])] = \
                c1_coef    
        else:
            for i in range(self.num_materials):
                c1_coef[i] = ( (2 * sp.epsilon_0 * self.set_material[i][0] - \
                    self.malla.dt() * self.set_material[i][1]) \
                    / self.coef_aux()[i] ) \
                    * math.sqrt(sp.mu_0 / sp.epsilon_0) 

                c1_malla[int(self.set_material[i][2]) : int(self.set_material[i][3])]= c1_coef[i]    
                             
    
        return c1_malla
    
    def c2_StDe(self):
        c2_coef = np.empty(self.num_materials)
        c2_malla = np.ones(self.malla.ncells+1) * 0.5 * math.sqrt(sp.mu_0 / sp.epsilon_0)

        if self.num_materials==1:
            c2_coef= 1.0 / (sp.c * self.coef_aux())

            c2_malla[int(self.set_material[2]) : int(self.set_material[3])] = \
                c2_coef
        else:
            for i in range(self.num_materials): 
                c2_coef[i]= 1.0 / (sp.c * self.coef_aux()[i]) 

                c2_malla[int(self.set_material[i][2]) : int(self.set_material[i][3])]= c2_coef[i]     
                      

        return c2_malla

    def c3_StDe(self):
        c3_coef = np.empty(self.num_materials)
        c3_malla =  np.zeros(self.malla.ncells+1)

        if self.num_materials==1:
            c3_coef = 4 * self.malla.dt()*\
            (self.set_material[1]*self.set_s_material[2]*self.set_s_material[0]-\
            self.set_material[0]*self.set_s_material[3]*self.set_s_material[1])\
            / (sp.c * np.power(self.coef_aux(),2))

            c3_malla[int(self.set_material[2]) : int(self.set_material[3])] = \
                c3_coef
        else:
            for i in range(self.num_materials):
                c3_coef[i] = 4 * self.malla.dt()*\
                (self.set_material[i][1]*self.set_s_material[i][2]*self.set_s_material[i][0]-\
                self.set_material[i][0]*self.set_s_material[i][3]*self.set_s_material[i][1])\
                / (sp.c * np.power(self.coef_aux()[i],2))    

                c3_malla[int(self.set_material[i][2]) : int(self.set_material[i][3])]= c3_coef[i]      

        return c3_malla

    def c4_StDe(self):
        c4_coef = np.empty(self.num_materials)
        c4_malla = np.zeros(self.malla.ncells+1)

        if self.num_materials==1:
            c4_coef = (1.0/(sp.c*self.coef_aux()))\
                *((2*sp.epsilon_0*self.set_s_material[0]*\
                self.set_s_material[4] + self.malla.dt()*self.set_s_material[1] * \
                self.set_s_material[5])/self.coef_aux())

            c4_malla[int(self.set_material[2]) : int(self.set_material[3])] = \
                c4_coef    
        else:
            for i in range(self.num_materials):
                c4_coef[i] = (1.0/(sp.c*self.coef_aux()[i]))\
                    *((2*sp.epsilon_0*self.set_s_material[i][0]*\
                    self.set_s_material[i][4] + self.malla.dt()*self.set_s_material[i][1] * \
                    self.set_s_material[i][5])/self.coef_aux()[i])        

                c4_malla[int(self.set_material[i][2]) : int(self.set_material[i][3])]= c4_coef[i]

        return c4_malla

               









