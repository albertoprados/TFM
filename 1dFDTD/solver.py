from mesh import Mesh
import numpy as np
import copy
import math
import scipy.constants as sp


class FDTD:
    def __init__(self, mesh, par, s_par, pulse):
        self.mesh=mesh
        self.pulse=pulse
        self.par=par
        self.s_par=s_par

    def boundarymur(self, ex, ex_old):
        ncells, dt, ddx= self.mesh.ncells, self.mesh.dt(), self.mesh.ddx

        c_bound=(sp.c*dt-ddx)/(sp.c*dt+ddx)

        ex[0]=ex_old[1] + c_bound * (ex[1]-ex_old[0])
        ex[ncells]=ex_old[ncells-1] + c_bound * (ex[ncells-1]-ex_old[ncells])


    def Include_SFDTD_Analysis(self,std_h,std_e,e,h,std_e_old):
        Stochastic_FDTD(self.mesh,self.par,self.s_par).StandardDeviation_E(std_h,std_e,e,h)
        Stochastic_FDTD(self.mesh,self.par,self.s_par).BoundaryCondition(std_e,std_e_old)

        Stochastic_FDTD(self.mesh,self.par,self.s_par).StandardDeviation_H(std_h,std_e)
        
        
        
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

        ca, cb, cc = self.mesh.materials()
    

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
    def __init__(self, mesh, par, s_par):
        self.mesh=mesh
        self.par=par
        self.s_par=s_par
               
    def StandardDeviation_H(self,std_h,std_e):
        std_h[:] = std_h[:] - 0.5 * (std_e[:-1] - std_e[1:])

    
    def StandardDeviation_E(self,std_h,std_e,e,h):   
        std_e[1:-1]=self.c1_StDe()[1:-1] * std_e[1:-1]+ \
                self.c2_StDe()[1:-1] * (std_h[1:] - std_h[:-1]) + \
                self.c3_StDe()[1:-1] * e[1:-1] + \
                self.c4_StDe()[1:-1] * (h[:-1] - h[1:])               


    def BoundaryCondition(self,std_e,std_e_old):
        ncells, dt, ddx= self.mesh.ncells, self.mesh.dt(), self.mesh.ddx

        c_bound=(sp.c*dt-ddx)/(sp.c*dt+ddx)

        std_e[0]=std_e_old[1] + c_bound * (std_e[1]-std_e_old[0])
        std_e[ncells]=std_e_old[ncells-1] + c_bound * \
             (std_e[ncells-1]-std_e_old[ncells])    



    def coef_aux(self):
        c_aux = np.empty(self.par.num_materials)

        for i in range(self.par.num_materials):     
            c_aux[i] = 2 * sp.epsilon_0 * self.par.epsilon_r()[i] + \
                self.mesh.dt() * self.par.sigma()[i]   
        
        return c_aux

    def c1_StDe(self):
        c1_coef = np.empty(self.par.num_materials)
        c1_malla = np.ones(self.mesh.ncells+1)

        for i in range(self.par.num_materials):
            c1_coef[i] = ( (2 * sp.epsilon_0 * self.par.epsilon_r()[i] - \
                self.mesh.dt() * self.par.sigma()[i]) \
                / self.coef_aux()[i] ) 
                

            c1_malla[self.par.start_m()[i] : self.par.end_m()[i]]= c1_coef[i]    
                             
        return c1_malla
    
    def c2_StDe(self):
        c2_coef = np.empty(self.par.num_materials)
        c2_malla = np.ones(self.mesh.ncells+1) * (self.mesh.dt()/ \
            (self.mesh.ddx*sp.epsilon_0)) * math.sqrt(sp.epsilon_0/sp.mu_0)
                
        for i in range(self.par.num_materials): 
            c2_coef[i]= (2.0 * self.mesh.dt() / (self.mesh.ddx * self.coef_aux()[i])) \
                        * math.sqrt(sp.epsilon_0/sp.mu_0)

            c2_malla[self.par.start_m()[i] : self.par.end_m()[i]]= c2_coef[i]     
                      
        return c2_malla

    def c3_StDe(self):
        c3_coef = np.empty(self.par.num_materials)
        c3_malla =  np.zeros(self.mesh.ncells+1)

        for i in range(self.par.num_materials):
            c3_coef[i] = 4 * self.mesh.dt() * sp.epsilon_0 * \
            (self.par.sigma()[i]*self.s_par.c_eps_E()[i]*self.s_par.std_eps_r()[i]-\
            self.par.epsilon_r()[i]*self.s_par.c_sigma_E()[i]*self.s_par.std_sigma()[i])\
            / (np.power(self.coef_aux()[i],2))    

            c3_malla[self.par.start_m()[i] : self.par.end_m()[i]]= c3_coef[i]      

        return c3_malla

    def c4_StDe(self):
        c4_coef = np.empty(self.par.num_materials)
        c4_malla = np.zeros(self.mesh.ncells+1)

        for i in range(self.par.num_materials):
            c4_coef[i] = ((2.0 * self.mesh.dt()/(self.mesh.ddx*self.coef_aux()[i]))\
                *((2*sp.epsilon_0*self.s_par.std_eps_r()[i]*\
                self.s_par.c_eps_H()[i] + self.mesh.dt()*self.s_par.std_sigma()[i] * \
                self.s_par.c_sigma_H()[i])/self.coef_aux()[i])) * \
                math.sqrt(sp.epsilon_0/sp.mu_0)           

            c4_malla[self.par.start_m()[i] : self.par.end_m()[i]]= c4_coef[i]

        return c4_malla         
    
   

class MonteCarlo:
    def __init__(self, mesh, mc_steps):
        self.mesh=mesh
        self.mc_steps=mc_steps


    def gaussiandistribution(self):
        n_materials=self.mesh.par.num_materials

        rnd_epsilon_r=np.zeros((n_materials,self.mc_steps))    
        rnd_sigma=np.zeros((n_materials,self.mc_steps))    

        for i in range(n_materials):
            rnd_epsilon_r[i]=np.random.normal(self.mesh.par.epsilon_r()[i], \
                self.mesh.par.std_eps_r()[i],self.mc_steps)

            rnd_sigma[i]=np.random.normal(self.mesh.par.sigma()[i], \
                self.mesh.par.std_sigma()[i],self.mc_steps)

        return rnd_epsilon_r, rnd_sigma


    def FDTDrun(self):
        n_materials=self.mesh.par.num_materials


        parameters=Materials(set_m, set_s_m)
        #Creo las instancias de malla necesarias
        for i in range(self.mc_steps):

            malla[i]=Mesh(self.mesh.ncells,self.mesh.ddx,parameters[i])









            v