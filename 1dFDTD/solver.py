from mesh import Mesh, Materials, S_Materials
from utilities import Utilities
import numpy as np
import copy
import math
import sys
import cmath
import scipy.constants as sp
import warnings
warnings.filterwarnings("error")


class FDTD:
    def __init__(self, mesh, pulse, time):
        self.mesh=mesh
        self.pulse=pulse
        self.s_par=self.mesh.s_par
        self.time=time

    def boundarymur(self,ex,ex_old):
        ncells, dt, ddx= self.mesh.ncells, self.mesh.dt(), self.mesh.ddx

        c_bound=(sp.c*dt-ddx)/(sp.c*dt+ddx)

        ex[0]=ex_old[1] + c_bound * (ex[1]-ex_old[0])
        ex[ncells]=ex_old[ncells-1] + c_bound * (ex[ncells-1]-ex_old[ncells])


    def Include_SFDTD_Analysis(self,std_h,std_e,e,h,std_e_old):
        std_e_old=copy.deepcopy(std_e) 
        
        Stochastic_FDTD(self.mesh).StandardDeviation_E(std_h,std_e,e,h)
        Stochastic_FDTD(self.mesh).BoundaryCondition(std_e,std_e_old)
        Stochastic_FDTD(self.mesh).StandardDeviation_H(std_h,std_e)
        
        
    def nsteps(self):
        return int(self.time / self.mesh.dt())  


    def FDTDLoop(self, Computation_Method, distribution: str = ""):
        
        ex=np.zeros(self.mesh.ncells+1)
        hy=np.zeros(self.mesh.ncells)
        ex_old=np.zeros(self.mesh.ncells+1)

        #Saving values for film
        #ex_film=np.zeros((self.nsteps()+1, self.mesh.ncells+1))

        #Fourier transform
        ex_k1=np.zeros(self.nsteps()+1)
        ex_k2=np.zeros(self.nsteps()+1)
        std_e_k1=np.zeros(self.nsteps()+1)
        std_e_k2=np.zeros(self.nsteps()+1)

        #Standard Deviation
        std_e=np.zeros(self.mesh.ncells+1)
        #std_e_film=np.zeros((self.nsteps()+1,self.mesh.ncells+1))    
        if Computation_Method == 'SFDTD':
            std_h=np.zeros(self.mesh.ncells)
            std_e_old=np.zeros(self.mesh.ncells+1)

        if Computation_Method == 'MFDTD':
            ca, cb, cc = self.mesh.cellsproperties(distribution)
        else:    
            ca, cb, cc = self.mesh.materials()

        k1, k2 = self.mesh.FFTpoints()

        for time_step in range(self.nsteps()+1):
            ex_old=copy.deepcopy(ex)

            ex[self.pulse.k_ini] += self.pulse.pulse(time_step)  

            ex[1:-1] = ca[1:-1] * ex[1:-1] + cb[1:-1] * (hy[:-1] - hy[1:])
            
            #Guardo los valores a representar
            #ex_film[time_step][:]=ex[:]
            
            #Guardo los valores para calcular la transformada
            ex_k1[time_step]=ex[k1]
            ex_k2[time_step]=ex[k2]
        
            self.boundarymur(ex,ex_old)  

            t= time_step + 1/2
            hy[self.pulse.k_ini] += 0.5 * self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.5 * self.pulse.pulse(t)  
             
            hy[:] = hy[:] + cc * (ex[:-1] - ex[1:])
            
            if Computation_Method == 'SFDTD':
                self.Include_SFDTD_Analysis(std_h,std_e,ex_old,hy,std_e_old)
                std_e_k1[time_step]=std_e[k1]
                std_e_k2[time_step]=std_e[k2]
                #std_e_film[time_step][:]=std_e[:]

        return ex_k1, ex_k2, std_e_k1, std_e_k2, ex, np.power(std_e ,2)



def draw_samples_from_distribution(mu: float, sigma:float, 
                                   n_samples: int, distribution: str):
    any_negative = True
    if distribution == "gaussian":
        while any_negative:
            samples = np.random.normal(mu, sigma, n_samples)
            any_negative = np.any(samples < 0)
            if any_negative == False:
                break
            
    elif distribution == "uniform":
        a = mu - np.sqrt(3) * sigma
        b = mu + np.sqrt(3) * sigma
        while any_negative:
            samples = np.random.uniform(a, b, n_samples)
            any_negative = np.any(samples < 0)
            if any_negative == False:
                break

    elif distribution == "gumbel":
        location = mu - 0.57721 * (np.sqrt(6)/np.pi) * sigma
        scale = np.sqrt(6) * sigma / np.pi
        while any_negative:
            samples = np.random.gumbel(location, scale, n_samples)
            any_negative = np.any(samples < 0)
            if any_negative == False:
                break

    return samples



class MonteCarlo:
    def __init__(self, mesh, pulse, time, mc_steps):
        self.mesh=mesh
        self.s_par=self.mesh.s_par
        self.pulse=pulse
        self.time=time
        self.mc_steps=mc_steps
        self.n_materials=self.mesh.par.num_materials
        self.materiales=self.mesh.par.materiales


    def nsteps(self):
        return int(self.time / self.mesh.dt())   
        

    def Pdf_layer(self, distribution: str):
        rnd_epsilon_r=np.zeros((self.n_materials, self.mc_steps))       
        rnd_sigma=np.zeros((self.n_materials, self.mc_steps))  

        for i in range(self.n_materials):
            rnd_epsilon_r[i] = draw_samples_from_distribution(self.mesh.par.epsilon_r()[i], \
                self.s_par.std_eps_r()[i], self.mc_steps, distribution)
            
            rnd_sigma[i] = draw_samples_from_distribution(self.mesh.par.sigma()[i], \
                self.s_par.std_sigma()[i], self.mc_steps, distribution)
            
        return rnd_epsilon_r, rnd_sigma
    

    def Layer_Method(self, k, distribution: str):
        ncells=self.mesh.ncells
        rnd_epsilon_r, rnd_sigma = self.Pdf_layer(distribution)

        aux_materiales = copy.deepcopy(self.materiales)

        for i in range(self.n_materials):
            aux_materiales[i][0]=rnd_epsilon_r[i][k]
            aux_materiales[i][1]=rnd_sigma[i][k]

        malla = Mesh(ncells, self.mesh.ddx, Materials(aux_materiales), self.s_par)
        ex_k1, ex_k2, _, _, ex, _= FDTD(malla, self.pulse, self.time).FDTDLoop('FDTD')

        return ex_k1, ex_k2, ex


    def MC_FDTD_Efield(self, layer_or_cell: str, distribution: str):

        ex_avg=np.zeros(self.nsteps()+1)
        ex_var=np.zeros(self.nsteps()+1)

        if layer_or_cell == "layer":
            for k in range(self.mc_steps):
                ex = self.Layer_Method(k, distribution)[0]
                ex_avg += ex  
                ex_var += np.power(ex,2)

                if k%100==0:
                    print(k)

        elif layer_or_cell == "cell":
            for k in range(self.mc_steps):
                ex = FDTD(self.mesh, self.pulse, 
                            self.time).FDTDLoop('MFDTD', distribution)[0]
                ex_avg += ex  
                ex_var += np.power(ex,2)

                if k%100==0:
                    print(k)

        coef_std= (self.mc_steps)/(self.mc_steps-1)
        ex_avg = ex_avg / self.mc_steps
        ex_var = (ex_var /(self.mc_steps-1)) - coef_std * (np.power(ex_avg,2))
        ex_var[ex_var < 0] = 0

        try:
            std_ex = np.sqrt(ex_var)
        except RuntimeWarning:
            print("Negative value")
        
        return ex_avg, ex_var, std_ex


    def M_FDTD(self,layer_or_cell):
        
        nsteps=FDTD(self.mesh, self.pulse, self.time).nsteps()
        ncells=self.mesh.ncells

        if layer_or_cell == 'layer':
            rnd_epsilon_r, rnd_sigma=self.Gaussian_Pdf_layer()
        
        #Definicion de vectores medios
        #ex_film_avg=np.zeros((nsteps+1,ncells+1))
        #ex_film_var=np.zeros((nsteps+1,ncells+1))

        #E_x last time step average
        ex_avg=np.zeros(ncells+1)
        ex_var=np.zeros(ncells+1)

        #E auxiliar FFT
        void=[[1,0, self.mesh.par.start_m()[0], self.mesh.par.end_m()[-1]]]
        malla_aux=Mesh(ncells, self.mesh.ddx, Materials(void), self.s_par)
        ex2_k1, ex2_k2, _, _, _,_,_,_=FDTD(malla_aux, self.pulse, self.time).FDTDLoop('FDTD')       

        ex_k1, ex_k2,_,_,_,_,_,_=FDTD(self.mesh, self.pulse, self.time).FDTDLoop('FDTD')

        freq =  Utilities().FFT(ex_k1,ex2_k1,ex_k2,ex2_k2,self.time)[2]

        #All coefficients R&T
        #R=np.zeros((len(freq),self.mc_steps))
        #T=np.zeros((len(freq),self.mc_steps))
    
        R_avg=np.zeros(len(freq))
        T_avg=np.zeros(len(freq))
        R_sum=np.zeros(len(freq))
        T_sum=np.zeros(len(freq))
        R_sum2=np.zeros(len(freq))
        T_sum2=np.zeros(len(freq))
        R_std=np.zeros(len(freq))
        T_std=np.zeros(len(freq))
        
        for k in range(self.mc_steps):
            """
            ex_k1=np.zeros(nsteps+1)
            ex_k2=np.zeros(nsteps+1)
            """
          
            if layer_or_cell == 'layer':
                ex_k1, ex_k2, ex = self.Layer_Method(k, rnd_epsilon_r, rnd_sigma)
            if layer_or_cell == 'cell':
                ex_k1, ex_k2, _, _, ex, _ = FDTD(self.mesh, self.pulse, self.time).FDTDLoop('MFDTD')      


            #Film
            ex_avg += ex  
            ex_var += np.power(ex,2)

            #RyT all coeficients
            #R[:,k] = Utilities().FFT(ex_k1,ex2_k1,ex_k2,ex2_k2,self.time)[0]
            #T[:,k] = Utilities().FFT(ex_k1,ex2_k1,ex_k2,ex2_k2,self.time)[1]

            R, T, _ = Utilities().FFT(ex_k1,ex2_k1,ex_k2,ex2_k2,self.time)
           
            R_sum += R
            T_sum += T
            R_sum2 += np.power(R,2)
            T_sum2 += np.power(T,2)

            if k%1000==0:
                print(k)

        coef_std= self.mc_steps/(self.mc_steps-1)

        ex_avg = ex_avg / self.mc_steps
        ex_var = (ex_var /(self.mc_steps-1)) - coef_std * (np.power(ex_avg,2))
        
        R_avg = R_sum / self.mc_steps
        T_avg = T_sum / self.mc_steps
        
        R_std = np.sqrt((R_sum2 /(self.mc_steps-1)) - coef_std * (np.power(R_avg,2)))
        T_std = np.sqrt((T_sum2 /(self.mc_steps-1)) - coef_std * (np.power(T_avg,2)))
        
        return R_avg, T_avg, R_std, T_std, ex_avg, ex_var   

       

class Stochastic_FDTD:
    def __init__(self, mesh):
        self.mesh=mesh
        self.s_par=self.mesh.s_par
               
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
        c_aux = np.empty(self.mesh.par.num_materials)

        for i in range(self.mesh.par.num_materials):     
            c_aux[i] = 2 * sp.epsilon_0 * self.mesh.par.epsilon_r()[i] + \
                self.mesh.dt() * self.mesh.par.sigma()[i]   
        
        return c_aux

    def c1_StDe(self):
        c1_coef = np.empty(self.mesh.par.num_materials)
        c1_malla = np.ones(self.mesh.ncells+1)

        for i in range(self.mesh.par.num_materials):
            c1_coef[i] = ( (2 * sp.epsilon_0 * self.mesh.par.epsilon_r()[i] - \
                self.mesh.dt() * self.mesh.par.sigma()[i]) \
                / self.coef_aux()[i] ) 
                

            c1_malla[self.mesh.par.start_m()[i] : self.mesh.par.end_m()[i]]= c1_coef[i]    
                             
        return c1_malla
    
    def c2_StDe(self):
        c2_coef = np.empty(self.mesh.par.num_materials)
        c2_malla = np.ones(self.mesh.ncells+1) * (self.mesh.dt()/ \
            (self.mesh.ddx*sp.epsilon_0)) * math.sqrt(sp.epsilon_0/sp.mu_0)
                
        for i in range(self.mesh.par.num_materials): 
            c2_coef[i]= (2.0 * self.mesh.dt() / (self.mesh.ddx * self.coef_aux()[i])) \
                        * math.sqrt(sp.epsilon_0/sp.mu_0)

            c2_malla[self.mesh.par.start_m()[i] : self.mesh.par.end_m()[i]]= c2_coef[i]     
                      
        return c2_malla

    def c3_StDe(self):
        c3_coef = np.empty(self.mesh.par.num_materials)
        c3_malla =  np.zeros(self.mesh.ncells+1)

        for i in range(self.mesh.par.num_materials):
            c3_coef[i] = 4 * self.mesh.dt() * sp.epsilon_0 * \
            (self.mesh.par.sigma()[i]*self.s_par.c_eps_E()[i]*self.s_par.std_eps_r()[i]-\
            self.mesh.par.epsilon_r()[i]*self.s_par.c_sigma_E()[i]*self.s_par.std_sigma()[i])\
            / (np.power(self.coef_aux()[i],2))    

            c3_malla[self.mesh.par.start_m()[i] : self.mesh.par.end_m()[i]]= c3_coef[i]      

        return c3_malla

    def c4_StDe(self):
        c4_coef = np.empty(self.mesh.par.num_materials)
        c4_malla = np.zeros(self.mesh.ncells+1)

        for i in range(self.mesh.par.num_materials):
            c4_coef[i] = ((2.0 * self.mesh.dt()/(self.mesh.ddx*self.coef_aux()[i]))\
                *((2*sp.epsilon_0*self.s_par.std_eps_r()[i]*\
                self.s_par.c_eps_H()[i] + self.mesh.dt()*self.s_par.std_sigma()[i] * \
                self.s_par.c_sigma_H()[i])/self.coef_aux()[i])) * \
                math.sqrt(sp.epsilon_0/sp.mu_0)           

            c4_malla[self.mesh.par.start_m()[i] : self.mesh.par.end_m()[i]]= c4_coef[i]

        return c4_malla         
    
   











            