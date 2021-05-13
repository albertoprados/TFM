from mesh import Mesh, Materials, S_Materials
from utilities import Utilities
import numpy as np
import copy
import math
import scipy.constants as sp


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

    def FDTDLoop(self, Stochastic_Analysis):
        self.Stochastic_Analysis=Stochastic_Analysis

        ex=np.zeros(self.mesh.ncells+1)
        hy=np.zeros(self.mesh.ncells)
        ex_old=np.zeros(self.mesh.ncells+1)

        #Saving values for film
        #ex_film=np.empty((self.nsteps()+1,self.mesh.ncells+1))

        #Fourier transform
        ex_k1=np.empty(self.nsteps()+1)
        ex_k2=np.empty(self.nsteps()+1)
        

        #Standard Deviation
        if Stochastic_Analysis == 'yes':
            std_e=np.zeros(self.mesh.ncells+1)
            std_h=np.zeros(self.mesh.ncells)
            std_e_old=std_e=np.zeros(self.mesh.ncells+1)

        #std_e_film=np.empty((self.nsteps()+1,self.mesh.ncells+1))    


        ca, cb, cc = self.mesh.materials()
        k1, k2=self.mesh.FFTpoints()

        for time_step in range(self.nsteps()):
            ex_old=copy.deepcopy(ex)
            ex[1:-1] = ca[1:-1] * ex[1:-1] + cb[1:-1] * (hy[:-1] - hy[1:])
            
            #Guardo los valores a representar
            #ex_film[time_step][:]=ex[:]
            
            #Guardo los valores para calcular la transformada
            ex_k1[time_step]=ex[k1]
            ex_k2[time_step]=ex[k2]
           
            ex[self.pulse.k_ini] += 0.5 * self.pulse.pulse(time_step)  
            
            self.boundarymur(ex,ex_old)  
            
            hy[:] = hy[:] + cc * (ex[:-1] - ex[1:])
               

            if Stochastic_Analysis == 'yes':
                self.Include_SFDTD_Analysis(std_h,std_e,ex_old,hy,std_e_old)
                #std_e_film[time_step][:]=std_e[:]
            
            
            t= time_step + 1/2
            hy[self.pulse.k_ini] += 0.25 * self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.25 * self.pulse.pulse(t)   

            

        return ex_k1, ex_k2, #ex_film, np.power(std_e_film ,2)




class MonteCarlo:
    def __init__(self, mesh, pulse, time, mc_steps):
        self.mesh=mesh
        self.s_par=self.mesh.s_par
        self.pulse=pulse
        self.time=time
        self.mc_steps=mc_steps
        self.n_materials=self.mesh.par.num_materials
        

    def Gaussian_Pdf(self):
        rnd_epsilon_r=np.zeros((self.n_materials,self.mc_steps))       
        rnd_sigma=np.zeros((self.n_materials,self.mc_steps))  

        for i in range(self.n_materials):
            rnd_epsilon_r[i]=np.random.normal(self.mesh.par.epsilon_r()[i], \
                self.s_par.std_eps_r()[i],self.mc_steps)
            rnd_sigma[i]=np.random.normal(self.mesh.par.sigma()[i], \
                self.s_par.std_sigma()[i],self.mc_steps)
    

        return rnd_epsilon_r, rnd_sigma


    def pdf(self,x,mu,var):
        #Distribución eponencial
        return (1.0/mu)*np.exp(-(1.0/mu)*x)

    def Metropolis(self,delta,markov_0):
        n=self.mc_steps
        #Valores medios
        eps_r=self.mesh.par.epsilon_r()
        sigma=self.mesh.par.sigma()
        #Desviaciones
        std_eps=self.s_par.std_eps_r()
        std_sigma=self.s_par.std_sigma()

        rnd_eps=np.zeros((self.n_materials,n))
        rnd_sigma=np.zeros((self.n_materials,n))
        #Inicio de la cadena de Markov
        for i in range(self.n_materials):
            rnd_eps[i][0]=markov_0[i][0]
            rnd_sigma[i][0]=markov_0[i][1]

        for i in range(self.n_materials):
            for k in range(n):
                #Propuesta 
                y=rnd_eps[i][k]+np.random.uniform(-delta[i][0],delta[i][0])
                z=rnd_sigma[i][k]+np.random.uniform(-delta[i][1],delta_sigma[i][1])

                if np.random.rand() < min(1,self.pdf(y,eps_r[i],std_eps[i]) \
                    /self.pdf(rnd_eps[i][k],eps_r[i],std_eps[i])):
                    rnd_eps[i][k+1]=y
                else:
                    rnd_eps[i][k+1]=rnd_eps[i][k]    

                if np.random.rand() < min(1,self.pdf(z,sigma[i],std_sigma[i]) \
                    /self.pdf(rnd_sigma[i][k],sigma[i],std_sigma[i])):
                    rnd_sigma[i][k+1]=z
                else:
                    rnd_sigma[i][k+1]=rnd_sigma[i][k]      

        return rnd_eps, rnd_sigma


    def MC(self, materiales, void):
        
        nsteps=FDTD(self.mesh, self.pulse, self.time).nsteps()
        ncells=self.mesh.ncells

        rnd_epsilon_r, rnd_sigma=self.Gaussian_Pdf()
        
        #Definicion de vectores medios
        #ex_film_avg=np.zeros((nsteps+1,ncells+1))
        #ex_film_var=np.zeros((nsteps+1,ncells+1))

        #E auxiliar FFT
        malla_aux=Mesh(ncells, self.mesh.ddx, Materials(void), self.s_par)
        ex2_k1, ex2_k2=FDTD(malla_aux, self.pulse, self.time).FDTDLoop('no')       

        ex_k1, ex_k2=FDTD(self.mesh, self.pulse, self.time).FDTDLoop('no')

        freq =  Utilities().FFT(ex_k1,ex2_k1,ex_k2,ex2_k2,self.time)[2]
        R=np.zeros((len(freq),self.mc_steps))
        T=np.zeros((len(freq),self.mc_steps))

        """
        R_avg=np.zeros(len(freq))
        T_avg=np.zeros(len(freq))
        R_std=np.zeros(len(freq))
        T_std=np.zeros(len(freq))
        """

        for k in range(self.mc_steps):
            for i in range(self.n_materials):
                materiales[i][0]=rnd_epsilon_r[i][k]
                materiales[i][1]=rnd_sigma[i][k]
        
            malla=Mesh(ncells,self.mesh.ddx,Materials(materiales),self.s_par)

            ex_k1, ex_k2=FDTD(malla, self.pulse, self.time).FDTDLoop('no')
            
            #Film
            #ex_film_avg += ex_film    
            #ex_film_var += np.power(ex_film,2)

            #RyT
            R[:,k] = Utilities().FFT(ex_k1,ex2_k1,ex_k2,ex2_k2,self.time)[0]
            T[:,k] = Utilities().FFT(ex_k1,ex2_k1,ex_k2,ex2_k2,self.time)[1]

            """
            R_avg += R
            T_avg += T
            R_std += np.power(R,2)
            T_std += np.power(T,2)
            """

            if (k % 10)==0:
                print(k)

        
        #ex_film_avg = ex_film_avg / self.mc_steps
        #ex_film_var = (ex_film_var /self.mc_steps) - (np.power(ex_film_avg,2))
        """
        R_avg = R_avg / self.mc_steps
        T_avg = T_avg / self.mc_steps
        R_std = (R_std /self.mc_steps) - (np.power(R_avg,2))
        T_std = (T_std /self.mc_steps) - (np.power(T_avg,2))
        """

        return   R, T, freq, #ex_film_avg, ex_film_var



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
    
   











            