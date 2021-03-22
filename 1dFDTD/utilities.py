import numpy as np
import scipy.constants as sp
import math


class Variables_FDTD:
    def __init__(self,ncells,ddx,time,epsilon_r,sigma,start_m,end_m):
        self.ncells=ncells
        self.ddx=ddx
        self.time=time
        self.epsilon_r= epsilon_r
        self.sigma= sigma
        self.start_m= start_m
        self.end_m= end_m 

    def dt(self):
        return self.ddx/(2*sp.c)    
    
    def FFTpoints(self):
        return self.start_m - 20, self.end_m + 20


class Variables_SFDTD:
    def __init__(self,std_eps,std_sigma,c_eps_E,c_sigma_E,c_eps_H,c_sigma_H):   
        self.std_eps=std_eps
        self.std_sigma=std_sigma
        self.c_eps_E= c_eps_E
        self.c_sigma_E= c_sigma_E
        self.c_eps_H= c_eps_H
        self.c_sigma_H= c_sigma_H




class Source:
    def __init__(self, sourcetype, t_0, s_0, k_ini):
        self.sourcetype=sourcetype
        self.t_0=t_0
        self.s_0=s_0
        self.k_ini=k_ini

    def pulse(self, time):
        
        self.time=time
        
        if self.sourcetype == 'gauss':
            pulse = math.exp(-0.5*( (self.t_0 - time) / self.s_0 )**2)
        
        return pulse



class Utilities:

    def FFT(self,e1tk1_total,e2tk1,e1tk2,e2tk2):
        
        #Hay que cancelar la parte incidente
        e1tk1_reflected = e1tk1_total - e2tk1  
        
        e1wk1=np.fft.fft(e1tk1_reflected)
        e2wk1=np.fft.fft(e2tk1)

        e1wk2=np.fft.fft(e1tk2)
        e2wk2=np.fft.fft(e2tk2)
    
        R=np.abs(e1wk1) / np.abs(e2wk1)
        T=np.abs(e1wk2) / np.abs(e2wk2)
        
        
        return  R, T
    

    def frequency(self,time,e1tk1):
        self.time=time
        
        N=len(e1tk1)

        freq= (1.0/self.time.time) * np.arange(N)         

        return freq
