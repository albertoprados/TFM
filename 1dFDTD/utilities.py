import numpy as np
import scipy.constants as sp
from scipy.constants import speed_of_light, epsilon_0, mu_0
import math
from mesh import Mesh, Materials

class Source:
    def __init__(self, sourcetype, delay, spread, freq, malla, k_ini):
        self.sourcetype=sourcetype
        self.delay=delay
        self.spread=spread
        self.k_ini=k_ini
        self.freq=freq
        self.malla=malla

    def pulse(self, time):
        
        self.time=time
        
        if self.sourcetype == 'gauss':
            pulse = math.exp(-0.5*( (self.delay - time) / self.spread )**2)
        
        if self.sourcetype == 'sin':
            pulse = math.sin(2.0*np.pi*self.freq*self.malla.dt()*time)

        return pulse
    


class Utilities:
    
    def FFT(self,e1tk1_total,e2tk1,e1tk2,e2tk2,time):
        #Frequency
        freq_min=0
        freq_max=1.2e11

        N=len(e1tk1_total)
        freq= ((2*np.pi)/time) * np.arange(0,N)  

        #Hay que cancelar la parte incidente
        e1tk1_reflected = e1tk1_total - e2tk1  
        
        e1wk1=np.fft.fft(e1tk1_reflected)
        e2wk1=np.fft.fft(e2tk1)

        e1wk2=np.fft.fft(e1tk2)
        e2wk2=np.fft.fft(e2tk2)

        e1wk1=e1wk1[(freq_min <= freq) & (freq < freq_max)]
        e2wk1=e2wk1[(freq_min <= freq) & (freq < freq_max)]
        e1wk2=e1wk2[(freq_min <= freq) & (freq < freq_max)]
        e2wk2=e2wk2[(freq_min <= freq) & (freq < freq_max)]
        freq=freq[(freq_min <= freq) & (freq < freq_max)]

        R=np.abs(e1wk1) / np.abs(e2wk1)
        T=np.abs(e1wk2) / np.abs(e2wk2)
        
        return  R, T, freq
    
    def FFT_std(self, stde_k1, stde_k2, e2_k1, e2_k2, time):
        #Frequency
        freq_min=0
        freq_max=1.2e11

        N=len(stde_k1)
        freq= ((2*np.pi)/time) * np.arange(0,N)  
        """
        stdew_k1=sum(np.power(stde_k1,2))
        stdew_k2=sum(np.power(stde_k2,2))
       
        for i in range(len(stde_k1)):
            if (i%100)==0:
                        print(i)
            for j in range(len(stde_k1)):
                if i != j:
                    stdew_k1 += stde_k1[i]*stde_k1[j]
                    stdew_k2 += stde_k2[i]*stde_k2[j]
        stdew_k1=np.sqrt(stdew_k1)       
        stdew_k2=np.sqrt(stdew_k2)          
        """
        
        
        stdew_k1=np.fft.fft(stde_k1)
        stdew_k2=np.fft.fft(stde_k2)
        
        """
        stdew_k1=np.zeros(N,dtype=complex)
        stdew_k2=np.zeros(N,dtype=complex)
        fourier_coef=np.zeros(N,dtype=complex)
        coef=np.exp((-complex(0,1)*2*np.pi)/N)

        for k in range(N):
            fourier_corr_k1=0 
            fourier_corr_k2=0
            if (k%100)==0:
                print(k)
            for i in range(N):
                if (i%100)==0:
                    print(i)
                a_i= np.power(coef,k*i)   
                fourier_coef[i]=np.power(a_i,2)
                for j in range(N):
                    if j>i :
                        a_j= np.power(coef,k*j) 
                        fourier_corr_k1 += 2 * a_i * a_j * stde_k1[i] * stde_k1[j]
                        fourier_corr_k2 += 2 * a_i * a_j * stde_k2[i] * stde_k2[j]


            stdew_k1[k]=np.dot(fourier_coef,np.power(stde_k1,2)) + fourier_corr_k1    
            stdew_k2[k]=np.dot(fourier_coef,np.power(stde_k2,2)) + fourier_corr_k2    
            
        stdew_k1=np.sqrt(stdew_k1)
        stdew_k2=np.sqrt(stdew_k2)
        """

        e2w_k1=np.fft.fft(e2_k1)
        e2w_k2=np.fft.fft(e2_k2)

        stdew_k1=stdew_k1[(freq_min <= freq) & (freq < freq_max)]
        stdew_k2=stdew_k2[(freq_min <= freq) & (freq < freq_max)]
        e2w_k1=e2w_k1[(freq_min <= freq) & (freq < freq_max)]
        e2w_k2=e2w_k2[(freq_min <= freq) & (freq < freq_max)]
        freq=freq[(freq_min <= freq) & (freq < freq_max)]

        Std_R= np.abs(stdew_k1) / np.abs(e2w_k1)
        Std_T= np.abs(stdew_k2) / np.abs(e2w_k2)   

        return Std_R, Std_T 


class MultiPanel: 
    def __init__(self, materials, mesh):
        self.materials = materials
        self.mesh=mesh
        self.eta_0 = np.sqrt(mu_0/epsilon_0)   

    def Phi(self, omega):
        uno=np.ones(len(omega))
        zero=np.zeros(len(omega))
        phi_total=np.array([[uno,zero],[zero,uno]])

        for i in range(len(self.materials)):
            phi=Panel(self.materials[i],self.mesh).phi(omega)

            phi_total=np.einsum('ijn,jkn->ikn', phi_total, phi)

        return phi_total        

    def _den(self, omega):
        Phi = self.Phi(omega)
        return Phi[0,0]*self.eta_0 + Phi[0,1] + Phi[1,0]*self.eta_0**2 + Phi[1,1]*self.eta_0
        
    def RyT(self, omega):
        Phi = self.Phi(omega)
        T= 2*self.eta_0 / self._den(omega)
        R=(Phi[0,0]*self.eta_0 + Phi[0,1] - Phi[1,0]*self.eta_0**2 - Phi[1,1]*self.eta_0) / \
            self._den(omega)

        return  np.abs(R), np.abs(T)

    


class Panel: 
    def __init__(self, material, mesh, mu_r = 1.0):
        self.material=material
        self.mesh=mesh
        self.mu_r = mu_r
        self.par=Materials([material])
        self.eta_0 = np.sqrt(mu_0/epsilon_0)

    def thickness(self):
        return (self.par.end_m()-self.par.start_m())\
                *self.mesh.ddx       

    def epsilon_c(self, omega):
        return self.par.epsilon_r()*epsilon_0 - complex(0,1)*self.par.sigma()/omega

    def mu_c(self):
        return self.mu_r * mu_0

    def gamma(self, omega):
        return complex(0,1) * omega * \
            np.sqrt(self.epsilon_c(omega) * self.mu_c())

    def eta(self, omega):
        return np.sqrt(self.mu_c() / self.epsilon_c(omega))

    def phi(self, omega):
        gd  = self.gamma(omega) * self.thickness()
        eta = self.eta(omega)
        return np.array([[np.cosh(gd),      np.sinh(gd) * eta], \
                         [np.sinh(gd) /eta, np.cosh(gd)      ]])

    def _den(self, omega):
        phi = self.phi(omega)
        return phi[0,0]*self.eta_0 + phi[0,1] + phi[1,0]*self.eta_0**2 + phi[1,1]*self.eta_0
        
    def T(self, omega):
        return  2*self.eta_0 / self._den(omega)

    def R(self, omega): 
        phi = self.phi(omega)
        return \
            (phi[0,0]*self.eta_0 + phi[0,1] - phi[1,0]*self.eta_0**2 - phi[1,1]*self.eta_0) / \
            self._den(omega)