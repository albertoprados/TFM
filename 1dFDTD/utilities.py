import numpy as np
import scipy.constants as sp
from scipy.constants import speed_of_light, epsilon_0, mu_0
import math
from mesh import Mesh

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
        N=len(e1tk1)

        w= ((2*np.pi)/time) * np.arange(1,N+1)         

        return w

class Panels: 
    def __init__(self, mesh, mu_r = 1.0):
        self.mesh=mesh
        self.mu_r = mu_r
        self.n_materials=self.mesh.par.num_materials

    def thickness(self):
        thickness=np.zeros(self.n_materials)
        for i in range(self.n_materials):
            thickness[i]=(self.mesh.par.end_m()[i]-self.mesh.par.start_m()[i])\
                *self.mesh.ddx

        return thickness        

    def eta_0(self):
        return np.sqrt(mu_0/epsilon_0)    

    def epsilon_c(self, omega):
        eps_c=np.zeros(self.n_materials,dtype=complex)

        for i in range(self.n_materials):
            eps_c[i]= self.mesh.par.epsilon_r()[i] * epsilon_0 - \
                (complex(0,1)*self.mesh.par.sigma()[i] / omega)

        return eps_c

    def mu_c(self, omega):
        return self.mu_r * mu_0

    def gamma(self, omega):
        gamma=np.zeros(self.n_materials,dtype=complex)
        for i in range(self.n_materials):
            gamma[i]= complex(0,1) * omega * \
            np.sqrt(self.epsilon_c(omega)[i] * self.mu_c(omega))

        return gamma

    def eta(self, omega):
        eta=np.zeros(self.n_materials,dtype=complex)
        for i in range(self.n_materials):
            eta[i]= np.sqrt(self.mu_c(omega) / self.epsilon_c(omega)[i])

        return eta

    def phi(self, omega):
        phi=np.zeros((self.n_materials,2,2),dtype=complex)
        phi_total=np.array([[1,0],[0,1]])

        gd=np.zeros(self.n_materials,dtype=complex)

        eta=self.eta(omega)
        for i in range(self.n_materials):
            gd[i]= self.gamma(omega)[i] * self.thickness()[i]

            phi[i]=np.array([[np.cosh(gd[i]),      np.sinh(gd[i]) * eta[i]], \
                         [np.sinh(gd[i]) /eta[i], np.cosh(gd[i])      ]])
            
            phi_total=np.matmul(phi_total,phi[i])


        return phi_total

    def _den(self, omega):
        phi = self.phi(omega)
        return phi[0,0]*self.eta_0() + phi[0,1] + phi[1,0]*self.eta_0()**2 + phi[1,1]*self.eta_0()
        
    def T(self, omega):
        return  2*self.eta_0() / self._den(omega)

    def R(self, omega): 
        phi = self.phi(omega)
        return \
            (phi[0,0]*self.eta_0() + phi[0,1] - phi[1,0]*self.eta_0()**2 - phi[1,1]*self.eta_0()) / \
            self._den(omega)

    
    def RyT_freq(self, freq):
        R_freq=np.zeros(len(freq),dtype=complex)
        T_freq=np.zeros(len(freq),dtype=complex)

        for i in range(len(freq)):
            R_freq[i]=self.R(freq[i])
            T_freq[i]=self.T(freq[i])

        return np.abs(R_freq), np.abs(T_freq)        
    

