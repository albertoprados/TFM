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

        w= ((2*np.pi)/time) * np.arange(N)         

        return w

class Panel: 
    def __init__(self, thickness, epsilon_r = 1.0, sigma = 0.0, mu_r = 1.0):
        self.thickness = thickness
        self.epsilon_r = epsilon_r
        self.mu_r = mu_r
        self.sigma = sigma

    def eta_0(self):
        return np.sqrt(mu_0/epsilon_0)    

    def epsilon_c(self, omega):
        return self.epsilon_r*epsilon_0 - complex(0,1)*self.sigma/omega

    def mu_c(self, omega):
        return self.mu_r * mu_0

    def gamma(self, omega):
        return complex(0,1) * omega * \
            np.sqrt(self.epsilon_c(omega) * self.mu_c(omega))

    def eta(self, omega):
        return np.sqrt(self.mu_c(omega) / self.epsilon_c(omega))

    def phi(self, omega):
        gd  = self.gamma(omega) * self.thickness
        eta = self.eta(omega)
        return np.array([[np.cosh(gd),      np.sinh(gd) * eta], \
                         [np.sinh(gd) /eta, np.cosh(gd)      ]])

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
