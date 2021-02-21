import numpy as np
import copy
from math import pi, sin, exp
import scipy.constants as sp
import doctest

class FDTD:
    def __init__(self, mesh, pulse, time):
        self.mesh=mesh
        self.pulse=pulse
        self.time=time

    def boundarymur(self, ex, bl, bh): 
        ex[0] = bl.pop(0)
        bl.append(ex[1])      

        ex[self.mesh.ncells] = bh.pop(0)
        bh.append(ex[self.mesh.ncells-1])


    def FDTDLoop(self,k1,k2):
        
        dt=self.mesh.ddx / (2*sp.c)
        nsteps= int(self.time  / dt)

        # COMENTAR: Mejor quitar nsteps, no guardar siempre todo...
        ex=np.zeros(self.mesh.ncells+1)
        hy=np.zeros(self.mesh.ncells+1)
        
        ex_save_k1=np.empty(nsteps+1)
        ex_save_k2=np.empty(nsteps+1)

        #ex_save_film=np.empty((nsteps+1,self.mesh.ncells+1))
        
        ca=self.mesh.material()[0][1:-1]
        cb=self.mesh.material()[1][1:-1]

        bl = [0, 0]
        bh = [0, 0]
       
        for time_step in range(1, nsteps + 1):

            # Calculate the Ex field, for cycle using slice notation
            ex[1:-1] = ca * ex[1:-1] + cb * (hy[:-2] - hy[1:-1])
            
            #Guardo los valores a representar
            #ex_save_film[time_step][:]=copy.deepcopy(ex[:])
            
            #Guardo los valores para calcular la transformada
            ex_save_k1[time_step]=copy.deepcopy(ex[k1])
            ex_save_k2[time_step]=copy.deepcopy(ex[k2])
           
            ex[self.pulse.k_ini] +=  0.5*self.pulse.pulse(time_step) 
            
            #Condiciones de contorno
            self.boundarymur(ex,bl,bh)  
            
            # Update h field
            hy[:-1] = hy[:-1] + 0.5 * (ex[:-1] - ex[1:])   

            t= time_step+1/2
            hy[self.pulse.k_ini] += 0.25* self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.25* self.pulse.pulse(t)   

                       

        return ex_save_k1, ex_save_k2,  #ex_save_film




class Source:
    def __init__(self, sourcetype, t_0, s_0, k_ini):
        self.sourcetype=sourcetype
        self.t_0=t_0
        self.s_0=s_0
        self.k_ini=k_ini

    def pulse(self, time):
        
        self.time=time
        
        if self.sourcetype == 'gauss':
            pulse = exp(-0.5*( (self.t_0 - time) / self.s_0 )**2)
        
        return pulse




#Clase para la Trasformada Rápida de Fourier
# COMENTAR: Esto es mas un namespace que una clase. 
# COMENTAR: Cuanto menos estado, mejor
class Utilities:

    def FFT(self,e1tk1,e2tk1,e1tk2,e2tk2):
        
        """Devuelve la transformada de fourier

        >>> Utilities().FFT(e1tk1,e2tk1,e1tk2,e2tk2)
        True
        """  
        #Hay que cancelar la parte incidente
        e1tk1 = e1tk1 - e2tk1
        
        e1wk1=np.fft.fft(e1tk1)
        e2wk1=np.fft.fft(e2tk1)

        e1wk2=np.fft.fft(e1tk2)
        e2wk2=np.fft.fft(e2tk2)
    
        R=np.abs(e1wk1) / np.abs(e2wk1)
        T=np.abs(e1wk2) / np.abs(e2wk2)

        if all(i < 1 for i in R):
            boolvar= True
        else:
            boolvar= False

        return boolvar, R, T
    

    def frequency(self,time,e1tk1):

        N=len(e1tk1)

        freq= (1.0/time) * np.arange(N)         

        return freq

#doctest.testmod()   