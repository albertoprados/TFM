import numpy as np
import scipy.constants as sp
from math import pi, sin, exp

class Mesh:
    def __init__(self, var):    
        self.var=var
    
    def material(self):

        ca = np.ones(self.var.ncells+1)
        cb = np.ones(self.var.ncells+1) * 0.5
        
        eaf = self.var.dt() * self.var.sigma / (2 * sp.epsilon_0 * self.var.epsilon_r)
        ca[self.var.start_m : self.var.end_m] = (1 - eaf ) / (1 + eaf )
        cb[self.var.start_m : self.var.end_m] = 0.5 / (self.var.epsilon_r * (1 + eaf ))

        return ca, cb

