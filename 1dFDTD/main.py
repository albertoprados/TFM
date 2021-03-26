from utilities import Source, Utilities, Variables_SFDTD
from mesh import Mesh, Materials
from solver import FDTD # Stochastic_FDTD
from viewer import Animator
import numpy as np




material1=Materials().material(4,0.04,110,140)
material2=Materials().material(3,0.02,80,95)
material3=Materials().material(6,0.5,60,70)
set_m=[material1, material2, material3]

set_materials=Materials().material_matrix(set_m)

aire=Materials().material(1,0,0,200)



"""
set_s_var1=Variables_SFDTD(1,0.005,0.8,0.8,0.5,0.5)
set_s_var2=Variables_SFDTD(0.1,0,1,1,1,1)
"""


malla1=Mesh(200,0.0001,set_materials)
malla2=Mesh(200,0.0001,aire)
pulso=Source('gauss',40,12,20)

#cambiar k_ini al pulso
et1k1= FDTD(malla1,pulso).FDTDLoop(5e-9)[0]
e2tk1= FDTD(malla2,pulso).FDTDLoop(5e-9)[0]
et1k2= FDTD(malla1,pulso).FDTDLoop(5e-9)[1]
e2tk2= FDTD(malla2,pulso).FDTDLoop(5e-9)[1]

#std_e= FDTD(malla1,pulso,set_var1,set_s_var1).FDTDLoop()[2]
#std_h= FDTD(malla1,pulso,set_var1,set_s_var1).FDTDLoop()[3]

ex_film=FDTD(malla1,pulso).FDTDLoop(5e-9)[2]
Animator().animationex(ex_film,malla1)

"""
std_film=FDTD(malla1,pulso,set_var1,set_s_var1).FDTDLoop()[4]
Animator().animationex(std_film,malla1,set_var1)
"""

r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(5e-9,et1k1)

Animator().fftgraph(freq,r,t)


#Animator().standard_deviation(std_e,std_h)

