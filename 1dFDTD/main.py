from utilities import Source, Utilities
from mesh import Mesh, Materials
from solver import FDTD # Stochastic_FDTD
from viewer import Animator
import numpy as np




material1=Materials().material(4,0.04,110,140)
s_material1=Materials().s_material(1,0.005,1,1,1,1)
material2=Materials().material(3,0.5,80,109)
s_material2=Materials().s_material(0.75,0.1,1,1,1,1)
#material3=Materials().material(6,0.5,60,70)
set_m=[material1, material2]
set_s_m=[s_material1,s_material2]

set_materials=Materials().material_matrix(set_m)
set_s_materials=Materials().s_material_matrix(set_s_m)

aire=Materials().material(1,0,0,200)
aire_s=Materials().s_material(0,0,0,0,0,0)


"""
set_s_var1=Variables_SFDTD(1,0.005,0.8,0.8,0.5,0.5)
set_s_var2=Variables_SFDTD(0.1,0,1,1,1,1)
"""


malla1=Mesh(200,0.0001,set_materials)
malla2=Mesh(200,0.0001,aire)
pulso=Source('gauss',40,12,20)

#cambiar k_ini al pulso
et1k1= FDTD(malla1,pulso,set_materials,set_s_materials).FDTDLoop(5e-9)[0]
e2tk1= FDTD(malla2,pulso,aire,aire_s).FDTDLoop(5e-9)[0]
et1k2= FDTD(malla1,pulso,set_materials,set_s_materials).FDTDLoop(5e-9)[1]
e2tk2= FDTD(malla2,pulso,aire,aire_s).FDTDLoop(5e-9)[1]



std_e_film= FDTD(malla1,pulso,set_materials,set_s_materials).FDTDLoop(5e-9)[3]
ex_film=FDTD(malla1,pulso,set_materials,set_s_materials).FDTDLoop(5e-9)[2]

Animator().animationex(std_e_film,malla1)

"""
Animator().animationex(std_film,malla1,set_var1)
"""

r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(5e-9,et1k1)

Animator().fftgraph(freq,r,t)


#Animator().standard_deviation(std_e,std_h)

