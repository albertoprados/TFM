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
set_m=[material1,material2]
set_s_m=[s_material1,s_material2]

set_materials=Materials().material_matrix(set_m)
set_s_materials=Materials().s_material_matrix(set_s_m)

aire=Materials().material(1,0,0,200)
aire_s=Materials().s_material(0,0,0,0,0,0)
airedef=Materials().s_material_matrix([aire])
aire_sdef=Materials().s_material_matrix([aire_s])

malla1=Mesh(200,0.0001,set_materials, set_s_materials)
malla2=Mesh(200,0.0001,airedef,aire_sdef)
pulso=Source('gauss',40,12,20)




et1k1= FDTD(malla1,pulso).FDTDLoop(5e-9)[0]
e2tk1= FDTD(malla2,pulso).FDTDLoop(5e-9)[0]
et1k2= FDTD(malla1,pulso).FDTDLoop(5e-9)[1]
e2tk2= FDTD(malla2,pulso).FDTDLoop(5e-9)[1]



std_e_film= FDTD(malla1,pulso).FDTDLoop(5e-9)[3]
ex_film=FDTD(malla1,pulso).FDTDLoop(5e-9)[2]

Animator().animationex(ex_film,malla1)
Animator().animationex(std_e_film,malla1)


r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(5e-9,et1k1)

Animator().fftgraph(freq,r,t)


#Animator().standard_deviation(std_e,std_h)

