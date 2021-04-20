from utilities import Source, Utilities
from mesh import Mesh, Materials
from solver import FDTD # Stochastic_FDTD
from viewer import Animator
import numpy as np


#Permitivity,Conductivity,Start_Point,End_Point
material_1=[4,0.04,110,140]
material_2=[3,0.07,80,110]

fat=[16.2,0.214,335,443]
skin=[39,0.43,443,551]
muscle=[55,0.87,551,659]

set_m=[fat,skin,muscle]

#Std_Permitivity,Std_Conductivity,Corr_Eps_E,Corr_Sigma_E,Corr_Eps_H,Corr_Sigma_H
s_material_1=[1,0.002,1,1,1,1]
s_material_2=[0.75,0.003,1,1,1,1]

s_skin=[3.4,0.1,1,1,1,1]
s_fat=[2.7,0.06,1,1,1,1]
s_muscle=[4.6,0.1,1,1,1,1]

set_s_m=[s_fat,s_skin,s_muscle]

parameters=Materials(set_m,set_s_m)

air=[[1,0,0,1059]]
s_air=[[0,0,0,0,0,0]]
parameters_air=Materials(air,s_air)

malla1=Mesh(1059,0.0005,parameters)
malla2=Mesh(1059,0.0005,parameters_air)
pulso=Source('sin',40,12,2e9,malla1,20)


"""
et1k1= FDTD(malla1,pulso).FDTDLoop(5e-9)[0]
e2tk1= FDTD(malla2,pulso).FDTDLoop(5e-9)[0]
et1k2= FDTD(malla1,pulso).FDTDLoop(5e-9)[1]
e2tk2= FDTD(malla2,pulso).FDTDLoop(5e-9)[1]
"""


std_e_film= FDTD(malla1,pulso).FDTDLoop(5e-9)[3]
ex_film=FDTD(malla1,pulso).FDTDLoop(5e-9)[2]

Animator().animationex(ex_film,malla1,'ex')
Animator().animationex(std_e_film,malla1,'std')

"""
r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(5e-9,et1k1)

Animator().fftgraph(freq,r,t)
"""


