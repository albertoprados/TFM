from utilities import Source, Utilities
from mesh import Mesh, Materials
from solver import FDTD # Stochastic_FDTD
from viewer import Animator
import numpy as np


#Permitivity,Conductivity,Start_Point,End_Point
material_1=[4,0.04,110,140]
material_2=[3,0.07,80,110]
set_m=[material_1,material_2]
#Std_Permitivity,Std_Conductivity,Corr_Eps_E,Corr_Sigma_E,Corr_Eps_H,Corr_Sigma_H
s_material_1=[1,0.002,1,1,1,1]
s_material_2=[0.75,0.003,1,1,1,1]
set_s_m=[s_material_1,s_material_2]

parameters=Materials(set_m,set_s_m)

air=[[1,0,0,200]]
s_air=[[0,0,0,0,0,0]]
parameters_air=Materials(air,s_air)

malla1=Mesh(200,0.001,parameters)
malla2=Mesh(200,0.001,parameters_air)
pulso=Source('gauss',40,12,20)



et1k1= FDTD(malla1,pulso).FDTDLoop(5e-9)[0]
e2tk1= FDTD(malla2,pulso).FDTDLoop(5e-9)[0]
et1k2= FDTD(malla1,pulso).FDTDLoop(5e-9)[1]
e2tk2= FDTD(malla2,pulso).FDTDLoop(5e-9)[1]



std_e_film= FDTD(malla1,pulso).FDTDLoop(5e-9)[3]
ex_film=FDTD(malla1,pulso).FDTDLoop(5e-9)[2]

Animator().animationex(ex_film,malla1,'ex')
Animator().animationex(std_e_film,malla1,'std')


r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(5e-9,et1k1)

Animator().fftgraph(freq,r,t)



