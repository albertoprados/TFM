from utilities import Source, Utilities
from mesh import Mesh, Materials, S_Materials
from solver import FDTD, MonteCarlo 
from viewer import Animator



#Permitivity,Conductivity,Start_Point,End_Point
skin=[39,0.43,443,551]
fat=[16.2,0.214,335,443]
muscle=[55,0.87,551,659]

set_m=[fat,skin,muscle]
parameters=Materials(set_m)

#Std_Permitivity,Std_Conductivity,Corr_Eps_E,Corr_Sigma_E,Corr_Eps_H,Corr_Sigma_H
s_skin=[3.4,0.1,1,1,1,1]
s_fat=[2.7,0.06,1,1,1,1]
s_muscle=[4.6,0.1,1,1,1,1]

set_s_m=[s_fat,s_skin,s_muscle]
s_parameters=S_Materials(set_s_m)


void=[[1,0,0,1059]]
s_void=[[0,0,0,0,0,0]]
par_void=Materials(void)
s_par_void=S_Materials(s_void)


malla1=Mesh(1059,0.0005,parameters)
malla2=Mesh(1059,0.0005,par_void)
pulso=Source('sin',40,12,2e9,malla1,20)


"""
et1k1= FDTD(malla1,pulso).FDTDLoop(5e-9)[0]
e2tk1= FDTD(malla2,pulso).FDTDLoop(5e-9)[0]
et1k2= FDTD(malla1,pulso).FDTDLoop(5e-9)[1]
e2tk2= FDTD(malla2,pulso).FDTDLoop(5e-9)[1]
"""



_, _, ex_film, std_e_film= FDTD(malla1,parameters,s_parameters,pulso).FDTDLoop()

Animator().animationex(ex_film,malla1,'ex')
Animator().animationex(std_e_film,malla1,'std')



"""
r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(5e-9,et1k1)

Animator().fftgraph(freq,r,t)
"""
"""
prueba=MonteCarlo(malla1, 100).gaussiandistribution()[1]

print(prueba)
"""