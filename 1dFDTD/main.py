from utilities import Source, Utilities, Panel
from mesh import Mesh, Materials, S_Materials
from solver import FDTD, MonteCarlo 
from viewer import Animator



#Permitivity,Conductivity,Start_Point,End_Point
skin=[39,0.43,443,551]
fat=[16.2,0.214,335,443]
muscle=[55,0.87,551,659]

set_m=[fat,skin,muscle]
parameters=Materials(set_m)
#--------------------------------
material1=[4,0.04,110,140]
par_m1=Materials([material1])

#Std_Permitivity,Std_Conductivity,Corr_Eps_E,Corr_Sigma_E,Corr_Eps_H,Corr_Sigma_H
s_skin=[3.4,0.1,1,1,1,1]
s_fat=[2.7,0.06,1,1,1,1]
s_muscle=[4.6,0.1,1,1,1,1]

set_s_m=[s_fat,s_skin,s_muscle]
s_parameters=S_Materials(set_s_m)
#---------------------------------
s_material1=[0.5,0.01,1,1,1,1]
s_par_m1=S_Materials([s_material1])



void=[[1,0,0,200]]
s_void=[[0,0,0,0,0,0]]
par_void=Materials(void)
s_par_void=S_Materials(s_void)


#Tiempo de simulacion
time=5e-9

#Parametros de la malla
#Malla 3 materiales
"""
malla1=Mesh(800,0.0005,parameters)
malla2=Mesh(800,0.0005,par_void)
"""
#-------------------------------------------
#Malla 1 material
malla1=Mesh(200,0.0005,par_m1)
malla2=Mesh(200,0.0005,par_void)

#Parametros del pulso
pulso=Source('gauss',40,12,2e9,malla1,20)

#--------------------------------------------
#--------------------------------------------

#Ejecucion y visualizacion
"""
ex1_k1, ex1_k2, ex_film, std_e_film= FDTD(malla1,s_par_m1,pulso,time).FDTDLoop('yes')

ex2_k1, ex2_k2 , _ , _ =FDTD(malla2,s_par_void,pulso,time).FDTDLoop('no')

r, t= Utilities().FFT(ex1_k1,ex2_k1,ex1_k2,ex2_k2)
freq=Utilities().frequency(time,ex1_k1)

Animator().fftgraph(freq,r,t)


Animator().animationex(ex_film1,malla1,'ex')
Animator().animationex(std_e_film,malla1,'std')
"""
#--------------------------------------------
#--------------------------------------------


#Monte Carlo

ex1_k1_mc, ex1_k2_mc, ex_film_mc, var_e_film_mc = MonteCarlo(malla1, [material1], s_par_m1, pulso, time, 10000).FDTDrun()
ex2_k1_mc, ex2_k2_mc, _ , _ = MonteCarlo(malla2, void, s_par_void, pulso, time, 1).FDTDrun()
#Film
Animator().animationex(ex_film_mc,malla1,'ex')
Animator().animationex(var_e_film_mc,malla1,'var')
#FFT

r_mc, t_mc= Utilities().FFT(ex1_k1_mc,ex2_k1_mc, ex1_k2_mc,ex2_k2_mc)
freq_mc=Utilities().frequency(time,ex1_k1_mc)

Animator().fftgraph(freq_mc,r_mc,t_mc)


