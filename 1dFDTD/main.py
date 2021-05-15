from utilities import Source, Utilities, MultiPanel
from mesh import Mesh, Materials, S_Materials
from solver import FDTD, MonteCarlo 
from viewer import Animator
import numpy as np


#Permitivity,Conductivity,Start_Point,End_Point
skin=[39,0.43,443,551]
fat=[16.2,0.214,335,443]
muscle=[55,0.87,551,659]

materiales=[fat,skin,muscle]
par_materiales=Materials(materiales)
#--------------------------------
material=[[4,0.04,110,140]]
par_material=Materials(material)

#Std_Permitivity,Std_Conductivity,Corr_Eps_E,Corr_Sigma_E,Corr_Eps_H,Corr_Sigma_H
s_skin=[3.4,0.1,1,1,1,1]
s_fat=[2.7,0.06,1,1,1,1]
s_muscle=[4.6,0.1,1,1,1,1]

s_materiales=[s_fat,s_skin,s_muscle]
par_s_materiales=S_Materials(s_materiales)
#---------------------------------
s_material=[[0.5,0.01,1,1,1,1]]
par_s_material=S_Materials(s_material)


s_void=[[0,0,0,0,0,0]]
par_s_void=S_Materials(s_void)


#Tiempo de simulacion
time=2e-8


#Parametros de la malla
#Malla 3 materiales

void=[[1,0,335,659]]
par_void=Materials(void)
malla1=Mesh(800,0.0005,par_materiales,par_s_materiales)
malla2=Mesh(800,0.0005,par_void,par_s_void)

#-------------------------------------------
#Malla 1 material
"""
void=[[1,0,110,140]]
par_void=Materials(void)
malla1=Mesh(200,0.001,par_material,par_s_material)
malla2=Mesh(200,0.001,par_void,par_s_void)
"""

#Parametros del pulso
pulso=Source('gauss',40,12,2e9,malla1,20)

#--------------------------------------------
#--------------------------------------------

#Ejecucion y visualizacion

ex1_k1, ex1_k2, stde_k1, stde_k2 = FDTD(malla1,pulso,time).FDTDLoop('yes')
ex2_k1, ex2_k2, _, _ = FDTD(malla2,pulso,time).FDTDLoop('no')

r, t, freq= Utilities().FFT(ex1_k1,ex2_k1,ex1_k2,ex2_k2,time)
std_r, std_t= Utilities().FFT_std(stde_k1,stde_k2,ex2_k1,ex2_k2,time)
#Resultado Analitico
r_panel, t_panel=MultiPanel(materiales, malla1).RyT(freq+1)

#Visualizacion
Animator().fftgraph(freq,r,t,std_r,std_t,r_panel,t_panel)

#Animator().animationex(ex_film,malla1,'ex')
#Animator().animationex(std_e_film,malla1,'std')

#--------------------------------------------
#--------------------------------------------


#Monte Carlo
"""
#ex_film_mc
#var_e_film_mc
r_mc, t_mc, freq_mc = MonteCarlo(malla1, pulso, time, 10).MC(material,void)
#Film
#Animator().animationex(ex_film_mc,malla1,'ex')
#Animator().animationex(var_e_film_mc,malla1,'var')

#Resultado Analitico
r_panel, t_panel=MultiPanel(material, malla1).RyT(freq_mc+1)

Animator().fftgraph_mc(freq_mc,r_mc,t_mc,r_panel,t_panel)
"""