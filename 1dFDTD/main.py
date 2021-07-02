from utilities import Source, Utilities, MultiPanel
from mesh import Mesh, Materials, S_Materials
from solver import FDTD, MonteCarlo 
from viewer import Animator
import pickle


#Permitivity,Conductivity,Start_Point,End_Point
skin=[39,0.43,443,551]
fat=[16.2,0.214,335,443]
muscle=[55,0.87,551,659]

skin_max1=[39.5,0.45,443,551]
fat_max1=[16.4,0.224,335,443]
muscle_max1=[55.6,0.9,551,659]

skin_min1=[38.5,0.41,443,551]
fat_min1=[16.0,0.204,335,443]
muscle_min1=[54.4,0.84,551,659]

materiales=[fat,skin,muscle]

materiales_max1=[fat_max1,skin_max1,muscle_max1]
materiales_min1=[fat_min1,skin_min1,muscle_min1]

par_materiales=Materials(materiales)

par_materiales_max1=Materials(materiales_max1)
par_materiales_min1=Materials(materiales_min1)
#--------------------------------
material=[[4,0.04,110,140]]

material_max1=[[4.5,0.05,110,140]]
material_min1=[[3.5,0.03,110,140]]
par_material=Materials(material)
par_material_max=Materials(material_max1)
par_material_min=Materials(material_min1)
#Std_Permitivity,Std_Conductivity,Corr_Eps_E,Corr_Sigma_E,Corr_Eps_H,Corr_Sigma_H
s_skin=[3.4,0.10,1,1,1,1]
s_fat=[2.7,0.06,1,1,1,1]
s_muscle=[4.6,0.10,1,1,1,1]

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

malla_max=Mesh(800,0.0005,par_materiales_max1,par_s_materiales)
malla_min=Mesh(800,0.0005,par_materiales_min1,par_s_materiales)

#-------------------------------------------
#Malla 1 material
"""
void=[[1,0,110,140]]
par_void=Materials(void)
malla1=Mesh(200,0.001,par_material,par_s_material)
malla2=Mesh(200,0.001,par_void,par_s_void)

malla_max=Mesh(200,0.001,par_material_max,par_s_material)
malla_min=Mesh(200,0.001,par_material_min,par_s_material)
"""
#Parametros del pulso
pulso=Source('gauss',40,12,2e9,malla1,20)

#--------------------------------------------
#--------------------------------------------

#Ejecucion y visualizacion

ex1_k1, ex1_k2, stde_k1, stde_k2, ex_film, var_e_film = FDTD(malla1,pulso,time).FDTDLoop('SFDTD')
ex2_k1, ex2_k2, _, _, _,_ = FDTD(malla2,pulso,time).FDTDLoop('FDTD')

r, t, freq= Utilities().FFT(ex1_k1,ex2_k1,ex1_k2,ex2_k2,time)
std_r, std_t= Utilities().FFT_std(stde_k1,stde_k2,ex2_k1,ex2_k2,time)

#Resultado Analitico
r_panel, t_panel=MultiPanel(materiales, malla1).RyT(freq+1)


#--------------------------------------------
#--------------------------------------------

#Monte Carlo

#ex_film_mc
#var_e_film_mc
r_mc, t_mc, std_r_mc, std_t_mc, ex_film_mc, var_e_film_mc = MonteCarlo(malla1, pulso, time, 10).M_FDTD('cell')

#Film
#Animator().animationex(ex_film_mc,malla1,'ex')
#Animator().animationex(var_e_film_mc,malla1,'var')
#Resultado Analitico
#r_panel, t_panel=MultiPanel(materiales, malla1).RyT(freq_mc+1)
#Animator().fftgraph(freq_mc, r_mc, t_mc, r_std_mc, t_std_mc, r_panel, t_panel)
#All coeficients
#Animator().fftgraph_mc(freq_mc,r_mc,t_mc,r_panel,t_panel)

#------------------------------------------
#------------------------------------------
"""
layer10gauss=[freq, r, t, std_r, std_t, r_panel, t_panel, r_mc, t_mc, std_r_mc, std_t_mc]
fichero=open("layer10gauss","wb")
pickle.dump(layer10gauss,fichero)
fichero.close()
del(fichero)
"""
Animator().Reflectance_simple(freq, r_mc, std_r_mc)

#Visualizacion
"""
Animator().Reflectance_graph(freq, r, std_r, r_mc, std_r_mc, r_panel)
Animator().Transmittance_graph(freq, t, std_t, t_mc, std_t_mc, t_panel)

Animator().animationex(ex_film,ex_film_mc,malla1,'ex')
Animator().animationex(var_e_film,var_e_film_mc,malla1,'std')
"""