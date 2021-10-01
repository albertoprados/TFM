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
time=8.6e-9


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
pulse_type='sin'
pulso=Source(pulse_type,40,12,2e9,malla1,20)

#--------------------------------------------
#--------------------------------------------

#Ejecucion y visualizacion

ex1_k1, ex1_k2, stde_k1, stde_k2, ex, var_ex = FDTD(malla1,pulso,time).FDTDLoop('SFDTD')
ex2_k1, ex2_k2, _, _, _,_ = FDTD(malla2,pulso,time).FDTDLoop('FDTD')

r, t, freq= Utilities().FFT(ex1_k1,ex2_k1,ex1_k2,ex2_k2,time)
std_r_old,std_t_old= Utilities().FFT_std(stde_k1,stde_k2,ex2_k1,ex2_k2,time)

r, t, freq= Utilities().FFT(ex1_k1,ex2_k1,ex1_k2,ex2_k2,time)
r_max, t_max, freq= Utilities().FFT(ex1_k1+stde_k1,ex2_k1,ex1_k2+stde_k2,ex2_k2,time)
r_min, t_min, freq= Utilities().FFT(ex1_k1-stde_k1,ex2_k1,ex1_k2-stde_k2,ex2_k2,time)
std_r, std_t= Utilities().FFT_std2(stde_k1,stde_k2,ex1_k1,ex1_k2,ex2_k1,ex2_k2,time)

#Resultado Analitico
r_panel, t_panel=MultiPanel(materiales, malla1).RyT(freq+1)

#Primera validacion r^2+t^2
#Animator().Reflectance_simple1(freq, r, t)
#Segunda validacion clase Panel
#Animator().Reflectance_simple2(freq, r, t, r_panel, t_panel)
#Tercera validacion 3 materiales con la clase panel
#Animator().Reflectance_simple2(freq, r, t, r_panel, t_panel)
#Cuarta validacion stdR pulso gaussiano y sinusoidal
Animator().Reflectance_simple3(freq, r, std_r,r_max,r_min)
"""
validation05=[ex, var_ex]
fichero=open("validation05","wb")
pickle.dump(validation05,fichero)
fichero.close()
del(fichero)
"""
"""
validation12=[ex, var_ex]
fichero=open("validation12","wb")
pickle.dump(validation12,fichero)
fichero.close()
del(fichero)
"""

"""
validation10y11=[r, t, std_r, std_t]
fichero=open("validation10y11","wb")
pickle.dump(validation10y11,fichero)
fichero.close()
del(fichero)
"""

"""
Animator().Reflectance_simple(freq, r, std_r, r_panel,r_max,r_min)
Animator().Reflectance_simple(freq, t, std_t, t_panel,t_max,t_min)
Animator().Reflectance_simple(freq, r, std_r2, r_panel,r_max,r_min)
Animator().Reflectance_simple(freq, t, std_t2, t_panel,t_max,t_min)
"""
#Animator().animationex(ex_film,malla1,'ex')
#Animator().animationex(var_e_film,malla1,'ex')
#--------------------------------------------
#--------------------------------------------

#Monte Carlo

#r_mc, t_mc, std_r_mc, std_t_mc, ex_mc, var_ex_mc = MonteCarlo(malla1, pulso, time, 10000).M_FDTD('layer')

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
layer10000gauss=[freq, r, t, std_r, std_t, r_panel, t_panel, r_mc, t_mc, std_r_mc, std_t_mc,ex,ex_mc,var_ex,var_ex_mc]
fichero=open("layer10000gauss","wb")
pickle.dump(layer10000gauss,fichero)
fichero.close()
del(fichero)
"""

#Animator().Reflectance_simple(freq, r_mc, std_r_mc)
#Animator().animationex(ex_film,ex_film_mc,malla1,'ex')
#Animator().animationex(var_e_film,var_e_film_mc,malla1,'std')

#Visualizacion
"""
Animator().Reflectance_graph(freq, r, std_r2, r_mc, std_r_mc, r_panel, pulse_type)
Animator().Transmittance_graph(freq, t, std_t2, t_mc, std_t_mc, t_panel,pulse_type)
Animator().lastsnapshot(ex,ex_mc,malla1,'ex')

"""
validacion05=open("validation05","rb")
var05=pickle.load(validacion05)

validacion12mc=open("layer10000sin","rb")
mc12=pickle.load(validacion12mc)

Animator().lastsnapshot2(var_ex, var05[1], mc12[14] ,malla1,'std')