from viewer import Animator
import pickle
from mesh import Mesh, Materials, S_Materials

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
time=8.33e-9


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
#layer10000gauss=[freq, r, t, std_r, std_t, r_panel, t_panel, r_mc, t_mc, std_r_mc, std_t_mc,ex,ex_mc,var_ex,var_ex_mc]
#Validacion6y7
#validation6y7=[r, t, std_r, std_t]
"""
validacion10y11mc=open("cell10000sin","rb")
mc10y11=pickle.load(validacion10y11mc)

validacion10y11sfdtd=open("validation10y11","rb")
sfdtd10y11=pickle.load(validacion10y11sfdtd)
"""

#validation12=[ex, var_ex] #Corr 1
#validation13=[ex, var_ex] #Corr 0.5
validacion12mc=open("layer10000sin","rb")
mc12=pickle.load(validacion12mc)

validacion12sfdtd=open("validation12","rb")
sfdtd12=pickle.load(validacion12sfdtd)

validacion13sfdtd=open("validation13","rb")
sfdtd13=pickle.load(validacion13sfdtd)

Animator().lastsnapshot2(sfdtd12[1],sfdtd13[1],mc12[13],malla1,'std')


#Animator().Reflectance_graph(mc10y11[0], sfdtd10y11[0], sfdtd10y11[2], mc10y11[7], mc10y11[9])
"""
Animator().Reflectance_graph(data[0], data[1], data[3], data[7], data[9], data[5],pulso)
Animator().Transmittance_graph(data[0], data[2], data[4], data[8], data[10], data[6],pulso)
Animator().lastsnapshot(data[11],data[12],malla1,'ex')
Animator().lastsnapshot(data[13],data[14],malla1,'std')
"""