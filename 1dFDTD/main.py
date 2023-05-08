from utilities import Source, Utilities, MultiPanel
from mesh import Mesh, Materials, S_Materials
from solver import FDTD, MonteCarlo 
from viewer import Animator, two_subplots, distribution_plot
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np

#Permitivity,Conductivity,Start_Point,End_Point #0.43, 0.214, 0.87
skin=[39, 0.43, 443, 551]
fat=[16.2, 0.214, 335, 443]
muscle=[55, 0.87, 551, 659]

materiales = [fat,skin,muscle]
par_materiales = Materials(materiales)

#Std_Permitivity,Std_Conductivity,Corr_Eps_E,Corr_Sigma_E,Corr_Eps_H,Corr_Sigma_H
s_skin = [3.4, 0.1, 0.55, 0.55, 0.55, 0.55]
s_fat = [2.7, 0.03, 0.55, 0.55, 0.55, 0.55]
s_muscle = [4.6, 0.1, 0.55, 0.55, 0.55, 0.55]

s_materiales = [s_fat,s_skin,s_muscle]
par_s_materiales = S_Materials(s_materiales)

#--------------------------------
material=[[4,0.04,110,140]]
par_material=Materials(material)

s_material=[[0.5,0.0,1,1,1,1]]
par_s_material=S_Materials(s_material)
#---------------------------------

void=[[1,0,335,659]]
par_void=Materials(void)
s_void=[[0,0,0,0,0,0]]
par_s_void=S_Materials(s_void)

#Tiempo de simulacion
time=11e-9

#--------------------------------
#Parametros de la malla
#Malla 3 materiales
malla1 = Mesh(800, 0.001, par_materiales, par_s_materiales)
malla2 = Mesh(800, 0.001, par_void, par_s_void)

#-------------------------------------------
#Parametros del pulso
pulse_type = 'gauss'
pulso = Source(pulse_type, 700, 160, 2e9, malla1, 20)

#--------------------------------------------
#--------------------------------------------

#Ejecucion y visualizacion
""" ex1_k1, ex1_k2, stde_k1, stde_k2, ex, var_ex, \
    ex_film, var_e_film = FDTD(malla1,pulso,time).FDTDLoop('SFDTD')
ex2_k1, ex2_k2, _, _, _,_,_,_ = FDTD(malla2,pulso,time).FDTDLoop('FDTD')

r, t, freq= Utilities().FFT(ex1_k1,ex2_k1,ex1_k2,ex2_k2,time)
std_r_old,std_t_old= Utilities().FFT_std(stde_k1,stde_k2,ex2_k1,ex2_k2,time)

#r, t, freq= Utilities().FFT(ex1_k1,ex2_k1,ex1_k2,ex2_k2,time)
r_max, t_max, freq= Utilities().FFT(ex1_k1+stde_k1,ex2_k1,ex1_k2+stde_k2,ex2_k2,time)
r_min, t_min, freq= Utilities().FFT(ex1_k1-stde_k1,ex2_k1,ex1_k2-stde_k2,ex2_k2,time)
std_r, std_t= Utilities().FFT_std2(stde_k1,stde_k2,ex1_k1,ex1_k2,ex2_k1,ex2_k2,time)

#Resultado Analitico
r_panel, t_panel=MultiPanel(materiales, malla1).RyT(freq+1) 

Animator().animationex(ex_film,malla1,'ex')
Animator().animationex(var_e_film,malla1,'std')   """


#--------------------------------------------
#--------------------------------------------

#Monte Carlo
mc_steps = 10000
ex_fdtd, _, std_ex_sfdtd= FDTD(malla1, pulso, time).FDTDLoop('SFDTD')[0:3]
mc_simulation = MonteCarlo(malla1, pulso, time, mc_steps)
ex_mc_mean_gauss, ex_var_mc_gauss, _, ex_fixed_tc_gauss = mc_simulation.MC_FDTD_Efield("layer", "gaussian")
ex_mc_mean_uniform, ex_var_mc_uniform, _, ex_fixed_tc_uniform = mc_simulation.MC_FDTD_Efield("layer", "uniform")
ex_mc_mean_gumbel, ex_var_mc_gumbel, _, ex_fixed_tc_gumbel = mc_simulation.MC_FDTD_Efield("layer", "gumbel")
time_array = np.linspace(0, time, mc_simulation.nsteps()+1)   
 
df = pd.DataFrame({"time": time_array, "E_FDTD": ex_fdtd, "Mean_E_Gauss": ex_mc_mean_gauss,
                   "Mean_E_Uniform": ex_mc_mean_uniform, "Mean_E_Gumbel": ex_mc_mean_gumbel,
                   "Var_E_SFDTD": np.square(std_ex_sfdtd), "Var_E_Gaussian": ex_var_mc_gauss, "Var_E_Uniform": ex_var_mc_uniform,
                   "Var_E_Gumbel": ex_var_mc_gumbel})

save_path = r'C:\Users\alber\workspace\Stochastic_FDTD\2023_Stoc_papers\2023_non_normal_distributions\resultados'

file_name = "\{}mc.csv".format(mc_steps)
df.to_csv(save_path+file_name, index=False)

distribution_plot("Electric field for a given cell and instant", "E value", "Pdf",
                  [ex_fixed_tc_gauss, ex_fixed_tc_uniform, ex_fixed_tc_gumbel], 
                  ["MC Gaussian", "MC Uniform", "MC Gumbel"], 
                  "save", save_path, r"\{}mc_histogram".format(mc_steps))
""" 
two_subplots("Random Parameters", "Time (ns)", ["Mean E field", "Var E field"],
             time_array, [ex_fdtd, ex_mc_mean_gauss, ex_mc_mean_uniform, ex_mc_mean_gumbel],
             [np.square(std_ex_sfdtd), ex_var_mc_gauss, ex_var_mc_uniform, ex_var_mc_gumbel],
            ["S-FDTD", "MC Gaussian", "MC Uniform", "MC Gumbel"],
            "save", save_path, r"\{}mc".format(mc_steps))  """


#------------------------------------------
#------------------------------------------




