from utilities import Source, Utilities, Variables_FDTD, Variables_SFDTD
from mesh import Mesh
from solver import FDTD, Stochastic_FDTD
from viewer import Animator



set_var1=Variables_FDTD(200,0.001,5e-9,4,0.04,110,140)
set_var2=Variables_FDTD(200,0.001,5e-9,1,0,110,140)

set_s_var1=Variables_SFDTD(1,0.005,0.8,0.8,0.5,0.5)
set_s_var2=Variables_SFDTD(0.1,0,1,1,1,1)



malla1=Mesh(set_var1)
malla2=Mesh(set_var2)
pulso=Source('gauss',40,12,20)

#cambiar k_ini al pulso
et1k1= FDTD(malla1,pulso,set_var1,set_s_var1).FDTDLoop()[0]
e2tk1= FDTD(malla2,pulso,set_var2,set_s_var2).FDTDLoop()[0]
et1k2= FDTD(malla1,pulso,set_var1,set_s_var1).FDTDLoop()[1]
e2tk2= FDTD(malla2,pulso,set_var2,set_s_var2).FDTDLoop()[1]

std_h= FDTD(malla1,pulso,set_var1,set_s_var1).FDTDLoop()[2]
std_e= FDTD(malla1,pulso,set_var1,set_s_var1).FDTDLoop()[3]

#ex_film=FDTD(malla1,pulso,set_var1).FDTDLoop()[2]
#Animator().animationex(ex_film,malla1,set_var1)

r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(set_var1,et1k1)

Animator().fftgraph(freq,r,t)


Animator().standard_deviation(std_e,std_h)

