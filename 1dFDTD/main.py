from utilities import Source, Utilities, Variables_FDTD, Variables_SFDTD
from mesh import Mesh
from solver import FDTD
from viewer import Animator



set_1_var=Variables_FDTD(200,0.001,5e-9,4,0.04,110,140)
set_2_var=Variables_FDTD(200,0.001,5e-9,1,0,110,140)


malla1=Mesh(set_1_var)
malla2=Mesh(set_2_var)
pulso=Source('gauss',40,12,20)

#cambiar k_ini al pulso
et1k1= FDTD(malla1,pulso,set_1_var).FDTDLoop()[0]
e2tk1= FDTD(malla2,pulso,set_2_var).FDTDLoop()[0]
et1k2= FDTD(malla1,pulso,set_1_var).FDTDLoop()[1]
e2tk2= FDTD(malla2,pulso,set_2_var).FDTDLoop()[1]

#ex_film=FDTD(malla1,pulso,set_1_var).FDTDLoop()[2]
#Animator().animationex(ex_film,malla1,set_1_var)

r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(set_1_var,et1k1)




Animator().fftgraph(freq,r,t)

