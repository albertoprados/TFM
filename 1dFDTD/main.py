from mesh import Mesh
from solver import FDTD, Utilities, Source
from viewer import Animator
import copy






malla1=Mesh(200,0.001,4,0,110,140)
malla2=Mesh(200,0.001,1,0,110,140)
pulso=Source('gauss',40,12,20)

#cambiar k_ini al pulso
et1k1= FDTD(malla1,pulso,5e-9).FDTDLoop(40,160)[0]
e2tk1= FDTD(malla2,pulso,5e-9).FDTDLoop(40,160)[0]
et1k2= FDTD(malla1,pulso,5e-9).FDTDLoop(40,160)[1]
e2tk2= FDTD(malla2,pulso,5e-9).FDTDLoop(40,160)[1]

#ex_film=FDTD(malla1,pulso,5e-9).FDTDLoop(40,160)[2]
#Animator().animationex(ex_film,malla1)

r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(5e-9,et1k1)




Animator().fftgraph(freq,r,t)

