from viewer import Animator
import pickle


fichero=open("cell10gauss","rb")
cell10gauss=pickle.load(fichero)

Animator().Reflectance_graph(cell10gauss[0], cell10gauss[1], cell10gauss[3], cell10gauss[7], cell10gauss[9], cell10gauss[5])
Animator().Transmittance_graph(cell10gauss[0], cell10gauss[2], cell10gauss[4], cell10gauss[8], cell10gauss[10], cell10gauss[6])
