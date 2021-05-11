import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from matplotlib import rc

class Animator:

    def animationex(self, exanimation, malla, field):
        self.exanimation=exanimation
        self.malla=malla
        self.field=field

        cb=malla.materials()[1]
       
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set(xlim=(0, malla.ncells), ylim=(-1.2, 1.2))

        x = np.linspace(0, malla.ncells, malla.ncells+1)

        line = ax.plot(x, exanimation[0, :], color='k', lw=2)[0]                

        def animate(i):
            line.set_ydata(exanimation[i, :])

        if field=="std":
            plt.ylabel('$\sigma^2$ (E$_x$)', fontsize='14')

        if field=="ex":
            plt.ylabel('E$_x$', fontsize='14')

        plt.plot((0.5 / cb - 1) / 3, 'k--',
                 linewidth=0.75) # The math on cb is just for scaling


        for i in range(malla.par.num_materials):
            plt.text(170 , 0.25 + 0.2 * i, 'Eps = {}'.format(malla.par.epsilon_r()[i]),
                    horizontalalignment='center')
            plt.text(170 , -0.75 + 0.2 * i, 'Cond = {}'.format(malla.par.sigma()[i]),
                    horizontalalignment='center')
        
        plt.xlabel('FDTD cells')

        plt.subplots_adjust(bottom=0.25, hspace=0.45)


        anim=FuncAnimation(fig, animate, interval=1, frames=20000)
        
        plt.draw()
        plt.show()    
   
    def fftgraph(self, freq, r, t, r_panel, t_panel):
        #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        #rc('text', usetex=True)  
        plt.plot(freq,r, label='R')
        plt.plot(freq,t, label='T')
        
        plt.plot(freq,r_panel, label='R Panel')
        plt.plot(freq,t_panel, label='T Panel')
        #plt.plot(freq,r*r+t*t, label='$R^2+T^2$')

        plt.ylim(-0.2,2.2)
        #plt.xlim(0, 1.5e11)
        
        plt.xlabel('Frequency w')
        plt.ylabel('R&T')
        plt.title('Reflected and transmitted E in frequency domain')

        plt.legend()
        plt.show()

    def fftgraph_mc(self, freq, r, t, r_std, t_std, r_panel, t_panel):

        plt.errorbar(freq, r, yerr= r_std)
        plt.errorbar(freq, t, yerr= t_std)
        
        plt.plot(freq,r_panel, label='R Panel')
        plt.plot(freq,t_panel, label='T Panel')
        #plt.plot(freq,r*r+t*t, label='$R^2+T^2$')

        plt.ylim(-0.2,2.2)
        #plt.xlim(0, 1.5e11)
        
        plt.xlabel('Frequency w')
        plt.ylabel('R&T')
        plt.title('Reflected and transmitted E in frequency domain')

        plt.legend()
        plt.show()
