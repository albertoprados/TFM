import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import rc

class Animator:

    def animationex(self, exanimation, malla, field):
        self.exanimation=exanimation
        self.malla=malla
        self.field=field

        cb=malla.materials()[1]
       
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set(xlim=(0, 200), ylim=(-1.2, 1.2))

        x = np.linspace(0, 200, 201)

        line = ax.plot(x, exanimation[0, :], color='k', lw=2)[0]                

        def animate(i):
            line.set_ydata(exanimation[i, :])

        if field=="std":
            plt.ylabel('StdE$_x$', fontsize='14')

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
    

        anim=FuncAnimation(fig, animate, interval=1, frames=1500)
        
        plt.draw()
        plt.show()    

    def fftgraph(self, freq, r, t):
        self.freq=freq
        self.r=r
        self.t=t

        #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        #rc('text', usetex=True)

        plt.plot(freq,r*r, label='R')
        plt.plot(freq,t*t, label='T')
        plt.plot(freq,r*r+t*t, label='$R^2+T^2$')

        plt.ylim(0,2)
        plt.xlim(0, 5e10)
        
        plt.xlabel('Frequency w')
        plt.ylabel('R&T')
        plt.title('Reflected and transmitted E in frequency domain')

        plt.legend()
        plt.show()
    
    def standard_deviation(self, std_e, std_h):
        self.std_e=std_e
        self.std_h=std_h

        #dom=len(std_e)

        plt.plot(std_e,label='std_e')
        plt.plot(std_h,label='std_h')

        plt.legend()
        plt.show()

   