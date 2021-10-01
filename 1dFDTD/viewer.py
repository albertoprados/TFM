from matplotlib.lines import Line2D
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from matplotlib import rc

class Animator:

    def animationex(self, exanimation, exanimation_mc, malla, field):

        cb=malla.materials()[1]
       
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set(xlim=(0, malla.ncells), ylim=(-1.2, 1.2))

        x = np.linspace(0, malla.ncells, malla.ncells+1)

        line = ax.plot(x, exanimation[0, :], color='k', lw=2)[0]                
        line2 = ax.plot(x, exanimation_mc[0,:], color='r', lw=2)[0] 

        def animate(i):
            line.set_ydata(exanimation[i, :])
            line2.set_ydata(exanimation_mc[i, :])

        if field=="std":
            plt.ylabel('$\sigma^2$ (E$_x$)', fontsize='14')

        if field=="ex":
            plt.ylabel('E$_x$', fontsize='14')

        plt.plot((0.5 / cb - 1) / 3, 'k--',
                 linewidth=0.75) # The math on cb is just for scaling
        plt.axvline(x=malla.par.start_m()[1])
       

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
    
    def lastsnapshot2(self, ex, ex_2, ex_mc, malla, field):
        cb=malla.materials()[1]
       
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set(xlim=(0, malla.ncells), ylim=(-1.2, 1.2))

        x = np.linspace(0, malla.ncells, malla.ncells+1)

        if field=="std":
            plt.ylabel('$\sigma^2${E$_x$} $(V^2/m^2)$', fontsize='14')
            plt.ylim(0,0.03)
            plt.text(malla.par.start_m()[0] + 54 , 0.025 , 'Fat',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[1] + 54 , 0.025 , 'Skin',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[2] + 54 , 0.025, 'Muscle',
                    horizontalalignment='center')
            plt.plot(x, ex, label='S-FDTD Corr 1',color='limegreen',marker='.', lw='0.5',markersize=2)
            plt.plot(x, ex_2, label='S-FDTD Corr 0.5',color='cornflowerblue',marker='.', lw='0.5',markersize=2)

       

        if field=="ex":
            plt.ylabel('E$_x$ (V/m)', fontsize='14')
            plt.text(malla.par.start_m()[0] + 54 , 0.75 , 'Fat',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[1] + 54 , 0.75 , 'Skin',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[2] + 54 , 0.75, 'Muscle',
                    horizontalalignment='center')
            plt.plot(x, ex, label='FDTD',color='limegreen', marker='.', lw='0.5',markersize=2)
                 
        for i in range(malla.par.num_materials):
            plt.axvline(x=malla.par.start_m()[i], linestyle='--' ,linewidth=0.75, color='k')
        plt.axvline(x=malla.par.end_m()[-1], linestyle='--' ,linewidth=0.75, color='k')
        
       

        plt.xlim(malla.par.start_m()[0]-100,malla.par.end_m()[-1]+100)

        plt.xlabel('Distance')

        plt.subplots_adjust(bottom=0.25, hspace=0.45)

        
        plt.plot(x, ex_mc, label='Monte Carlo (layer)',color='red',marker='.', lw='0.5',markersize=2)

        plt.legend()
        plt.draw()
        plt.show()  

    def lastsnapshot(self, ex, ex_mc, malla, field):
        cb=malla.materials()[1]
       
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set(xlim=(0, malla.ncells), ylim=(-1.2, 1.2))

        x = np.linspace(0, malla.ncells, malla.ncells+1)

        if field=="std":
            plt.ylabel('$\sigma^2${E$_x$} $(V^2/m^2)$', fontsize='14')
            plt.ylim(0,0.03)
            plt.text(malla.par.start_m()[0] + 54 , 0.025 , 'Fat',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[1] + 54 , 0.025 , 'Skin',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[2] + 54 , 0.025, 'Muscle',
                    horizontalalignment='center')
            plt.plot(x, ex, label='S-FDTD',color='limegreen', marker='.', lw='0.5',markersize=2)

        if field=="ex":
            plt.ylabel('E$_x$ (V/m)', fontsize='14')
            plt.text(malla.par.start_m()[0] + 54 , 0.75 , 'Fat',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[1] + 54 , 0.75 , 'Skin',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[2] + 54 , 0.75, 'Muscle',
                    horizontalalignment='center')
            plt.plot(x, ex, label='FDTD',color='limegreen',marker='.', lw='0.5',markersize=3.5)
                 
        for i in range(malla.par.num_materials):
            plt.axvline(x=malla.par.start_m()[i], linestyle='--' ,linewidth=0.75, color='k')
        plt.axvline(x=malla.par.end_m()[-1], linestyle='--' ,linewidth=0.75, color='k')
        
       

        plt.xlim(malla.par.start_m()[0]-100,malla.par.end_m()[-1]+100)

        plt.xlabel('Distance')

        plt.subplots_adjust(bottom=0.25, hspace=0.45)

        
        plt.plot(x, ex_mc, label='Monte Carlo (cell)',color='red',marker='.', lw='0.5',markersize=2)

        plt.legend()
        plt.draw()
        plt.show()    

    def Transmittance_graph(self, freq, t, std_t, t_mc, std_t_mc):
        plt.errorbar(freq, t, yerr= std_t, label='$|T|$ (S-FDTD)', color='darkviolet', capsize=5)
        plt.errorbar(freq, t_mc, yerr= std_t_mc, label='$|T|$ (Monte Carlo (cell))', color='green', capsize=5)

        #plt.ylim(0,1.2)
        
        plt.xlabel('$\omega$ (Hz)')
        #plt.ylabel('R')
        #plt.title('Reflectance in frequency domain')

        plt.legend()
        plt.grid(True)
        plt.show()

    def Reflectance_graph(self, freq, r, std_r, r_mc, std_r_mc):
        plt.errorbar(freq, r, yerr= std_r, label='$|R|$ (S-FDTD)', color='orangered', capsize=5)
        plt.errorbar(freq, r_mc, yerr= std_r_mc, label='$|R|$ (Monte Carlo (cell))', color='blue', capsize=5)

        #plt.ylim(0,1.2)
        
        plt.xlabel('$\omega$ (Hz)')
        #plt.ylabel('R')
        #plt.title('Reflectance in frequency domain')

        plt.legend()
        plt.grid(True)
        plt.show()

    def Reflectance_simple(self, freq, r, std_r, r_panel,r_max,r_min):
        plt.errorbar(freq, r, yerr= std_r)
        plt.plot(freq,r_max, label='rmax')
        plt.plot(freq,r_min, label='rmin')
        plt.plot(freq,r_panel, label='rpanel')
        
        
        plt.xlabel('Frequency w')
        plt.ylabel('R')
        plt.title('Reflectance in frequency domain')

        plt.legend()
        plt.grid(True)
        plt.show()

    def Reflectance_simple1(self, freq, r, t):
        plt.plot(freq,r, label='$|R|$', color='orangered',marker='.', lw='0.5',markersize=2)
        plt.plot(freq,t, label='$|T|$', color='slateblue',marker='.', lw='0.5',markersize=2)
        plt.plot(freq,r**2+t**2, label='$|R|^2+|T|^2$',color='black',marker='.', lw='0.5',markersize=2)
        
        
        plt.xlabel('$\omega$ (Hz)')
        #plt.ylabel('R')
        #plt.title('Reflectance in frequency domain')

        plt.legend()
        plt.grid(True)
        plt.show()

    def Reflectance_simple2(self, freq, r, t, rpanel, tpanel):
        plt.plot(freq,r, label='$|R|$', color='orangered',marker='.', lw='1',markersize=5)
        plt.plot(freq,t, label='$|T|$', color='slateblue',marker='.', lw='1',markersize=5)
        plt.plot(freq,rpanel, label='$Analytic |R|$', color='mediumblue',marker='.', lw='1',markersize=2)
        plt.plot(freq,tpanel, label='$Analytic |T|$', color='darkorange',marker='.', lw='1',markersize=2)
        #plt.plot(freq,r**2+t**2, label='$|R|^2+|T|^2$',color='black',marker='.', lw='1',markersize=2)
        
        
        plt.xlabel('$\omega$ (Hz)')
        #plt.ylabel('R')
        #plt.title('Reflectance in frequency domain')

        plt.legend()
        plt.grid(True)
        plt.show()

    def Reflectance_simple3(self, freq, r, std_r,r_max,r_min):
        plt.errorbar(freq, r, yerr= std_r, label='$|R|$', color='orangered', capsize=5)
        plt.plot(freq,r_max, label='$|R_{min}|$', color='blue',marker='.', lw='0.75',markersize=2)
        plt.plot(freq,r_min, label='$|R_{max}|$', color='black',marker='.', lw='0.75',markersize=2)
        
        
        
        plt.xlabel('$\omega$ (Hz)')
        #plt.ylabel('R')
        #plt.title('Reflectance in frequency domain')

        plt.legend()
        plt.grid(True)
        plt.show()

