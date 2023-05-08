from matplotlib.lines import Line2D
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import (LinearLocator, FormatStrFormatter,
MultipleLocator, AutoMinorLocator)
from matplotlib import rc
from scipy.stats import gumbel_r, norm

class Animator:

    def animationex(self, exanimation, malla, field):

        cb=malla.materials()[1]
       
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set(xlim=(0, malla.ncells), ylim=(-1.2, 1.2))

        x = np.linspace(0, malla.ncells, malla.ncells+1)

        line = ax.plot(x, exanimation[0, :], color='k', lw=2)[0]                
        #line2 = ax.plot(x, exanimation_mc[0,:], color='r', lw=2)[0] 

        def animate(i):
            line.set_ydata(exanimation[i, :])
            #line2.set_ydata(exanimation_mc[i, :])

        if field=="std":
            plt.ylabel(r'$\sigma^2$ (E$_x$)', fontsize='14')

        if field=="ex":
            plt.ylabel(r'E$_x$', fontsize='14')

        plt.plot((0.99 / cb - 1) / 3, 'k--',
                 linewidth=0.75) # The math on cb is just for scaling
        plt.axvline(x=malla.par.start_m()[1], ymin=0.5, linestyle = '--', 
                                            color = "k", linewidth = 0.75)
        plt.axvline(x=malla.par.start_m()[2], ymin=0.5, linestyle = '--', 
                                            color = "k", linewidth = 0.75)

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
            plt.ylabel(r'$\sigma^2${E$_x$} $(V^2/m^2)$', fontsize='14')
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
            plt.ylabel(r'E$_x$ (V/m)', fontsize='14')
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
            plt.ylabel(r'$\sigma^2${E$_x$} $(V^2/m^2)$', fontsize='14')
            plt.ylim(0,0.03)
            plt.text(malla.par.start_m()[0] + 54 , 0.025 , 'Fat',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[1] + 54 , 0.025 , 'Skin',
                    horizontalalignment='center')
            plt.text(malla.par.start_m()[2] + 54 , 0.025, 'Muscle',
                    horizontalalignment='center')
            plt.plot(x, ex, label='S-FDTD',color='limegreen', marker='.', lw='0.5',markersize=2)

        if field=="ex":
            plt.ylabel(r'E$_x$ (V/m)', fontsize='14')
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

def two_subplots(title, xlabel, ylabel, time_array: list, mean_ex: list,
                  var_ex: list, labels: list, 
                  show_or_save: str = "save", output_path: str = "",
                  file_name: str = ""):

    figure, axs = plt.subplots(2, sharex= True, figsize = (9, 8))
    figure.suptitle(title, fontsize = 18)
    
    for i in range(2):
        axs[i].set_ylabel(ylabel[i], fontsize = 15)
        axs[i].xaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i].tick_params(axis='both', which='major', labelsize=10)
        axs[i].xaxis.offsetText.set_fontsize(10)
        axs[i].grid(which='both')
        if i == 1:
            axs[i].set_xlabel(xlabel, fontsize = 15)
    
    for i in range(len(mean_ex)):
        axs[0].plot(time_array * 1e9, mean_ex[i], 
            label = labels[i])
        axs[0].legend(fontsize = 8) 
        axs[1].plot(time_array * 1e9, var_ex[i], 
            label = labels[i])
        axs[1].legend(fontsize = 8) 

    if show_or_save == "show":
        plt.tight_layout(pad=0.7)
        plt.show()
    else:
        figure.tight_layout(pad=0.7)
        figure.savefig(output_path + file_name)

def distribution_plot(title, xlabel, ylabel, ex_fixed_time_cell: list,
                      labels: list, 
                      show_or_save: str = "save", output_path: str = "",
                      file_name: str = ""):
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    for i in range(len(ex_fixed_time_cell)):
        plt.hist(ex_fixed_time_cell[i], bins = 200, density = True, label = labels[i])
        plt.legend()

    plt.tight_layout(pad=0.7)
    
    if show_or_save == "show":
        plt.show()
    else:
        plt.savefig(output_path + file_name)