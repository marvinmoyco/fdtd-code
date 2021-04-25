from scipy import constants
import numpy as np
import math
import matplotlib.pyplot as plt
#Defining the constants
c_0 = constants.speed_of_light
mu_0 = constants.mu_0
epsilon_0 = constants.epsilon_0

def gaussian_source(f_max,t_prop,delta_t,delta_z,c_0):
    """
    Gaussian Pulse Source that is used for injecting source in 1D FDTD.
    """
    #Set source permittivity and permeability
    mu_src =1 #these parameters should be the material permittivity/permeability at the grid position of the source injection
    epsilon_src = 1

    #Computing source parameters
    tau = 0.5/f_max
    t_0 = 6*tau
    T = 12*tau + 5*t_prop #Total time of simulation
    N_t = math.ceil(T/delta_t) #Number of time steps
    t = np.linspace(0,N_t*delta_t,N_t)

    n_src = np.sqrt(epsilon_src*mu_src)
    deltaT = (n_src*delta_z/2*c_0)+(delta_t/2)

    A = -np.sqrt(epsilon_src/mu_src)
    x_E = (t - t_0)/tau
    x_H = (t-t_0 + deltaT)/tau

    Esrc = np.exp(-np.power(x_E,2))
    Hsrc = A*np.exp(-np.power(x_H,2))
    return Esrc,Hsrc,t,N_t

def plot_single(x=0,y1=0,y2=0,size=[20,13],labels=["X-Axis","Y-Axis","Title"]):
    """
    Plot a vector in a single time step (snapshot in time)
    """
    plt.figure(figsize = size)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.title(labels[2])
    plt.plot(x,y1)
    plt.plot(x,y2)
    plt.show()

def plot_fields(E_plot,H_plot,N_t,injection_point,title="",save=False):
    plt.ion()
    fig = plt.figure(1,[10,6])
    ax = fig.add_subplot(111)
    plt.ylabel("Value")
    
    _,x = E_plot.shape
    for x in range(x):
        plt.axvline(x=x,ymin=0.3,ymax=0.7,color="forestgreen",linestyle = "-")
    plt.axhline(y=2,color="forestgreen",linestyle = "-")
    plt.axhline(y=-2,color="forestgreen",linestyle = "-")
    plt.axvline(x=injection_point, color = "grey", linestyle="--")
    plt.xlabel("z (Space)")
    plt.ylim((-5,5))
    plt.xlim((0,x))
    lineE, = ax.plot(E_plot[0,:])
    lineH, = ax.plot(H_plot[0,:])
    fig.show()

    for i in range(1,N_t):
        plt.ylim(-5,5)
        plt.legend(handles = [lineE,lineH],labels=["Electric Field","Magnetic Field"])
        print(f"Currently plotting Iteration step: {i}/{N_t}")
        plt.title(f'FDTD 1-D {title} | Iteration step: {i}/{N_t}')
        if save == True:
            plt.savefig(f"photos/1d-fdtd{i}.jpeg")
        lineE.set_ydata(E_plot[i,:])
        lineH.set_ydata(H_plot[i,:])
        fig.canvas.draw()
        fig.canvas.flush_events()


#The FDTD Algorithm used here are in the Ey/Hx Modes
def algo_no_source(Ey,Hx,N_z,m_Ey,m_Hx):
    """
    This Basic FDTD Algorithm uses the Dirichlet Boundary Condition and has no source excitation.
    """

    #Update Hx from Ey (loop in space)
    for k in range(N_z -1): #Leave out the last cell @ index=N_z-1 for the boundary condition
        Hx[k] = Hx[k] + m_Hx*(Ey[k+1] - Ey[k])
    #Dirichlet Boundary Condition for Hx at the end of the grid
    Hx[N_z-1] = Hx[N_z-1] + m_Hx*(0 - Ey[N_z-1])

    #Dirichlet Boundary Condition for Ey at the start of the grid
    Ey[0] = Ey[0] + m_Ey*(Hx[0]-0)
    #Update Ey from Hx (loop in space)
    for k in range(1,N_z):
        Ey[k] = Ey[k] + m_Ey*(Hx[k]-Hx[k-1])

    return Ey,Hx


def algo_soft_pabc(Ey,Hx,N_z,m_Ey,m_Hx,z_low,z_high,e2,h2):
    """
    This Basic FDTD Algorithm uses the Perfect Absorbing Boundary Condition and has Soft Source excitation.
    """

    #Record H at Boundary
    h2 = z_low.pop(0)
    z_low.append(Hx[1])

    #Update Hx from Ey (loop in space)
    for k in range(N_z -1): #Leave out the last cell @ index=N_z-1 for the boundary condition
        Hx[k] = Hx[k] + m_Hx*(Ey[k+1] - Ey[k])
    # Perfect Absorbing Boundary Condition for Hx at the end of the grid
    Hx[N_z-1] = Hx[N_z-1] + m_Hx*(e2 - Ey[N_z-1])

    #Record E at Boundary
    e2 = z_high.pop(0)
    z_high.append(Ey[N_z-1])

    # Perfect Absorbing Boundary Condition for Ey at the start of the grid
    Ey[0] = Ey[0] + m_Ey*(Hx[0]-h2)
    #Update Ey from Hx (loop in space)
    for k in range(1,N_z):
        Ey[k] = Ey[k] + m_Ey*(Hx[k]-Hx[k-1])

    return Ey,Hx,z_low,z_high

