#!/usr/bin/env python3
#Importing of necessary libraries
from functions import *

#User Input
f_max = float(input("Enter the max frequency of the source excitation (in Hz): "))
mode = int(input("Enter the mode of FDTD Algorithm: \n(1) Basic Algorithm, No source excitation\n(2) Algorithm with Hard Source\n(3) Algorithm with Soft Source\n(4) Algorithm with Soft Source (with PABC)\n(5) Algorithm with TF/SF (with PABC):"))
#NOTE: The assumption in this code is that the computational domain is composed of free-space
#There is no other material present in the domain so the material properties are that of the free space

#Computing Grid Resolution
n_max = 1 #Due to free space
lambda_min = c_0/(f_max*n_max)
N_lambda = 10
delta_lambda = lambda_min/N_lambda
d_min = 0.01
delta_d = d_min/4 #The d_min is assumed to be 0.01 since there is no device in the domain
delta_init = min(delta_lambda,delta_d)
#Since there is currently no critical dimension in the domain, assume it is 0.7
spacers = 15
d_critical = 0.7
N_z = math.ceil(d_critical/delta_init) + spacers
delta_z = d_critical/N_z
print("=====================================================================")
print(f"Number of Yee cells: {N_z} cells\nLength of each cell (Delta_z): {delta_z} m")

#Computing material properties
mu_r = np.ones((1,N_z))    #Due to free space
epsilon_r = np.ones((1,N_z))
n_r = np.sqrt(mu_r[:,0]*epsilon_r[:,0])

#Computing Time step and source excitation
n_bc = 1 #Refractive index at the boundaries (assume free space)
delta_t = (n_bc*delta_z)/(2*c_0)
t_prop = (n_r*N_z*delta_z)/c_0 #time it takes to propagate in the domain
Esrc,Hsrc,t,N_t = gaussian_source(f_max,t_prop,delta_t,delta_z,c_0)
injection_point = 100 #Set this before the device/model in the domain
print("=====================================================================")
print(f"Time step: {delta_t} seconds")
print(f"Number of iterations: {N_t} steps")
print(f"Time vector: {t.shape} [Shape]")
print(f"E-Field (Source):{Esrc.shape}, H-Field (Source): {Hsrc.shape}")
print(f"E={Esrc}, H={Hsrc}")
plot_single(t,Esrc,Hsrc,labels=["Time","Magnitude","Gaussian Pulse Source"])

# Computing the update coefficients
m_E = c_0*delta_t/(epsilon_r*delta_z) #This is assuming that every cell is in free space
m_H = c_0*delta_t/(mu_r*delta_z)

#Field initialization
E = np.zeros((1,N_z))
H = np.zeros((1,N_z))
print("=====================================================================")
print(f"Update coefficients: m_E = {m_E.shape}, m_H = {m_H.shape}")
print(f"Field vectors: Ey: {E.shape}, Hx: {H.shape}")

#Initialize Boundary Terms (For Perfect Absorbing Boundary Conditions)
z_low = [0,0]
z_high = [0,0]
e2 = 0
h2 = 0

#Initialize save matrix
E_plot = np.zeros((N_t,N_z))
H_plot = np.zeros((N_t,N_z))

plot_title = ""

#Check here what mode to use...
if mode == 1: #Basic Algorithm, No source excitation
    #Loop in time for Algorithm
    plot_title = "No Source"
    for i in range(N_t):

        #Update H from E (loop in space)
        for k in range(N_z -1): #Leave out the last cell @ index=N_z-1 for the boundary condition
            H[:,k] = H[:,k] + m_H[:,k]*(E[:,k+1] - E[:,k])
        #Dirichlet Boundary Condition for H at the end of the grid
        H[:,N_z-1] = H[:,N_z-1] + m_H[:,N_z-1]*(0 - E[:,N_z-1])

        #Dirichlet Boundary Condition for E at the start of the grid
        E[:,0] = E[:,0] + m_E[:,0]*(H[:,0]-0)
        #Update E from H (loop in space)
        for k in range(1,N_z):
            E[:,k] = E[:,k] + m_E[:,k]*(H[:,k]-H[:,k-1])

        #Save into matrix
        E_plot[i,:] = E
        H_plot[i,:] = H
        print("=====================================================================")
        print(f"FDTD Algorithm {plot_title}: Successfully computed field values! iteration: {i}/{N_t}")

elif mode == 2: #Algorithm with Hard Source
    plot_title = "Hard Source"

    #Loop in time for Algorithm
    for i in range(N_t):
        #Update H from E (loop in space)
        for k in range(N_z -1): #Leave out the last cell @ index=N_z-1 for the boundary condition
            H[:,k] = H[:,k] + m_H[:,k]*(E[:,k+1] - E[:,k])
        #Dirichlet Boundary Condition for H at the end of the grid
        H[:,N_z-1] = H[:,N_z-1] + m_H[:,N_z-1]*(0 - E[:,N_z-1])

        #Dirichlet Boundary Condition for E at the start of the grid
        E[:,0] = E[:,0] + m_E[:,0]*(H[:,0]-0)
        #Update E from H (loop in space)
        for k in range(1,N_z):
            E[:,k] = E[:,k] + m_E[:,k]*(H[:,k]-H[:,k-1])

        #Insert Source Excitation (Hard Source)
        E[:,injection_point] = Esrc[i]

        #Save into matrix
        E_plot[i,:] = E
        H_plot[i,:] = H
        print("=====================================================================")
        print(f"FDTD Algorithm {plot_title}: Successfully computed field values! iteration: {i}/{N_t}")

elif mode == 3: #Algorithm with Soft Source
    plot_title = "Soft Source"
    #Loop in time for Algorithm
    for i in range(N_t):
        #Update H from E (loop in space)
        for k in range(N_z -1): #Leave out the last cell @ index=N_z-1 for the boundary condition
            H[:,k] = H[:,k] + m_H[:,k]*(E[:,k+1] - E[:,k])
        #Dirichlet Boundary Condition for H at the end of the grid
        H[:,N_z-1] = H[:,N_z-1] + m_H[:,N_z-1]*(0 - E[:,N_z-1])

        #Dirichlet Boundary Condition for E at the start of the grid
        E[:,0] = E[:,0] + m_E[:,0]*(H[:,0]-0)
        #Update E from H (loop in space)
        for k in range(1,N_z):
            E[:,k] = E[:,k] + m_E[:,k]*(H[:,k]-H[:,k-1])

        #Inserting source excitation (Soft Source)
        E[:,injection_point] = E[:,injection_point] + Esrc[i]

        #Save into matrix
        E_plot[i,:] = E
        H_plot[i,:] = H
        print("=====================================================================")
        print(f"FDTD Algorithm {plot_title}: Successfully computed field values! iteration: {i}/{N_t}")

elif mode == 4: #Algorithm with Soft Source (with PABC)
    plot_title = "Soft Source with PABC"
    #Loop in time for Algorithm
    for i in range(N_t):
        #Record H at Boundary
        h2 = z_low.pop(0)
        z_low.append(H[1])

        #Update H from E (loop in space)
        for k in range(N_z -1): #Leave out the last cell @ index=N_z-1 for the boundary condition
            H[:,k] = H[:,k] + m_H[:,k]*(E[:,k+1] - E[:,k])
        # Perfect Absorbing Boundary Condition for H at the end of the grid
        H[:,N_z-1] = H[:,N_z-1] + m_H[:,N_z-1]*(e2 - E[:,N_z-1])

        #Record E at Boundary
        e2 = z_high.pop(0)
        z_high.append(E[:,N_z-1])

        # Perfect Absorbing Boundary Condition for E at the start of the grid
        E[:,0] = E[:,0] + m_E[:,0]*(H[:,0]-h2)
        #Update E from H (loop in space)
        for k in range(1,N_z):
            E[:,k] = E[:,k] + m_E[:,k]*(H[:,k]-H[:,k-1])

        #Inserting source excitation (Soft Source)
        E[:,injection_point] = E[:,injection_point] + Esrc[i]

        #Save into matrix
        E_plot[i,:] = E
        H_plot[i,:] = H
        print("=====================================================================")
        print(f"FDTD Algorithm {plot_title}: Successfully computed field values! iteration: {i}/{N_t}")

elif mode == 5: #Algorithm with TF/SF (with PABC)
    plot_title = "TF/SF with PABC"
    #Loop in time for Algorithm
    for i in range(N_t):

         #Record H at Boundary
        h2 = z_low.pop(0)
        z_low.append(H[1])

        #Update H from E (loop in space)
        for k in range(N_z -1): #Leave out the last cell @ index=N_z-1 for the boundary condition
            H[:,k] = H[:,k] + m_H[:,k]*(E[:,k+1] - E[:,k])
        # Perfect Absorbing Boundary Condition for H at the end of the grid
        H[:,N_z-1] = H[:,N_z-1] + m_H[:,N_z-1]*(e2 - E[:,N_z-1])

        #Handling H source
        H[:,injection_point-1] = H[:,injection_point-1] - m_H[:,injection_point-1]*Esrc[injection_point]


        #Record E at Boundary
        e2 = z_high.pop(0)
        z_high.append(E[:,N_z-1])

        # Perfect Absorbing Boundary Condition for E at the start of the grid
        E[:,0] = E[:,0] + m_E[:,0]*(H[:,0]-h2)
        #Update E from H (loop in space)
        for k in range(1,N_z):
            E[:,k] = E[:,k] + m_E[:,k]*(H[:,k]-H[:,k-1])

        #Handling E source
        E[:,injection_point] = E[:,injection_point] - m_E[:,injection_point]*Hsrc[:,injection_point-1]

        #Save into matrix
        E_plot[i,:] = E
        H_plot[i,:] = H
        print("=====================================================================")
        print(f"FDTD Algorithm {plot_title}: Successfully computed field values! iteration: {i}/{N_t}")


#Plotting the field values....
plot_fields(E_plot,H_plot,N_t,injection_point,title=plot_title,save=False)













