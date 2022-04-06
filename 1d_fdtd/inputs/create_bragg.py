import numpy as np
import csv
from scipy import constants

# In creating the csv file for the bragg model, we need to set 2 things:
# 1. The materials with different refractive index (1 high and 1 low)
# 2. The wavelength that we want to filter out.

# In this case, we are selecting Plexiglass as the material with a low 
# refractive index and Concrete with moisture (wet concrete) as the 
# material with a high refractive index

f = 5E9
wavelength = constants.c/f
fmax = 5E9
assert f <= fmax
print(f"Frequency: {f} Hz | Wavelength: {wavelength} m")

# Create the necessary parameters for the materials
# For the Concrete 
mu_H = 1
epsilon_H = 10
n_H = np.sqrt(mu_H*epsilon_H)

# For the Plexiglas
mu_L = 1
epsilon_L = 3.4
n_L = np.sqrt(mu_L*epsilon_L)

# Calculating the desired widths of the materials
d_L = wavelength/(4*n_L)
d_H = wavelength/(4*n_H)

print(f"For material with low refractive index")
print(f"n_L = {n_L} | d_L = {d_L} m")
print(f"For material with high refractive index")
print(f"n_H = {n_H} | d_H = {d_H} m")


# Save into csv file


#The number of layers must be even (since the materials are stacked as a pair)
num_layers = 10

with open('bragg_3GHz.csv', 'w', newline='') as csvfile:
    fieldnames = ["layer size","magnetic permeability","electric permittivity","simulation parameters"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames,delimiter=',')

    writer.writeheader()
    for row in range(1,num_layers+1):
        
        #Write the data of low refractive index in odd number while high refractive index in even
        if row % 2 == 0:
            writer.writerow({'layer size': d_H, 'magnetic permeability': mu_H, 'electric permittivity': epsilon_H, 'simulation parameters': 0})
        else:
            #if odd....

            if row == 1:
                # put the freq
                writer.writerow({'layer size': d_L, 'magnetic permeability': mu_L, 'electric permittivity': epsilon_L, 'simulation parameters': fmax})
            elif row == 3:
                # put the number of layers
                writer.writerow({'layer size': d_L, 'magnetic permeability': mu_L, 'electric permittivity': epsilon_L, 'simulation parameters': num_layers})
            else:
                writer.writerow({'layer size': d_L, 'magnetic permeability': mu_L, 'electric permittivity': epsilon_L, 'simulation parameters': 0})