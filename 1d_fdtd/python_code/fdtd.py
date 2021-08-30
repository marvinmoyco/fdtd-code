#Library imports
import csv
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt

c_0 = constants.speed_of_light
mu_0 = constants.mu_0
epsilon_0 = constants.epsilon_0



# Class definitions for 1D FDTD



class Simulation:

    def __init__(self):
        """
        At object creation, the input file (csv) where the simulation parameters are set should be included.
        """
        #Initialization of relevant variables
        self.sim_param = {"dz":0,"dt":0,"Nz":0,"Nt":0,"fmax":0,"source_type":"no_source"}
        #Storage variable for the input file parameters
        self.input = {}
        self.comp_domain = {}
        self.spacer_cells = 0
        self.source = 0

    def init(self,input_path=None,spacers=0,injection_point=0): 
        """
        This function will read the input file and calculate dz and dt of the whole computational domain as well as the creation of the computational domain.
        """
        print('====================================================================================================================================')
        print("Initialization of an FDTD Simulation....")
        #Read the input file
        if input_path != None:
            with open(input_path,mode='r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    for key_val in row.items():
                        key,value = key_val
                        if float(value) != 0: #Skip if the value is zero
                            try:
                                self.input[key].append(float(value))
                            except KeyError:
                                self.input[key] = [float(value)]
                        elif key == 'simulation parameters': #Make sure that the simulation parameters are unaffected
                            try:
                                self.input[key].append(float(value))
                            except KeyError:
                                self.input[key] = [float(value)]
        else: 
            print('==========================================================================')
            print('Error: No input file in the argument.')
            return

        #Check if the shape of the input file are valid
        try:
            assert len(self.input['layer size']) == self.input['simulation parameters'][2]
            assert len(self.input['magnetic permeability']) == self.input['simulation parameters'][2]
            assert len(self.input['electric permittivity']) == self.input['simulation parameters'][2]
            assert self.input['simulation parameters'][1] == 0 or self.input['simulation parameters'][1] == 1
        except AssertionError:
            print("Error: Input file do not have a valid size (e.g. nmodel is not equal to the number of rows of either: (1) layer size, (2) magnetic permeability, or (3) electric permeability.\nOr the source type entered is invalid.")
            return
        print('====================================================================================================================================')
        print(f"Contents of the input file: {self.input}")

        #Computing dz and Nz
        n_max = np.sqrt(np.array(self.input['magnetic permeability'])*np.array(self.input['electric permittivity'])).max(axis=0)
        lambda_min = c_0/(self.input['simulation parameters'][0]*n_max)
        #Computing cell size based on smallest wavelength
        delta_lambda = lambda_min/25
        d_min = np.array(self.input['layer size']).min()
        #Compute cell size based on the smallest layer size
        delta_size = d_min/25 #Assume the smallest dimension is resolved by 10 cells
        #Compute the final cell size by getting the smaller one from the two (d_critical is ignored in this case) and getting half of it to make sure it is within acceptable
        self.sim_param['dz'] = min(delta_lambda,delta_size)/2

        #Compute the equivalent number of cells in the simulation model
        model_ncells = []
        for layer in self.input['layer size']:
            model_ncells.append(int(np.ceil(layer/self.sim_param['dz'])))
        model_cells = sum(model_ncells)
        self.sim_param['Nz'] = model_cells
        #Check if there is a specified spacer regions
        if spacers != 0: #If there is a specified spacer region
            self.sim_param['Nz'] += 2*spacers
            self.spacer_cells = spacers
        else: #If there is not, add half of the simulation model's size to both end of the computational domain.
            self.sim_param['Nz'] += model_cells
            self.spacer_cells = int(np.ceil(model_cells/2))
        
        #Check if it is divisible by two
        while self.sim_param['Nz'] % 2 != 0: #If the total cell size is not divisible by 2...
            self.sim_param['Nz'] += 1 #Add 1 cell until it is divisble by 2

        #Check to make sure that the injection point is located inside the spacer region (left side)
        if injection_point == 0:
            injection_point = int(np.ceil(self.spacer_cells/5))
            print(f"Injection point: {injection_point}-th cell")
        try:
            assert injection_point < self.spacer_cells and injection_point > 0
        except AssertionError:
            print("Error: Injection point must be located before the simulation model and within the computational domain")
            return

        print('====================================================================================================================================')
        print(f"Cell size based on wavelength: {delta_lambda} m, Cell size based on dimensions: {delta_size} m")
        print(f"Computed cell size (dz): {self.sim_param['dz']} m")
        print(f"Total number of cells (Nz): {self.sim_param['Nz']}")
        
        #print(f"{model_ncells[0]} ->{model_ncells[1]} ->{model_ncells[2]}  ")
        
        #Create computational domain
        #Every cell that is not in the simulation model is assumed to be vacuum
        self.comp_domain['z'] = np.linspace(0,self.sim_param['Nz']*self.sim_param['dz'],int(self.sim_param['Nz']))
        self.comp_domain['mu'] = np.ones(self.comp_domain['z'].size)
        self.comp_domain['epsilon'] = np.ones(self.comp_domain['z'].size)

        #Check the spacer region to determine where to place the simulation model
        start = self.spacer_cells
        end = self.spacer_cells
        for layer in range(int(self.input['simulation parameters'][2])):
            end += model_ncells[layer]
            self.comp_domain['mu'][start:end] = self.input['magnetic permeability'][layer]
            start = end
        #print(f"{self.comp_domain['mu'][0:self.spacer_cells]}->{self.comp_domain['mu'][self.spacer_cells:self.spacer_cells+model_ncells[0]]} -> {self.comp_domain['mu'][self.spacer_cells+model_ncells[0]:self.spacer_cells+model_ncells[0]+model_ncells[1]]}->{self.comp_domain['mu'][self.spacer_cells+model_ncells[0]+model_ncells[1]:self.spacer_cells+model_ncells[0]+model_ncells[1]+model_ncells[2]]} -> {self.comp_domain['mu'][self.spacer_cells + model_cells:]}")
        #print(f"{self.comp_domain['mu'][0:self.spacer_cells].shape}->{self.comp_domain['mu'][self.spacer_cells:self.spacer_cells+model_ncells[0]].shape} -> {self.comp_domain['mu'][self.spacer_cells+model_ncells[0]:self.spacer_cells+model_ncells[0]+model_ncells[1]].shape}->{self.comp_domain['mu'][self.spacer_cells+model_ncells[0]+model_ncells[1]:self.spacer_cells+model_ncells[0]+model_ncells[1]+model_ncells[2]].shape} -> {self.comp_domain['mu'][self.spacer_cells + model_cells:].shape}")
        start = self.spacer_cells
        end = self.spacer_cells
        for layer in range(int(self.input['simulation parameters'][2])):
            end += model_ncells[layer]
            self.comp_domain['epsilon'][start:end] = self.input['electric permittivity'][layer]
            start = end
        #print(f"{self.comp_domain['epsilon'][0:self.spacer_cells]}->{self.comp_domain['epsilon'][self.spacer_cells:self.spacer_cells+model_ncells[0]]} -> {self.comp_domain['epsilon'][self.spacer_cells+model_ncells[0]:self.spacer_cells+model_ncells[0]+model_ncells[1]]}->{self.comp_domain['epsilon'][self.spacer_cells+model_ncells[0]+model_ncells[1]:self.spacer_cells+model_ncells[0]+model_ncells[1]+model_ncells[2]]} -> {self.comp_domain['epsilon'][self.spacer_cells + model_cells:]}")
        #print(f"{self.comp_domain['epsilon'][0:self.spacer_cells].shape}->{self.comp_domain['epsilon'][self.spacer_cells:self.spacer_cells+model_ncells[0]].shape} -> {self.comp_domain['epsilon'][self.spacer_cells+model_ncells[0]:self.spacer_cells+model_ncells[0]+model_ncells[1]].shape}->{self.comp_domain['epsilon'][self.spacer_cells+model_ncells[0]+model_ncells[1]:self.spacer_cells+model_ncells[0]+model_ncells[1]+model_ncells[2]].shape} -> {self.comp_domain['epsilon'][self.spacer_cells + model_cells:].shape}")
        
        #Creating refractive index vector
        self.comp_domain['n'] = np.sqrt(self.comp_domain['mu']*self.comp_domain['epsilon'])


        print(f"Sizes of the vectors in computational domain: Length (in meters) of computational domain: {self.sim_param['Nz']*self.sim_param['dz']} m")
        print(f"Z vector: {self.comp_domain['z'].shape}, Mu vector: {self.comp_domain['mu'].shape}, Epsilon vector: {self.comp_domain['epsilon'].shape}, Refractive Index vector: {self.comp_domain['n'].shape}")
        #print(f"Refractive Index vector: {self.comp_domain['n']}")
        print('====================================================================================================================================')
        print("Layout of the Computational Domain:")
        print(f"|----spacer",end="")
        for layer in range(int(self.input['simulation parameters'][2])):
            print(f"---model_layer_{layer +1}",end="")
        print("---spacer----|")
        print()
        print(f"|----{self.spacer_cells} cells",end="")
        for layer in range(int(self.input['simulation parameters'][2])):
            print(f"---{model_ncells[layer]} cells",end="")
        print(f"---{self.spacer_cells} cells----|")
        #Computing CFL condition for the time step
        self.sim_param['dt'] = (self.comp_domain['n'][0]*self.sim_param['dz'])/(2*c_0)
        
        print('====================================================================================================================================')
        print(f"Time step (dt): {self.sim_param['dt']} seconds")

        #Recording of fmax and source type
        self.sim_param['fmax'] = self.input['simulation parameters'][0]
        if self.input['simulation parameters'][1] == 0: #If the source type is Gaussian...
            self.sim_param['source_type'] = "gaussian"
            
        elif self.input['simulation parameters'][1] == 0: #If it is sinusoidal...
            self.sim_param['source_type'] = "sinusoidal"

        return self.comp_domain,self.sim_param

    def computeSource(self):
        pass

    def __str__(self):
        return f"Value of x: {self.x}"


class Source:
    
    def __init__(self,sim_param = None):
        """
        Upon initialization, the source type will be checked to initialize the proper variables
        """
        try:
            assert sim_param != None
        except AssertionError:
            print("Error: Simulation parameters from init() is needed before proceeding")
            return False

        self.fmax = sim_param['fmax']
        self.dt = sim_param['dt']
        self.dz = sim_param['dz']
        self.Nz = sim_param['Nz']
        
    def GaussianSource(self):
        pass

    def SinusoidalSource(self):
        pass
        
       


def main():
    sim = Simulation()
    sim.init(input_path='sample_input.csv')
    

if __name__ == "__main__":
    main()