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
        self.sim_param = {"dz":0,
                          "dt":0,
                          "Nz":0,
                          "Nt":0,
                          "fmax":0,
                          "sim_time":0,
                          "source_type":"no_source",
                          "boundary_type":"dirichlet",
                          "excitation_method":"hard_source"
                          }
        #Storage variable for the input file parameters
        self.input = {}
        self.comp_domain = {}
        self.spacer_cells = 0
        self.source = None

    def initialize(self,input_path=None,spacers=0,injection_point=0): 
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
        print(f"nmax: {n_max}\n")
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
        print(f"model_cells: {model_cells}")
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
        print(f"Final Nz: {self.sim_param['Nz']}")
        #Check to make sure that the injection point is located inside the spacer region (left side)
        try:
            assert injection_point < self.spacer_cells and injection_point > 0
        except AssertionError:
            print("Error: Injection point must be located before the simulation model and within the computational domain")
            return
        
        if injection_point == 0:
            self.comp_domain['injection_point'] = int(np.ceil(self.spacer_cells/5))
            print(f"Injection point: {self.comp_domain['injection_point']}-th cell")
        else:
            self.comp_domain['injection_point'] = injection_point

        print('====================================================================================================================================')
        print(f"Cell size based on wavelength: {delta_lambda} m, Cell size based on dimensions: {delta_size} m")
        print(f"Computed cell size (dz): {self.sim_param['dz']} m")
        print(f"Total number of cells (Nz): {self.sim_param['Nz']}")
        
        #print(f"{model_ncells[0]} ->{model_ncells[1]} ->{model_ncells[2]}  ")
        
        #Create computational domain
        #Every cell that is not in the simulation model is assumed to be vacuum
        self.comp_domain['z'] = np.arange(0,self.sim_param['Nz']*self.sim_param['dz'],self.sim_param['dz'])
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
        print(f"z: {self.comp_domain['z']}")
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
        #Computing CFL condition for the time step and vector for the simulation time
        self.sim_param['dt'] = (self.comp_domain['n'][0]*self.sim_param['dz'])/(2*c_0)

        #Check if there is an inputted simulation time
        if self.input['simulation parameters'][3] == 0 or self.input['simulation parameters'][3] is None:
            self.sim_param['sim_time'] = (self.comp_domain['n'][0]*self.sim_param['Nz']*self.sim_param['dz'])/c_0
        else: 
            self.sim_param['sim_time'] = self.input['simulation parameters'][3]
        print('====================================================================================================================================')
        print(f"Time step (dt):  {self.sim_param['dt']} seconds | Total sim time (T): {self.sim_param['sim_time']} seconds")

        #Recording of fmax and source type
        self.sim_param['fmax'] = self.input['simulation parameters'][0]
        if self.input['simulation parameters'][1] == 0: #If the source type is Gaussian...
            self.sim_param['source_type'] = "gaussian"
            
        elif self.input['simulation parameters'][1] == 1: #If it is sinusoidal...
            self.sim_param['source_type'] = "sinusoidal"
        
        #Initialize source object
        self.computeSource()

        return self.comp_domain,self.sim_param

    def computeSource(self):
        self.source = Source(self.sim_param)
        if self.sim_param['source_type'] == 'gaussian':
            #Call GaussianSource() to calculate the needed values to insert a gaussian pulse
            self.source.GaussianSource(show_plot=True,comp_domain=self.comp_domain)
        elif self.sim_param['source_type'] == 'sinusoidal':
            #Call SinusoidalSource() to calculate the needed values to insert a sine wave
            self.source.SinusoidalSource(show_plot=True,comp_domain=self.comp_domain)
        else:
            return None

    def simulate(self,boundary_condition="dirichlet",source_excitation="hard"):
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

        self.source_param = {"tau":0,
                            "t0":0,
                            "Nt":0,
                            "type":"",
                            "fmax":sim_param['fmax'],
                            "dt":sim_param['dt'],
                            "dz":sim_param['dz'],
                            "Nz":sim_param['Nz'],
                            "sim_time":sim_param['sim_time'] #sim_time
                            }
        self.source_output = {"t":0,
                              "Esrc":0,
                              "Hsrc":0
                             }


    def plot(self,size=[20,13],labels=["X-Axis","Y-Axis","Title"],save=False,filename="source_plot"):
        plt.figure()
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
        plt.title(labels[2])
        plt.plot(self.source_output['t'],self.source_output['Esrc'])
        plt.plot(self.source_output['t'],self.source_output['Hsrc'])
        plt.show()
        if save == True:
            plt.savefig(filename + ".jpeg")

    def __initialize(self,type='',show_plot=False,t0_coeff=1,prop_coeff=1,tau_coeff=1,comp_domain=None):
        try:
            assert comp_domain != None
        except AssertionError:
            print("comp_domain dict is missing. ")

        if self.source_param['type'] == 'gaussian':

            #Computing the initial parameters for the source
            self.source_param['tau'] = 0.5/self.source_param['fmax'] #should be bandwidth
            self.source_param['t0'] = t0_coeff*self.source_param['tau']
            t_propagation = (np.amax(comp_domain['n'])*self.source_param['Nz']*self.source_param['dz'])/c_0
            total_time_initial = tau_coeff*self.source_param['tau'] + prop_coeff*t_propagation
            if self.source_param['sim_time'] < total_time_initial:
                total_time = self.source_param['sim_time']
            else:
                total_time = total_time_initial

        elif self.source_param['type'] == 'sinusoidal':

            #Computing the initial parameters for the source
            self.source_param['tau'] = 3/self.source_param['fmax'] #should be bandwidth
            self.source_param['t0'] = t0_coeff*self.source_param['tau']
            t_propagation = (np.amax(comp_domain['n'])*self.source_param['Nz']*self.source_param['dz'])/c_0
            total_time_initial = tau_coeff*self.source_param['tau'] + prop_coeff*t_propagation
            if self.source_param['sim_time'] < total_time_initial:
                total_time = self.source_param['sim_time']
            else:
                total_time = total_time_initial
        else:
            print("Error: initialize function independently access from Source class. Do not use this function, instead, use GaussianSource() or SinusoidalSource()")

        return total_time, t_propagation

    def GaussianSource(self,show_plot=False,t0_coeff=5,prop_coeff=3,tau_coeff=12,comp_domain=None):
        """
        Generates a Gaussian pulse that will be inserted into the computational domain.
        """

        self.source_param['type'] = "gaussian"

        total_time, t_propagation = self.__initialize(t0_coeff=5,
                                                       prop_coeff=3,
                                                       tau_coeff=12,
                                                       comp_domain=comp_domain
                                                     )
        print(f"Total time: {total_time} seconds")
        
        #Compute the total steps needed in the simulation 
        self.source_param['Nt'] = np.ceil(total_time/self.source_param['dt'])
        print(f"Nt_before={self.source_param['Nt']} with dtype={self.source_param['Nt'].dtype}")
        
        #Computing the refractive index at the injection point
        n_src = comp_domain['n'][comp_domain['injection_point']]
        print(f"dt={self.source_param['dt']}      Nt={self.source_param['Nt']}")
        prod = int(np.ceil(self.source_param['dt']*self.source_param['Nt']))
        print(f"prod={prod}")

        #Generating the time vector
        self.source_output['t'] = np.linspace(0,self.source_param['dt']*self.source_param['Nt'],int(np.ceil(self.source_param['Nt'])))
        print(f"t_0: {self.source_param['t0']}, tau: {self.source_param['tau']}")
        
        #Computing for the electric field component of the source
        t_E = (self.source_output['t']-self.source_param['t0'])/self.source_param['tau']
        print(f"t_E: {t_E}\n t: {self.source_output['t']}")
        self.source_output['Esrc'] = np.exp(-np.power(t_E,2))
        print(f"Esrc: {self.source_output['Esrc']}")
        
        #Computing the magnetic field component of the source
        #Computing the adjustment term for the time step in magnetic field (due to staggered nature of FDTD)
        adj_H = (n_src*self.source_param['dz'])/(2*c_0) + (self.source_param['dt']/2)
        t_H = ((self.source_output['t']-self.source_param['t0']) + adj_H)/self.source_param['tau']

        self.source_output['Hsrc'] = -np.exp(-np.power(t_H,2))


        if show_plot == True:
            self.plot(labels=["Time (in seconds)","Value","Gaussian Source"],save=True)

        return self.source_output['t'],self.source_output['Esrc'],self.source_output['Hsrc']

   

    def SinusoidalSource(self,show_plot=False,t0_coeff=.5,prop_coeff=3,tau_coeff=.5,comp_domain=None):
        """
        Generates a Sinusoidal signal enveloped initially with an exponential function that will be inserted into the computational domain.
        """

        self.source_param['type'] = "sinusoidal"

        total_time, t_propagation = self.__initialize(t0_coeff=.5,
                                                       prop_coeff=3,
                                                       tau_coeff=12,
                                                       comp_domain=comp_domain
                                                     )
        
        print(f"T: {total_time} seconds")
        #Compute the total steps needed in the simulation 
        self.source_param['Nt'] = np.ceil(total_time/self.source_param['dt'])
        
        #Computing the refractive index at the injection point
        n_src = comp_domain['n'][comp_domain['injection_point']]

        #Generating the time vector
        self.source_output['t'] = np.linspace(0,self.source_param['dt']*self.source_param['Nt'],int(np.ceil(self.source_param['Nt'])))
        print(f"t shape: {self.source_output['t'].shape} \n t: {self.source_output['t']}")
        #Initialize source vectors
        self.source_output['Esrc'] = np.zeros((self.source_output['t'].shape))
        self.source_output['Hsrc'] = np.zeros((self.source_output['t'].shape))
        print(f"Esrc shape: {self.source_output['Esrc'].shape}, Hsrc shape: {self.source_output['Hsrc'].shape}")
        
        
        
        #Selecting a sub-array inside time vector to apply exponential function
        t_condition = self.source_output['t'] <= self.source_param['t0']
        t_initial = self.source_output['t'][t_condition]
        t_initial_size = self.source_output['t'][t_condition].size

        #Computing input time vector for both source components
        t_E = (self.source_output['t']-self.source_param['t0'])/self.source_param['tau']
        t_E_initial = t_E[:t_initial_size]
        adj_H = (n_src*self.source_param['dz'])/(2*c_0) + (self.source_param['dt']/2)
        t_H = ((self.source_output['t']-self.source_param['t0'])+adj_H)/self.source_param['tau']
        t_H_initial = t_H[:t_initial_size]

        #Compute the source when the sine wave is enveloped
        self.source_output['Esrc'][:t_initial_size] = np.exp(-np.power(t_E_initial,2))*(np.sin(2*np.pi*self.source_param['fmax']*t_initial))
        self.source_output['Hsrc'][:t_initial_size] = -np.exp(-np.power(t_H_initial,2))*(np.sin(2*np.pi*self.source_param['fmax']*t_initial))

        #Compute the source components after the enveloping period
        self.source_output['Esrc'][t_initial_size:] = (np.sin(2*np.pi*self.source_param['fmax']*self.source_output['t'][t_initial_size:]))
        self.source_output['Hsrc'][t_initial_size:] = -(np.sin(2*np.pi*self.source_param['fmax']*self.source_output['t'][t_initial_size:]))
        print(f"Esrc shape: {self.source_output['Esrc'].shape}, Hsrc shape: {self.source_output['Hsrc'].shape}")
        

        if show_plot == True:
            self.plot(["Time (in seconds)","Value","Sinusoidal Source"],save=True,filename="source_plot_s")

        return self.source_output['t'],self.source_output['Esrc'],self.source_output['Hsrc']

   


       


def main():
    sim = Simulation()
    sim.initialize(input_path='sample_input.csv',injection_point=15)
    sim.computeSource()
    print(f"comp_domain: {sim.comp_domain.keys()}")
    for i in sim.input.keys():
        print(f"key: {i} value: {type(sim.input[i])}")
if __name__ == "__main__":
    main()