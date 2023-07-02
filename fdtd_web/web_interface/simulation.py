#Library imports
import csv
from re import S, sub
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import h5py
import os
import datetime
import time
from asyncio import new_event_loop
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from .utility import *



# Class definitions for 1D FDTD

def print_breaker(type=""):
        if type == "major":
            print("=============================================================")
        elif type == "minor":
            print("-------------------------------------------------------------")


class Source:

    c_0 = constants.speed_of_light

    def __init__(self,sim_param=None,comp_dom=None):
        self.source_param = {}
        self.source_output = {}
        self.source_comp_dom = {}

        # Check if the input arguments are properly entered
        try:
            assert sim_param != None
            assert comp_dom != None
        except AssertionError:
            print("Input arguments for creating a Source object not detected!")
            return None

        self.source_param = sim_param
        self.source_comp_dom = comp_dom

    def InitializeSource(self,
                         t0_coeff = 3,
                         prop_coeff = 1.0,
                         tau_coeff = 1.0,
                         nmax = 1):

        # Computing tau depending on the source type...
        if self.source_param['source_type'] == "gaussian":
            self.source_param['tau'] = 0.5/self.source_param['fmax']

        elif self.source_param['source_type'] == "sinusoidal":
            self.source_param['tau'] = 3/self.source_param['fmax']

        elif self.source_param['source_type'] == "rectangular":
            self.source_param['tau'] = 1/self.source_param['fmax']

        elif self.source_param['source_type'] == "modulated sine":
            self.source_param['tau'] = 0.5/self.source_param['fmax']

        else:
            print("ERROR: Source type selected is not valid.")
            return None

        self.source_param['t0'] = t0_coeff*self.source_param['tau']
        self.source_param['t_prop'] = (nmax*self.source_param['Nz']*self.source_param['dz'])/Source.c_0
        self.source_param['sim_time'] = (tau_coeff*self.source_param['tau']) + (prop_coeff*self.source_param['t_prop'])

        print("Finished computing related constants for creating the source excitation!")

    def GaussianSource(self,
                       t0_coeff = 6.0,
                       prop_coeff = 6.0,
                       tau_coeff = 12.0,
                       nmax = 1,
                       nsrc = 1):


        self.InitializeSource(t0_coeff,prop_coeff,tau_coeff,nmax)

        print_breaker("major")
        print(f"Creating the source excitation.. | Type: {self.source_param['source_type']}")
        print(f"t0: {self.source_param['t0']} | tau: {self.source_param['tau']} | T: {self.source_param['sim_time']}")

        # Calculating Nt
        self.source_param['Nt'] = np.ceil(self.source_param['sim_time']/self.source_param['dt'])
        print(f"dt: {self.source_param['dt']} seconds | Nt: {self.source_param['Nt']} iterations")


        self.source_output['t'] = np.linspace(0,self.source_param['Nt']*self.source_param['dt'],int(self.source_param['Nt']))

        # Computing the time vector specific to Electric field
        self.source_output['t_E'] = (self.source_output['t'] - self.source_param['t0'])/self.source_param['tau']

        # Computing the Electric field component of the source
        self.source_output['Esrc'] = np.exp(-(self.source_output['t_E']**2))

        # Computing the time vector specific to magnetic field
        # adj_H is an adjustment to the time vector bec. H is computed during half time steps
        adj_H = (nsrc*self.source_param['dz']/(2*Source.c_0)) + (self.source_param['dt']/2)
        self.source_output['t_H'] = ((self.source_output['t'] - self.source_param['t0'])+adj_H)/self.source_param['tau']
        
        # Computing the Magnetic field component of the source excitation
        self.source_output['Hsrc'] = -np.sqrt(self.source_comp_dom['n'][int(self.source_param['inj_point'])])*np.exp(-(self.source_output['t_H']**2)) 

    def SinusoidalSource(self,
                         t0_coeff = 3.0,
                         prop_coeff = 3.0,
                         tau_coeff = 3.0,
                         nmax = 1,
                         nsrc = 1):


        self.InitializeSource(t0_coeff,prop_coeff,tau_coeff,nmax)

        print_breaker("major")
        print(f"Creating the source excitation.. | Type: {self.source_param['source_type']}")
        print(f"t0: {self.source_param['t0']} | tau: {self.source_param['tau']} | T: {self.source_param['sim_time']}")

        # Calculating Nt
        self.source_param['Nt'] = np.ceil(self.source_param['sim_time']/self.source_param['dt'])
        print(f"dt: {self.source_param['dt']} seconds | Nt: {self.source_param['Nt']} iterations")

        self.source_output['t'] = np.linspace(0,self.source_param['Nt']*self.source_param['dt'],int(self.source_param['Nt']))

        # Creating a sliced arrays of the time vector
        t_condition = self.source_output['t'][self.source_output['t'] < self.source_param['t0']]

        

        # Computing the time vector specific to Electric field
        self.source_output['t_E'] = (self.source_output['t'] - self.source_param['t0'])/self.source_param['tau']
        t_E_initial = self.source_output['t_E'][:len(t_condition)]

        # Computing the Electric field component of the source
        self.source_output['Esrc'][:len(t_E_initial)] = np.exp(-(t_E_initial**2))*(np.sin(2*np.pi*self.source_param['fmax']*t_E_initial))
        self.source_output['Esrc'][len(t_E_initial)+1:] = (np.sin(2*np.pi*self.source_param['fmax']*self.source_output['t'][len(t_E_initial):]))

        # Computing the time vector specific to magnetic field
        # adj_H is an adjustment to the time vector bec. H is computed during half time steps
        adj_H = (nsrc*self.source_param['dz']/(2*Source.c_0)) + (self.source_param['dt']/2)
        self.source_output['t_H'] = ((self.source_output['t'] - self.source_param['t0'])+adj_H)/self.source_param['tau']
        t_H_initial = self.source_output['t_H'][:len(t_condition)]

        # Computing the Magnetic field component of the source excitation
        self.source_output['Hsrc'][:len(t_E_initial)] = -np.sqrt(self.source_comp_dom['n'][int(self.source_param['inj_point'])])*np.exp(-(t_H_initial**2))*(np.sin(2*np.pi*self.source_param['fmax']*t_H_initial))
        self.source_output['Hsrc'] = -np.sqrt(self.source_comp_dom['n'][int(self.source_param['inj_point'])])*(np.sin(2*np.pi*self.source_param['fmax']*self.source_output['t'][len(t_H_initial):]))


    def RectWaveSource(self,
                       delay = 0.0,
                       width = 0.0,
                       nmax = 1,
                       nsrc = 1):


        self.InitializeSource(delay,width,12.0,nmax)

        print_breaker("major")
        print(f"Creating the source excitation.. | Type: {self.source_param['source_type']}")
        print(f"t0: {self.source_param['t0']} | tau: {self.source_param['tau']} | T: {self.source_param['sim_time']}")

        # Calculating Nt
        self.source_param['Nt'] = np.ceil(self.source_param['sim_time']/self.source_param['dt'])
        print(f"dt: {self.source_param['dt']} seconds | Nt: {self.source_param['Nt']} iterations")


        self.source_output['t'] = np.linspace(0,self.source_param['Nt']*self.source_param['dt'],int(self.source_param['Nt']))

        # Getting the start and end index for the rectangular pulse
        start_index = np.floor(self.source_param['tau']/self.source_param['dt'])
        end_index = np.ceil((self.source_param['tau'] + self.source_param['t0']/4)/self.source_param['dt'])

        spacing = np.linspace(0.99,0,10)
        # Computing the Electric field component of the source
        self.source_output['Esrc'][start_index:end_index+1] = 1
        self.source_output['Esrc'][start_index-len(spacing):start_index] = np.flip(spacing)
     
        # Computing the Magnetic field component of the source excitation
        self.source_output['Hsrc'][start_index:end_index+1] = -1
        self.source_output['Esrc'][end_index:end_index+len(spacing)+1] = -1*spacing
     
    def ModulatedSineSource(self,
                       t0_coeff = 3.0,
                       prop_coeff = 3.0,
                       tau_coeff = 3.0,
                       nmax = 1,
                       nsrc = 1):


        self.InitializeSource(t0_coeff,prop_coeff,tau_coeff,nmax)

        print_breaker("major")
        print(f"Creating the source excitation.. | Type: {self.source_param['source_type']}")
        print(f"t0: {self.source_param['t0']} | tau: {self.source_param['tau']} | T: {self.source_param['sim_time']}")

        # Calculating Nt
        self.source_param['Nt'] = np.ceil(self.source_param['sim_time']/self.source_param['dt'])
        print(f"dt: {self.source_param['dt']} seconds | Nt: {self.source_param['Nt']} iterations")


        self.source_output['t'] = np.linspace(0,self.source_param['Nt']*self.source_param['dt'],int(self.source_param['Nt']))

        # Computing the time vector specific to Electric field
        self.source_output['t_E'] = (self.source_output['t'] - self.source_param['t0'])/self.source_param['tau']

        # Computing the Electric field component of the source
        self.source_output['Esrc'] = np.sin(2*np.pi*self.source_param['fmax']*self.source_output['t'])*np.exp(-(self.source_output['t_E']**2))

        # Computing the time vector specific to magnetic field
        # adj_H is an adjustment to the time vector bec. H is computed during half time steps
        adj_H = (nsrc*self.source_param['dz']/(2*Source.c_0)) + (self.source_param['dt']/2)
        self.source_output['t_H'] = ((self.source_output['t'] - self.source_param['t0'])+adj_H)/self.source_param['tau']
        
        # Computing the Magnetic field component of the source excitation
        self.source_output['Hsrc'] = -np.sqrt(self.source_comp_dom['n'][int(self.source_param['inj_point'])])*np.sin(2*np.pi*self.source_param['fmax']*self.source_output['t'])*np.exp(-(self.source_output['t_H']**2))


class Subdomain:

    c_0 = constants.speed_of_light

    def __init__(self,sim_param=None,sources=None,comp_domain=None,id=None):

        try:
            assert sim_param != None
            assert sources != None
            assert comp_domain != None
            assert isinstance(id,int)
        except AssertionError:
            print_breaker("minor")
            print("ERROR: Invalid input arguments detected")
            return None

        
        self.subdom_param = sim_param
        self.source = sources
        self.comp_dom = comp_domain
        self.subdom_param['id'] = id
        self.subdom_fields = {}
        self.E_boundary_terms = [0.0,0.0]
        self.H_boundary_terms = [0.0,0.0]
        self.output = {}
      
        #Initialize ghost cells
        self.l_ghost = 0.0
        self.r_ghost = 0.0

        # Initializing the field vectors to 0
        self.subdom_fields['E'] = np.zeros(self.comp_dom['mu'].shape)
        self.subdom_fields['H'] = np.zeros(self.comp_dom['mu'].shape)
        self.subdom_fields['m_E'] = np.zeros(self.comp_dom['mu'].shape)
        self.subdom_fields['m_H'] = np.zeros(self.comp_dom['mu'].shape) 

        #self.output['E'] = self.subdom_fields['E']
        ##self.output['H'] = self.subdom_fields['H']
        #self.output['m_E'] = self.subdom_fields['m_E']
        #self.output['m_H'] = self.subdom_fields['m_H']
  
        self.output['subdom_time'] = 0.0
        print_breaker('minor')
        print(f"Subdomain {self.subdom_param['id']}:")
        print(f"E shape: {self.subdom_fields['E'].shape} | H shape: {self.subdom_fields['H'].shape}")
        print(f"m_E shape: {self.subdom_fields['m_E'].shape} | m_H shape: {self.subdom_fields['m_H'].shape}")
        print(f"Mu: {self.comp_dom['mu'].shape} | Epsilon: {self.comp_dom['epsilon'].shape}")
        print(f"ID: {id}")
        
        
        if self.subdom_param['id'] == 0:
           
            # Check if the 1st subdomain has left spacer region
            if self.subdom_param['subdomain_size'] < self.subdom_param['inj_point']:
                self.subdom_param['inj_point'] = np.floor(self.subdom_param['subdomain_size']/4)

            # Store the sources into the 1st subdomain
            self.source = sources

            # Adjust the injection point (due to the padding)
            self.subdom_param['inj_point'] += self.subdom_param['overlap']

        start = int(self.subdom_param['overlap'])
        end = int(self.comp_dom['mu'].shape[0] - self.subdom_param['overlap'])
        # Compute the update coefficients
        if self.subdom_param['id'] == 0:
            # For the 1st subdomain (padding on the left)
            self.subdom_fields['m_E'][start:] = (Subdomain.c_0*self.subdom_param['dt'])/(self.comp_dom['epsilon'][start:]*self.subdom_param['dz'])
            self.subdom_fields['m_H'][start:] = (Subdomain.c_0*self.subdom_param['dt'])/(self.comp_dom['mu'][start:]*self.subdom_param['dz'])

        elif self.subdom_param['id'] == self.subdom_param['n_subdom'] - 1:
            # For the last subdomain (padding on the right)
            self.subdom_fields['m_E'][:end] = (Subdomain.c_0*self.subdom_param['dt'])/(self.comp_dom['epsilon'][:end]*self.subdom_param['dz'])
            self.subdom_fields['m_H'][:end] = (Subdomain.c_0*self.subdom_param['dt'])/(self.comp_dom['mu'][:end]*self.subdom_param['dz'])

        else:
            # For subdomains without padding
            self.subdom_fields['m_E'] = (Subdomain.c_0*self.subdom_param['dt'])/(self.comp_dom['epsilon']*self.subdom_param['dz'])
            self.subdom_fields['m_H'] = (Subdomain.c_0*self.subdom_param['dt'])/(self.comp_dom['mu']*self.subdom_param['dz'])

    def simulate_subdom(self,curr_iter = 0, boundary_condition="dirichlet",excitation_method="hard"):
        
        # Before starting the simulation, copy the boundary data 1st. The boundary condition will not be dirichlet for the internal boundary data.
        # "Ghost cells" will be used as a way to calculate the last cells in each subdomains, that way, the FDTD algorithm can completely finish
        # This is done after every iteration in the FDTD Time loop, meaning this will be the "Transfer of boundary data stage" of the Schwarz Alternating Method.

        # Check to make sure that the subdom_fields are not empty
        try:
            assert self.subdom_param != None
            assert self.subdom_fields != None
            assert self.source != None
            assert isinstance(curr_iter,int)

        except AssertionError:
            print_breaker("minor")
            print("ERROR: Invalid input arguments detected")
            return None

        #print(f"Simulating in Subdomain {self.subdom_param['id']}")
        # Store the parameters in the subdom_param dict
        self.subdom_param['boundary_condition'] = boundary_condition
        self.subdom_param['excitation_method'] = excitation_method
        #print(self.subdom_param)
        # Initialize array indices
        start = 0
        stop = 0
        
        if self.subdom_param['id'] == 0:
            
            start = int(self.subdom_param['overlap'])
            end = int(self.comp_dom['mu'].shape[0])

        elif self.subdom_param['id'] == self.subdom_param['n_subdom'] -1 :
            start = 0
            end = int(self.comp_dom['mu'].shape[0] - self.subdom_param['overlap'])

        else:
            start = 0
            end = int(self.comp_dom['mu'].shape[0])


        # Step 1: Store boundary data for the 1st subdomain (for the external boundary data)
        if self.subdom_param['id'] == 0: #For the left EXT boundary of 1st subdomain
            if self.subdom_param['boundary_condition'] == "pabc":
                E_value = self.E_boundary_terms.pop(0)
                self.subdom_fields['E'][0] = E_value
                self.E_boundary_terms.append(self.subdom_fields['E'][1]) 

            elif self.subdom_param['boundary_condition'] == "dirichlet":
                self.subdom_fields['E'][0] = 0.0

        else: # For the left INT boundary of all subdomains except the 1st
            # Use the ghost cells here by updating the leftmost index (0) using the update equation
            #print(f"Subdomain ID: {self.subdom_param['id']}")
            self.subdom_fields['E'][0] = self.subdom_fields['E'][0] + (self.subdom_fields['m_E'][0]*(self.subdom_fields['H'][0] - self.l_ghost))


        # Step 2: Update the H vector from E
        self.subdom_fields['H'][:-1] = self.subdom_fields['H'][:-1] + (self.subdom_fields['m_H'][:-1]*(self.subdom_fields['E'][1:] - self.subdom_fields['E'][:-1]))


        # Step 3: Update source excitation (applicable only when the subdom is the 1st one)
        if self.subdom_param['id'] == 0:
            #print("Updating Hsrc")
            if self.subdom_param['excitation_method'] == "tfsf":
                self.subdom_fields['H'][int(self.subdom_param['inj_point']) -1 ] -=  self.subdom_fields['m_H'][int(self.subdom_param['inj_point']) - 1]*self.source['Hsrc'][curr_iter]

        # Step 4: Store H boundary terms
        if self.subdom_param['id'] == 0: #For the left EXT boundary of last subdomain
            if self.subdom_param['boundary_condition'] == "pabc":
                H_value = self.H_boundary_terms.pop(0)
                self.subdom_fields['H'][-1] = H_value
                self.H_boundary_terms.append(self.subdom_fields['H'][-2]) 

            elif self.subdom_param['boundary_condition'] == "dirichlet":
                self.subdom_fields['H'][-1] = 0.0

        else: # For the right INT boundary of all subdomains except the last
            # Use the ghost cells here by updating the rightmost index (n) using the update equation

            self.subdom_fields['H'][-1] = self.subdom_fields['H'][-1] + (self.subdom_fields['m_H'][-1]*( self.r_ghost - self.subdom_fields['E'][-1]))


        # Step 5: Update E from H
        self.subdom_fields['E'][1:] = self.subdom_fields['E'][1:] + (self.subdom_fields['m_E'][1:]*(self.subdom_fields['H'][1:] - self.subdom_fields['H'][:-1]))

        # Step 6: Update E source excitaiton
        if self.subdom_param['id'] == 0:
            #print("Updating Esrc")
            if self.subdom_param['excitation_method'] == "hard":
                self.subdom_fields['E'][int(self.subdom_param['inj_point'])] = self.source['Esrc'][curr_iter]

            elif self.subdom_param['excitation_method'] == "soft":
                self.subdom_fields['E'][int(self.subdom_param['inj_point'])] += self.source['Esrc'][curr_iter]

            elif self.subdom_param['excitation_method'] == "tfsf":
                self.subdom_fields['E'][int(self.subdom_param['inj_point'])] -= (self.subdom_fields['m_E'][int(self.subdom_param['inj_point'])]*self.source['Esrc'][curr_iter])


        # Step 7: Store the data of the EXT boundaries for FFT computations in 1st and last subdomains


        # Step 8: Saving the current snapshot in time in the output matrices
        if curr_iter == 0:

            # Saving the field values for the first time
            self.output['E'] = self.subdom_fields['E']
            self.output['H'] = self.subdom_fields['H']
            self.output['m_E'] = self.subdom_fields['m_E']
            self.output['m_H'] = self.subdom_fields['m_H']
            #self.output['R'] = self.subdom_fields['R']
            #self.output['T'] = self.subdom_fields['T']
            #self.output['S'] = self.subdom_fields['S']
            #self.output['C'] = self.subdom_fields['C']

        else:
            # If this is not the first time, stack the results vertically to create 2D matrices
            # Each row is attributed to the current iteration in time and each column is each cell in the comp domain
            # 2D matrices shape: (Nt,Nz)

            self.output['E'] = np.vstack((self.output['E'],self.subdom_fields['E']))
            self.output['H'] = np.vstack((self.output['H'],self.subdom_fields['H']))
            #self.output['R'] = np.vstack((self.output['R'],self.subdom_fields['R']))
            #self.output['T'] = np.vstack((self.output['T'],self.subdom_fields['T']))
            #self.output['S'] = np.vstack((self.output['S'],self.subdom_fields['S']))
            #self.output['C'] = np.vstack((self.output['C'],self.subdom_fields['C']))





        return None






class Simulation:

    #Initializing class variables
    N_lambda = 10 #Amount of samples to resolve the smallest wavelength
    N_d = 1 #Amount of samples to resolve the smallest dimension in the computational domain
    c_0 = constants.speed_of_light
    mu_0 = constants.mu_0
    epsilon_0 = constants.epsilon_0
    
    def __init__(self,input_filepath="",json_data=None):
        """
        At object creation, the input file (csv) where the simulation parameters are set should be included.
        """

        # Initializing instance variables...
        self.sim_param = {}
        self.comp_domain = {}
        self.input_data = {}
        self.sim_fields = {}
        self.sim_source = {}
        self.output_data = {}
        self.subdomains = []
        self.date_str = ""
        self.user_data = {}

        if input_filepath == "manual":
            # Guide to indices of simulation_parameter vector
            #    -Index 0: fmax
            #    -Index 1: source type
            #    -Index 2: n_model
            #    -Index 3: t_sim
            print_breaker("major")
            print("Creating a simulation object...")
            initial_input = True
            while(initial_input):
                try:
                    fmax = float(input("Frequency of interest/Bandwidth (Hz): "))
                    source_type = input("Source Type (gaussian, sinusoidal, rectangular, modulated sine): ")
                    assert source_type == "gaussian" or source_type == "sinusoidal" or source_type == "rectangular" or source_type == "modulated sine"

                    n_model = int(input("Number of layers in the device model: "))

                    
                    self.sim_param["fmax"] = fmax
                    self.sim_param["source_type"] = source_type
                    self.sim_param["n_model"] = n_model
                    initial_input = False
                except (ValueError,AssertionError):
                    print("Invalid inputs. Try again...")
            n_count = 0 #Used to count the correct inputted layers
            self.input_data["layer_size"] = np.array([])
            self.input_data["magnetic_permeability"] = np.array([])
            self.input_data["electric_permittivity"] = np.array([])

            while(n_count < self.sim_param["n_model"]):
                print_breaker("minor")
                print(f"Layer {n_count+1}: ")
                try:
                    layer_size = float(input("Layer size (in meters): "))
                    magnetic_permeability = float(input("Relative magnetic permeability: "))
                    electric_permittivity = float(input("Relative electric permittivity: "))

                    self.input_data["layer_size"] = np.append(self.input_data["layer_size"],layer_size)
                    self.input_data["magnetic_permeability"] = np.append(self.input_data["magnetic_permeability"],magnetic_permeability)
                    self.input_data["electric_permittivity"] = np.append(self.input_data["electric_permittivity"],electric_permittivity)
                    n_count += 1

                except ValueError:
                    print("Invalid input. Try again....")
        elif json_data is not None:
            # When parsed json data is passed instead of csv file
            # json_data is a dictionary from the parsed data
            date_format = '%d/%m/%Y @ %H:%M:%S'
            self.date_str = datetime.datetime.strptime(json_data['creation_datetime'],date_format)
            
            self.input_data['layer_size'] = parseSimParam(json_data['layer_size'])
            self.input_data['magnetic_permeability'] = parseSimParam(json_data['mu'])
            self.input_data['electric_permittivity'] = parseSimParam(json_data['epsilon'])

            self.sim_param['fmax'] = float(json_data['fmax'])
            self.sim_param["source_type"] = json_data['source_type']
            self.sim_param["n_model"] = int(json_data['n_model'])

            self.sim_param['algo'] = json_data['algo']
            self.sim_param['boundary_cond'] = json_data['boundary_cond']
            self.sim_param['multithreading_flag'] = json_data['multithreading']
            self.sim_param['source_excitation'] = json_data['source_excitation']
            self.sim_param['output_type'] = json_data['output_type']

            self.user_data['username'] = json_data['username']
            self.user_data['user_email'] = json_data['user_email']
            self.user_data['sim_description'] = json_data['sim_description']
            self.user_data['custom_name'] = json_data['custom_name']
       

        else:
            # Guide to indices of simulation_parameter vector
            #    -Index 0: fmax
            #    -Index 1: source type
            #    -Index 2: n_model
            #    -Index 3: t_sim

            self.input_data['layer_size'] = np.array([])
            self.input_data["magnetic_permeability"] = np.array([])
            self.input_data["electric_permittivity"] = np.array([])
            
            layer_size = []
            magnetic_perm = []
            electric_perm = []
            sim_parameters = []
            print_breaker("major")
            print(f"Loading the csv file {input_filepath} ....")
            #When the filepath is a proper one, read the file using the csv module
            with open(input_filepath,mode='r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    layer_size.append(float(row['layer size']))
                    magnetic_perm.append(float(row['magnetic permeability']))
                    electric_perm.append(float(row['electric permittivity']))
                    sim_parameters.append(float(row['simulation parameters']))

            n_model = int(sim_parameters[2])
            #Based on the received data, store the necessary values...
            self.sim_param['fmax'] = float(sim_parameters[0])
            self.sim_param['n_model'] = n_model
            if sim_parameters[1] == 0:
                self.sim_param['source_type'] = "gaussian"
            elif sim_parameters[1] == 1:
                self.sim_param['source_type'] = "sinusoidal"
            elif sim_parameters[1] == 2:
                self.sim_param['source_type'] = "rectangular"
            elif sim_parameters[1] == 3:
                self.sim_param['source_type'] = "modulated sine"


            for i in range(n_model):
                self.input_data["layer_size"] = np.append(self.input_data["layer_size"],layer_size[i])
                self.input_data['magnetic_permeability'] = np.append(self.input_data['magnetic_permeability'],magnetic_perm[i])
                self.input_data['electric_permittivity'] = np.append(self.input_data['electric_permittivity'],electric_perm[i])


        #Print the loaded data
        print_breaker("major")
        print("Simulation Parameters:")
        print(self.sim_param)
        print("Device Model:")
        print("Layer # \t Layer size \t Rel. magnetic perm. \t Rel. electric perm.")

        for i in range(self.sim_param['n_model']):
            print(f"{i+1} \t {self.input_data['layer_size'][i]} \t {self.input_data['magnetic_permeability'][i]} \t {self.input_data['electric_permittivity'][i]}")

    def init_comp_domain(self,
                        spacer = 0,
                        inj_point = 0,
                        n_subdom = 1,
                        overlap = 0,
                        multithread = False,
                        algo="fdtd"):
        # This function is the "pre-processing" stage of the simulation. All of the needed pre-requisite computations are done before the actual
        # simulation is done. The simulate() can be used after this function successfully finished.
        # Input arguments:
        # 1. spacer_cells - in meters. amount of spacing between end of comp domain and the device model (will be adjusted in the computations)
        # 2. injection_point - Position (in cell index) of where the source is injected in the comp domain. This will always be within the spacer region.
        # 3. num_subdomains - Number of subdomains created in FDTD-Schwarz Algorithm.
        # 4. multithread - flag to determine whether the algorithm used is serial or parallel.
        # 5. overlap - amount of overlap used in Schwarz method. Range (0,1) 0% to 100% of spacer region. If the overlap size is greater than 1, it is number cells.
        # 6. algorithm - toggles between basic fdtd and fdtd-schwarz algorithm.

        # Double check that the sim_param and input_data dicts are not empty....
        print_breaker("major")
        print("Entering initialization of computational domain...")
       

        try:

            assert bool(self.sim_param)
            assert bool(self.input_data)
            if bool(self.user_data):
                #if the program is run from the web interface, use the values in the json file
                # javascript file cant parse subdomain in json data 
                spacer = 0
                inj_point = 0
                n_subdom = self.sim_param['n_subdom']

        except AssertionError:
            print("sim_param and input_data dictionaries are empty. Exiting program...")
            print_breaker("major")
            return None

        #Saving the number of subdomains
        self.sim_param['n_subdom'] = n_subdom

        # Computing for dz...
        n_max = np.amax(np.sqrt(self.input_data['magnetic_permeability']*self.input_data['electric_permittivity']))

        # Computing the cell size based on the smallest wavelength (fmax)
        lambda_min = Simulation.c_0/(n_max*self.sim_param['fmax'])
        delta_lambda = lambda_min/Simulation.N_lambda

        # Computing the cell size based on the smallest layer size (minimum size in meters)
        d_min = np.amin(self.input_data['layer_size'])
        delta_size = d_min/Simulation.N_d

        # Getting the minimum cell size of the two computed variables
        self.sim_param['dz'] = min(delta_lambda,delta_size)

        # Computing Nz and number of cells per layer size
        model_ncells = np.ceil(self.input_data['layer_size']/self.sim_param['dz'])
        self.sim_param['Nz'] = model_ncells.sum()

        print_breaker("minor")
        print("Number of cells before adjustment (if fdtd-schwarz will be used)")
        print("Cell amount per layer: ")
        for i in range(self.sim_param['n_model']):
            if i == self.sim_param['n_model']-1:
                print(f"Layer {i+1}")
            else:
                print(f"Layer {i+1}",end=' \t ')
        for i in range(self.sim_param['n_model']):
            if i == self.sim_param['n_model']-1:
                print(f"{model_ncells[i]}")
            else:
                print(f"{model_ncells[i]}",end=' \t ')
        print(f"Number of cells of the device model (Nz): {self.sim_param['Nz']} cells")
        print(f"Cell size (dz): {self.sim_param['dz']} m")        

        # Checking the spacer cells if it is valid
        if spacer != 0: #If the spacer arg is not 0, use it to create the amount of spacing in between the device model
            self.sim_param['Nz'] += 2*(spacer/self.sim_param['dz'])
            self.sim_param['left_spacer'] = np.ceil(spacer/2)
            self.sim_param['right_spacer'] = np.ceil(spacer/2)
     
        else: #Make spacer regions in the left and right side of the device model using Nz

            if self.sim_param['Nz'] % 2 != 0: 
                # Make sure to not round both spacer region (only the left side) to make Nz consistent with vector lengths
                self.sim_param['left_spacer'] = np.ceil(self.sim_param['Nz']/2)
                self.sim_param['right_spacer'] = np.floor(self.sim_param['Nz']/2)
                self.sim_param['Nz'] += self.sim_param['Nz']
            else:
                # No need to round the values here since Nz is divisible by 2
                self.sim_param['left_spacer'] = self.sim_param['Nz']/2
                self.sim_param['right_spacer'] = self.sim_param['Nz']/2
                self.sim_param['Nz'] += self.sim_param['Nz']

        # Adjust Nz and spacers if algo is fdtd-schwarz
        self.sim_param['algo'] = algo
        if self.sim_param['algo'] == "fdtd-schwarz":
            while self.sim_param['Nz'] % self.sim_param['n_subdom'] != 0: 
                # While Nz is not divisible by n_subdom, adjust the value by incrementing the left spacer region
                self.sim_param['Nz'] += 1
                self.sim_param['left_spacer'] += 1

            print_breaker('minor')
            print("Adjusted cell amounts")
            print(f"Number of cells of the device model (Nz): {self.sim_param['Nz']} cells")
            
        print(f"Left spacer region: {self.sim_param['left_spacer']} cells | Right spacer region: {self.sim_param['right_spacer']} cells")

        # Convert the injection point if necessary
        if inj_point < 0:
            print("ERROR: Invalid injection point detected! Exiting program...")
            return None
        elif inj_point > self.sim_param['left_spacer']:
            print('ERROR: Injection point is inside the device model! Exiting program...')
            return None
        else:
            if inj_point == 0:
                #self.sim_param['inj_point'] = np.floor(self.sim_param['left_spacer']/2)
                self.sim_param['inj_point'] = np.floor(self.sim_param['Nz']/2)
            else:
                self.sim_param['inj_point'] = inj_point

        print(f"Injection point (position index in the computational domain): {self.sim_param['inj_point']} index")
        print(f"Nz: {self.sim_param['Nz']}")

        # Creating comp domain vectors...
        self.comp_domain['z'] = np.arange(0,self.sim_param['Nz']*self.sim_param['dz'],self.sim_param['dz'])
        rows_z = self.comp_domain['z'].shape
        assert rows_z == self.sim_param['Nz']
        self.comp_domain['mu'] = np.ones(self.comp_domain['z'].shape)
        self.comp_domain['epsilon'] = np.ones(self.comp_domain['z'].shape)

        start_index = int(self.sim_param['left_spacer'])
        end_index = int(self.sim_param['left_spacer'])
      
        for i in range(len(model_ncells)):
            print(i)
            end_index += model_ncells[i]
            print(f"END INDEX: {end_index}")
            self.comp_domain['mu'][int(start_index):int(end_index)] = self.input_data['magnetic_permeability'][i]
            self.comp_domain['epsilon'][int(start_index):int(end_index)] = self.input_data['electric_permittivity'][i]

            start_index = end_index

        self.comp_domain['n'] = np.sqrt(self.comp_domain['mu']*self.comp_domain['epsilon'])

        print_breaker("minor")
        print(f"Z vector shape: {self.comp_domain['z'].shape} | Mu vector shape: {self.comp_domain['mu'].shape}")
        print(f"Epsilon vector shape: {self.comp_domain['epsilon'].shape} | N vector shape: {self.comp_domain['n'].shape}")

        print("Layout of Computational Domain:")
        print("|----spacer",end="")
        for i in range(self.sim_param['n_model']):
            print(f"----model_layer_{i+1}", end="")
        print("----spacer----|")

        print(f"|---{self.sim_param['left_spacer']} cells",end="")
        for i in model_ncells:
            print(f"---{i} cells",end="")
        print(f"---{self.sim_param['right_spacer']} cells---|")

        # Computing the CFL condition
        self.sim_param['dt'] = (1*self.sim_param['dz'])/(2*Simulation.c_0)

        self.sim_param['sim_time'] = (self.comp_domain['n'][0]*self.sim_param['Nz']*self.sim_param['dz'])/Simulation.c_0

        print_breaker("minor")
        print(f"Time step (dt): {self.sim_param['dt']} seconds | Sim time: {self.sim_param['sim_time']} seconds")

        # Save the initial dz and dt to the list
        self.sim_param['dz_list'] = [self.sim_param['dz']]
        self.sim_param['dt_list'] = [self.sim_param['dt']]
        

        # Computing the source excitation...

        self.compute_source()

        # At this point, the code needs to check if it is in serial or parallel mode. 
        # Meaning, all of the pre-processing for both serial and parallel versions are done in this method.
        # The only data transferred to the subdomain class are the simulation parameters and the
        # computational domain vectors

        self.sim_param['algo'] = algo
        self.sim_param['multithread'] = multithread

        print_breaker('minor')
        print(f"Algorithm: {self.sim_param['algo']} | Multithreading enabled: {self.sim_param['multithread']}")

        if self.sim_param['algo'] == "fdtd":

            #Initialize the fields
            self.sim_fields['E'] = np.zeros(self.comp_domain['z'].shape)
            self.sim_fields['H'] = np.zeros(self.comp_domain['z'].shape)
            # self.sim_fields['m_E'] = np.zeros(self.comp_domain['z'].shape)
            # self.sim_fields['m_H'] = np.zeros(self.comp_domain['z'].shape)

            # Compute the update coefficients
            self.sim_fields['m_E'] = (Simulation.c_0*self.sim_param['dt'])/(self.comp_domain['epsilon']*self.sim_param['dz'])
            self.sim_fields['m_H'] = (Simulation.c_0*self.sim_param['dt'])/(self.comp_domain['mu']*self.sim_param['dz'])

            print_breaker('minor')
            print("Field shapes:")
            print(f"Z shape: {self.comp_domain['z'].shape}")
            print(f"E shape: {self.sim_fields['E'].shape} | m_E shape: {self.sim_fields['m_E'].shape}")
            print(f"H shape: {self.sim_fields['H'].shape} | m_H shape: {self.sim_fields['m_H'].shape}")
        
        elif self.sim_param['algo'] == "fdtd-schwarz":
            
            # Check if the number of subdomains are valid
            print_breaker('minor')
            print(f"Number of subdomains: {self.sim_param['n_subdom']}")

            try:
                assert self.sim_param['n_subdom'] % 2 == 0
                assert self.sim_param['n_subdom'] <= 64 and self.sim_param['n_subdom'] >=1
            except AssertionError:
                print_breaker('minor')
                print("ERROR: Valid number of subdomains not detected.")
                return None

            # Layout of a subdomain:
            # ||--overlapping region--|--non-overlapping region--|--overlapping region--||

            # Computing the size of overlapping region
            if overlap > 0 and overlap < 1: # Assume this is percentage of left spacer
                self.sim_param['overlap'] = self.sim_param['left_spacer']*overlap

            elif overlap > 1: # If it is the number of cells
                self.sim_param['overlap'] = overlap

                #Adjust the overlap size to make sure that Nz + overlap is divisible by the number of subdomains 
                while (self.sim_param['Nz'] + self.sim_param['overlap']) % self.sim_param['n_subdom'] != 0:
                    self.sim_param['overlap'] += 1

            else:
                print_breaker('minor')
                print("ERROR: Valid algorithm not detected")
                return None

            # Computing the non overlapping region size
            self.sim_param['non_overlap'] = int(self.sim_param['Nz'] - ((self.sim_param['n_subdom']-1)*self.sim_param['overlap']))/self.sim_param['n_subdom']

            # The method used here is similar to how a manual STFT is done
            # frame_size = subdomain size
            # hop_size = how much the 'frame' moves in the computational domain.
            self.sim_param['overlap'] = int(self.sim_param['overlap'])
            frame_size = self.sim_param['non_overlap'] + 2*self.sim_param['overlap']
            self.sim_param['subdomain_size'] = int(frame_size)
            start = int(0)
            stop = int(frame_size)
            hop_size = self.sim_param['non_overlap'] + self.sim_param['overlap']

            print_breaker('minor')
            print(f"Non-overlap size: {self.sim_param['non_overlap']} cells | Overlap size: {self.sim_param['overlap']} cells")
            print(f"Frame size: {frame_size} cells | Hop size: {hop_size} cells")

            padded_mu = np.pad(self.comp_domain['mu'],self.sim_param['overlap'],'constant')
            padded_epsilon = np.pad(self.comp_domain['epsilon'],self.sim_param['overlap'],'constant')
            padded_z = np.pad(self.comp_domain['z'],self.sim_param['overlap'],'constant')

            print(f"Padded mu vector: {padded_mu.shape}")
            print(f"Padded epsilon vector: {padded_epsilon.shape}")
            print(f"Padded z vector: {padded_z.shape}")

            # Getting the subsets for each subdomain...
            
            for i in range(self.sim_param['n_subdom']):
                print_breaker('minor')
                print(f"Start: {start} | Stop: {stop}")

                if i == 0:
                    mu_2D = padded_mu[start:stop]
                    epsilon_2D = padded_epsilon[start:stop]
                    z_2D = padded_z[start:stop]
                else:
                    mu_2D = np.vstack((mu_2D,padded_mu[start:stop]))
                    epsilon_2D = np.vstack((epsilon_2D,padded_epsilon[start:stop]))
                    z_2D  = np.vstack((z_2D,padded_z[start:stop]))

                start += int(hop_size)
                stop += int(hop_size)

            print(f"Stacked MU: {mu_2D.shape}")
            print(f"Stacked EPSILON: {epsilon_2D.shape}")

            # At this point, we can just assume that the source will always be injected into the 1st subdomain (center position)

            # Call the preprocess subdomain to create the subdomain objects.
            self.init_subdomains(mu_2D,epsilon_2D,z_2D)

        # Computing df - Frequency steps in Frequency response
        self.sim_param['df'] = 1/(self.sim_param['dt']*self.sim_param['Nt'])

        # Initialize the frequency vectors for computing the freq response
        self.sim_param['n_freq'] = self.sim_param['Nt']
        self.sim_param['fmax_fft'] = 0.5*(1/self.sim_param['dt'])
        #self.output_data['Freq_range'] = np.linspace(0,self.sim_param['fmax_fft'],int(self.sim_param['n_freq']))
        self.output_data['Freq_range'] = np.linspace(0,self.sim_param['fmax'],int(self.sim_param['Nt']))

        print_breaker('minor')
        print(f"Number of freq. samples: {self.sim_param['n_freq']} | Freq vector shape: {self.output_data['Freq_range'].shape}")
        print(f"df: {self.sim_param['df']}")
        # Creating the kernel frequency for FFT computation
        self.sim_fields['Kernel_Freq'] = np.exp(-1j*2.0*np.pi*self.sim_param['dt']*self.output_data['Freq_range'])
        print(f"Kernel freq shape: {self.sim_fields['Kernel_Freq'].shape}")

        # Initialize the FFT vectors
        self.sim_fields['R'] = np.zeros(self.sim_fields['Kernel_Freq'].shape)
        self.sim_fields['T'] = np.zeros(self.sim_fields['Kernel_Freq'].shape)
        self.sim_fields['C'] = np.zeros(self.sim_fields['Kernel_Freq'].shape)
        self.sim_fields['S'] = np.zeros(self.sim_fields['Kernel_Freq'].shape)

        print("FFT Vector shapes:")
        print(f"Reflectance: {self.sim_fields['R'].shape} | Transmittance: {self.sim_fields['T'].shape}")
        print(f"Conservation of Energy: {self.sim_fields['C'].shape} | Source FFT: {self.sim_fields['S'].shape}")

        self.sim_param['preprocess_finished'] = True

    def update_sim_param(self,n_wavelength=0, n_dim = 0):
        print_breaker("major")
        print("Changing the simulation parameters")
        print(f"Current values: N_wavelength: {Simulation.N_lambda + n_wavelength} | N_dim: {Simulation.N_d + n_dim}")

       

        try:

            assert bool(self.sim_param)
            assert bool(self.input_data)
        except AssertionError:
            print("sim_param and input_data dictionaries are empty. Exiting program...")
            print_breaker("major")
            return None

        
        # Computing for dz...
        n_max = np.amax(np.sqrt(self.input_data['magnetic_permeability']*self.input_data['electric_permittivity']))

        # Computing the cell size based on the smallest wavelength (fmax)
        lambda_min = Simulation.c_0/(n_max*self.sim_param['fmax'])
        delta_lambda = lambda_min/(Simulation.N_lambda + n_wavelength)

        # Computing the cell size based on the smallest layer size (minimum size in meters)
        d_min = np.amin(self.input_data['layer_size'])
        delta_size = d_min/(Simulation.N_d + n_dim)

        # The final cell size is obtained by getting the smallest of delta_lambda 
        # and delta_size to make sure that the comp domain can resolve all the 
        # necessary features (wavelength or dimension)
            

        # Getting the minimum cell size of the two computed variables
        self.sim_param['dz'] = min(delta_lambda,delta_size)

        # Computing Nz and number of cells per layer size
        model_ncells = np.ceil(self.input_data['layer_size']/self.sim_param['dz'])
        self.sim_param['Nz'] = model_ncells.sum()

        # Check Nz to have more than 1 number of cells. Otherwise, end the function 
        # and continue to the next iteration..
        try:

            assert self.sim_param['Nz'] > 1

        except AssertionError:

            print(f"Nz value: {self.sim_param['Nz']} is not valid")
            return None

        print_breaker("minor")
        print("Number of cells before adjustment (if fdtd-schwarz will be used)")
        print("Cell amount per layer: ")
        for i in range(self.sim_param['n_model']):
            if i == self.sim_param['n_model']-1:
                print(f"Layer {i+1}")
            else:
                print(f"Layer {i+1}",end=' \t ')
        for i in range(self.sim_param['n_model']):
            if i == self.sim_param['n_model']-1:
                print(f"{model_ncells[i]}")
            else:
                print(f"{model_ncells[i]}",end=' \t ')
        print(f"Number of cells of the device model (Nz): {self.sim_param['Nz']} cells")
        print(f"Cell size (dz): {self.sim_param['dz']} m")        

        # Make the spacer region half of the comp domain (each side)
        self.sim_param['left_spacer'] = int(np.ceil(self.sim_param['Nz']/2))
        self.sim_param['right_spacer'] = int(np.ceil(self.sim_param['Nz']/2))
        self.sim_param['Nz'] += self.sim_param['Nz']


        # Adjust Nz and spacers if algo is fdtd-schwarz
        if self.sim_param['algo'] == "fdtd-schwarz":
            while self.sim_param['Nz'] % self.sim_param['n_subdom'] != 0: 
                # While Nz is not divisible by n_subdom, adjust the value by incrementing the left spacer region
                self.sim_param['Nz'] += 1
                self.sim_param['left_spacer'] += 1

            print_breaker('minor')
            print("Adjusted cell amounts")
            print(f"Number of cells of the device model (Nz): {self.sim_param['Nz']} cells")
            
        print(f"Left spacer region: {self.sim_param['left_spacer']} cells | Right spacer region: {self.sim_param['right_spacer']} cells")
        
        while self.sim_param['inj_point'] > (self.sim_param['left_spacer']/2):
            self.sim_param['inj_point'] -= 1

        # Convert the injection point if necessary
        if self.sim_param['inj_point'] < 0:
            print("ERROR: Invalid injection point detected! Exiting program...")
            return None
        elif self.sim_param['inj_point'] > self.sim_param['left_spacer']:
            print('ERROR: Injection point is inside the device model! Exiting program...')
            return None
        else:
         
            self.sim_param['inj_point'] = np.floor(self.sim_param['left_spacer']/2)
           

        print(f"Injection point (position index in the computational domain): {self.sim_param['inj_point']} index")
        print(f"Nz: {self.sim_param['Nz']}")

        # Creating comp domain vectors...
        self.comp_domain['z'] = np.arange(0,self.sim_param['Nz']*self.sim_param['dz'],self.sim_param['dz'])
        rows_z = self.comp_domain['z'].shape
        assert rows_z == self.sim_param['Nz']
        self.comp_domain['mu'] = np.ones(self.comp_domain['z'].shape)
        self.comp_domain['epsilon'] = np.ones(self.comp_domain['z'].shape)

        start_index = int(self.sim_param['left_spacer'])
        end_index = int(self.sim_param['left_spacer'])
      
        for i in range(len(model_ncells)):
            print(i)
            end_index += model_ncells[i]
            print(f"END INDEX: {end_index}")
            self.comp_domain['mu'][int(start_index):int(end_index)] = self.input_data['magnetic_permeability'][i]
            self.comp_domain['epsilon'][int(start_index):int(end_index)] = self.input_data['electric_permittivity'][i]

            start_index = end_index

        self.comp_domain['n'] = np.sqrt(self.comp_domain['mu']*self.comp_domain['epsilon'])

        print_breaker("minor")
        print(f"Z vector shape: {self.comp_domain['z'].shape} | Mu vector shape: {self.comp_domain['mu'].shape}")
        print(f"Epsilon vector shape: {self.comp_domain['epsilon'].shape} | N vector shape: {self.comp_domain['n'].shape}")

        print("Layout of Computational Domain:")
        print("|----spacer",end="")
        for i in range(self.sim_param['n_model']):
            print(f"----model_layer_{i+1}", end="")
        print("----spacer----|")

        print(f"|---{self.sim_param['left_spacer']} cells",end="")
        for i in model_ncells:
            print(f"---{i} cells",end="")
        print(f"---{self.sim_param['right_spacer']} cells---|")

        # Computing the CFL condition
        self.sim_param['dt'] = (1*self.sim_param['dz'])/(2*Simulation.c_0)

        self.sim_param['sim_time'] = (self.comp_domain['n'][0]*self.sim_param['Nz']*self.sim_param['dz'])/Simulation.c_0

        print_breaker("minor")
        print(f"Time step (dt): {self.sim_param['dt']} seconds | Sim time: {self.sim_param['sim_time']} seconds")

        # Save the computed dz and dt
        self.sim_param['dz_list'].append(self.sim_param['dz'])
        self.sim_param['dt_list'].append(self.sim_param['dt'])

        # Computing the source excitation...

        self.compute_source()

        # At this point, the code needs to check if it is in serial or parallel mode. 
        # Meaning, all of the pre-processing for both serial and parallel versions are done in this method.
        # The only data transferred to the subdomain class are the simulation parameters and the
        # computational domain vectors

        print_breaker('minor')
        print(f"Algorithm: {self.sim_param['algo']} | Multithreading enabled: {self.sim_param['multithread']}")

        if self.sim_param['algo'] == "fdtd":

            #Initialize the fields
            self.sim_fields['E'] = np.zeros(self.comp_domain['z'].shape)
            self.sim_fields['H'] = np.zeros(self.comp_domain['z'].shape)
            # self.sim_fields['m_E'] = np.zeros(self.comp_domain['z'].shape)
            # self.sim_fields['m_H'] = np.zeros(self.comp_domain['z'].shape)

            # Compute the update coefficients
            self.sim_fields['m_E'] = (Simulation.c_0*self.sim_param['dt'])/(self.comp_domain['epsilon']*self.sim_param['dz'])
            self.sim_fields['m_H'] = (Simulation.c_0*self.sim_param['dt'])/(self.comp_domain['mu']*self.sim_param['dz'])

            print_breaker('minor')
            print("Field shapes:")
            print(f"Z shape: {self.comp_domain['z'].shape}")
            print(f"E shape: {self.sim_fields['E'].shape} | m_E shape: {self.sim_fields['m_E'].shape}")
            print(f"H shape: {self.sim_fields['H'].shape} | m_H shape: {self.sim_fields['m_H'].shape}")
        
        elif self.sim_param['algo'] == "fdtd-schwarz":
            
            # Check if the number of subdomains are valid
            print_breaker('minor')
            print(f"Number of subdomains: {self.sim_param['n_subdom']}")

            try:
                assert self.sim_param['n_subdom'] % 2 == 0
                assert self.sim_param['n_subdom'] <= 64 and self.sim_param['n_subdom'] >=1
            except AssertionError:
                print_breaker('minor')
                print("ERROR: Valid number of subdomains not detected.")
                return None

            # Layout of a subdomain:
            # ||--overlapping region--|--non-overlapping region--|--overlapping region--||

            # Computing the size of overlapping region
            #Adjust the overlap size to make sure that Nz + overlap is divisible by the number of subdomains 
            while (self.sim_param['Nz'] + self.sim_param['overlap']) % self.sim_param['n_subdom'] != 0:
                self.sim_param['overlap'] += 1

            

            # Computing the non overlapping region size
            self.sim_param['non_overlap'] = int(self.sim_param['Nz'] - ((self.sim_param['n_subdom']-1)*self.sim_param['overlap']))/self.sim_param['n_subdom']

            # The method used here is similar to how a manual STFT is done
            # frame_size = subdomain size
            # hop_size = how much the 'frame' moves in the computational domain.
            self.sim_param['overlap'] = int(self.sim_param['overlap'])
            frame_size = self.sim_param['non_overlap'] + 2*self.sim_param['overlap']
            self.sim_param['subdomain_size'] = int(frame_size)
            start = int(0)
            stop = int(frame_size)
            hop_size = self.sim_param['non_overlap'] + self.sim_param['overlap']

            print_breaker('minor')
            print(f"Non-overlap size: {self.sim_param['non_overlap']} cells | Overlap size: {self.sim_param['overlap']} cells")
            print(f"Frame size: {frame_size} cells | Hop size: {hop_size} cells")

            padded_mu = np.pad(self.comp_domain['mu'],self.sim_param['overlap'],'constant')
            padded_epsilon = np.pad(self.comp_domain['epsilon'],self.sim_param['overlap'],'constant')
            padded_z = np.pad(self.comp_domain['z'],self.sim_param['overlap'],'constant')

            print(f"Padded mu vector: {padded_mu.shape}")
            print(f"Padded epsilon vector: {padded_epsilon.shape}")
            print(f"Padded z shape: {padded_z.shape}")

            # Getting the subsets for each subdomain...
            mu_2D = 0
            epsilon_2D = 0
            z_2D = 0
            for i in range(self.sim_param['n_subdom']):
                print_breaker('minor')
                print(f"Start: {start} | Stop: {stop}")

                if i == 0:
                    mu_2D = padded_mu[start:stop]
                    epsilon_2D = padded_epsilon[start:stop]
                    z_2D = padded_z[start:stop]
                else:
                    mu_2D = np.vstack((mu_2D,padded_mu[start:stop]))
                    epsilon_2D = np.vstack((epsilon_2D,padded_epsilon[start:stop]))
                    z_2D = np.vstack((z_2D,padded_z[start:stop]))

                start += int(hop_size)
                stop += int(hop_size)

            print(f"Stacked MU: {mu_2D.shape}")
            print(f"Stacked EPSILON: {epsilon_2D.shape}")

            # At this point, we can just assume that the source will always be injected into the 1st subdomain (center position)

            # Call the preprocess subdomain to create the subdomain objects.
            #   Delete the previous subdomains
            self.subdomains = []
            self.init_subdomains(mu_2D,epsilon_2D,z_2D)

        # This part is needed by both fdtd and fdtd-schwarz algo

        # Computing df - Frequency steps in Frequency response
        self.sim_param['df'] = 1/(self.sim_param['dt']*self.sim_param['Nt'])

        # Initialize the frequency vectors for computing the freq response
        self.sim_param['n_freq'] = self.sim_param['Nt']
        self.sim_param['fmax_fft'] = 0.5*(1/self.sim_param['dt'])
        self.output_data['Freq_range'] = np.linspace(0,self.sim_param['fmax_fft'],int(self.sim_param['n_freq']))
        
        print_breaker('minor')
        print(f"Number of freq. samples: {self.sim_param['n_freq']} | Freq vector shape: {self.output_data['Freq_range'].shape}")
        print(f"df: {self.sim_param['df']}")
        # Creating the kernel frequency for FFT computation
        self.sim_fields['Kernel_Freq'] = np.exp(-1j*2.0*np.pi*self.sim_param['dt']*self.output_data['Freq_range'])
        print(f"Kernel freq shape: {self.sim_fields['Kernel_Freq'].shape}")

        # Initialize the FFT vectors
        self.sim_fields['R'] = np.zeros(self.sim_fields['Kernel_Freq'].shape)
        self.sim_fields['T'] = np.zeros(self.sim_fields['Kernel_Freq'].shape)
        self.sim_fields['C'] = np.zeros(self.sim_fields['Kernel_Freq'].shape)
        self.sim_fields['S'] = np.zeros(self.sim_fields['Kernel_Freq'].shape)

        print("FFT Vector shapes:")
        print(f"Reflectance: {self.sim_fields['R'].shape} | Transmittance: {self.sim_fields['T'].shape}")
        print(f"Conservation of Energy: {self.sim_fields['C'].shape} | Source FFT: {self.sim_fields['S'].shape}")

        self.sim_param['preprocess_finished'] = True

        return True




    def compute_source(self):
        
        # Creating a new Source object
        self.sim_source['source'] = Source(self.sim_param,self.comp_domain)

        print_breaker("major")
        print(f"Computing source excitation. | Source type: {self.sim_param['source_type']}")

        if self.sim_param['source_type'] == "gaussian":
            self.sim_source['source'].GaussianSource(3,12,2,np.amax(self.comp_domain['n']),self.comp_domain['n'][int(self.sim_param['inj_point'])])

        elif self.sim_param['source_type'] == "sinusoidal":
            self.sim_source['source'].SinusoidalSource(3,12,2,np.amax(self.comp_domain['n']),self.comp_domain['n'][int(self.sim_param['inj_point'])])
            

        elif self.sim_param['source_type'] == "rectangular":
            self.sim_source['source'].RectWaveSource(6,1,np.amax(self.comp_domain['n']),self.comp_domain['n'][int(self.sim_param['inj_point'])])
            

        elif self.sim_param['source_type'] == "modulated sine":
            self.sim_source['source'].ModulatedSineSource(3,12,2,np.amax(self.comp_domain['n']),self.comp_domain['n'][int(self.sim_param['inj_point'])])
            

        # Transfer computed parameters to the Simulation object
        self.sim_param['Nt'] = self.sim_source['source'].source_param['Nt']
        self.sim_param['tau'] = self.sim_source['source'].source_param['tau']
        self.sim_param['t0'] = self.sim_source['source'].source_param['t0']

        self.sim_source['t'] = self.sim_source['source'].source_output['t']
        self.sim_source['t_E'] = self.sim_source['source'].source_output['t_E']
        self.sim_source['t_H'] = self.sim_source['source'].source_output['t_H']
        self.sim_source['Esrc'] = self.sim_source['source'].source_output['Esrc']
        self.sim_source['Hsrc'] = self.sim_source['source'].source_output['Hsrc']

        # Save the computed source to the output data dict
        self.output_data['t'] = self.sim_source['t']
        self.output_data['t_E'] = self.sim_source['t_E']
        self.output_data['t_H'] = self.sim_source['t_H']
        self.output_data['Esrc'] = self.sim_source['Esrc']
        self.output_data['Hsrc'] = self.sim_source['Hsrc']

        print_breaker('minor')
        print(f"Esrc shape: {self.sim_source['Esrc'].shape} | Hsrc shape: {self.sim_source['Hsrc'].shape}")


    def init_subdomains(self,mu_2D,epsilon_2D,z_2D):
        

        # Create subdomains based on the specified number of subdomains
        for i in range(int(self.sim_param['n_subdom'])):
            # Create a new dict for the subdomain
            subdomain = {}
            subdomain['mu'] = mu_2D[i,:]
            subdomain['epsilon'] = epsilon_2D[i,:]
            subdomain['z'] = z_2D[i,:]
            #print(f"Shapes: Mu: {subdomain['mu'].shape} | Epsilon: {subdomain['epsilon'].shape}")
            
            self.subdomains.append(Subdomain( self.sim_param.copy(),
                                              self.sim_source,
                                              subdomain,
                                              i
                                                ))
           # print(f"Initializing: Subdomain {self.subdomains[i].subdom_param['id']}")
            #print_breaker('minor')
            #print(f"Subdomain {i}:")
            #print(f"Mu shape: {self.subdomains[i].comp_dom['mu'].shape} | Epsilon shape: {self.subdomains[i].comp_dom['epsilon'].shape} ")

    def simulate(self,boundary_condition = "dirichlet", excitation_method = "hard"):
        
        # Save the input arguments
        try:
            assert boundary_condition == "dirichlet" or boundary_condition == "pabc"
            assert excitation_method == "hard" or excitation_method == "soft" or excitation_method == "tfsf"
        except AssertionError:
            print_breaker('minor')
            print(f"ERROR: Invalid input arguments detected")

        # Check if the comp_domain is processed properly
        try:
            assert self.sim_param['preprocess_finished']
        except AssertionError:
            print_breaker('minor')
            print(f"Preprocessing is not yet complete/finished incomplete. The simulation will not proceed.")
            return None

        self.sim_param['boundary_condition'] = boundary_condition
        self.sim_param['excitation_method'] = excitation_method

        # Check the algorithm selected
        if self.sim_param['algo'] == "fdtd":
            self.simulate_serial()
        elif self.sim_param['algo'] == "fdtd-schwarz":
            self.simulate_fdtd_schwarz()


    def simulate_serial(self):
        
        fdtd_start_time = time.perf_counter()

        print_breaker("major")
        print("Starting the simulation. Current simulation parameters: ")
        print(f"algorithm: {self.sim_param['algo']} | boundary condition: {self.sim_param['boundary_condition']}")
        print(f"source excitation method: {self.sim_param['excitation_method']} | source type: {self.sim_param['source_type']}")
        
        # Initializing temporary variables
        E_boundary_terms = [0.0,0.0]
        H_boundary_terms = [0.0,0.0]
        end_index = int(self.sim_param['Nz'])
        R = np.zeros(self.sim_fields['Kernel_Freq'].shape)
        T = np.zeros(self.sim_fields['Kernel_Freq'].shape)

      
        # FDTD Time Loop
        for curr_iteration in range(int(self.sim_param['Nt'])):

            # Step 1: Store E boundary terms/Set the boundary conditions
            if self.sim_param['boundary_condition'] == "pabc":

                # Remove the value from the list
                E_value = E_boundary_terms.pop(0)

                #Assign that value to the 1st element of the E field vector
                self.sim_fields['E'][0] = E_value

                # Add the 2nd element (index 1) to the list
                E_boundary_terms.append(self.sim_fields['E'][1])

            elif self.sim_param['boundary_condition'] == "dirichlet":

                # Set the 1st element to 0
                self.sim_fields['E'][0] = 0

            # Step 2: Update the H-field from the E-field (FDTD Space Loop)
            self.sim_fields['H'][:-1] = self.sim_fields['H'][:end_index-1] + \
                                        (self.sim_fields['m_H'][:-1])*(self.sim_fields['E'][1:]-self.sim_fields['E'][:-1])

            # Step 3: Insert the H-field source (applicable only to TFSF excitation method)
            if self.sim_param['excitation_method'] == "tfsf":
                self.sim_fields['H'][int(self.sim_param['inj_point'])-1] -= (self.sim_fields['m_H'][int(self.sim_param['inj_point'])-1])*\
                                                                            (self.sim_source['Esrc'][curr_iteration])

            # Step 4: Store H boundary term/set the boundary conditions
            if self.sim_param['boundary_condition'] == "pabc":
                
                # Remove the value from the list
                H_value = H_boundary_terms.pop(0)

                # Assign the value to the last element of the H field vector
                self.sim_fields['H'][-1] = H_value

                # Add the 2nd to the last element to the list
                H_boundary_terms.append(self.sim_fields['H'][-2])
            
            elif self.sim_param['boundary_condition'] == "dirichlet":

                self.sim_fields['H'][-1] = 0

            # Step 5: Update E field from H field (FDTD Space Loop)
            self.sim_fields['E'][1:] = self.sim_fields['E'][1:] + \
                                        self.sim_fields['m_E'][1:]*(self.sim_fields['H'][1:] - self.sim_fields['H'][:-1])

            # Step 6: Inject the E field component of the source 
            if self.sim_param['excitation_method'] == "hard":
                self.sim_fields['E'][int(self.sim_param['inj_point'])] = self.sim_source['Esrc'][curr_iteration]
            elif self.sim_param['excitation_method'] == "soft":
                self.sim_fields['E'][int(self.sim_param['inj_point'])] += self.sim_source['Esrc'][curr_iteration]
            elif self.sim_param['excitation_method'] == "tfsf":
                self.sim_fields['E'][int(self.sim_param['inj_point'])] -= (self.sim_fields['m_E'][int(self.sim_param['inj_point'])]*\
                                                                            self.sim_source['Hsrc'][curr_iteration])


            # Step 7: Compute the Fourier Transform of the current simulation iteration

            # Computing Fourier transform of Reflectance and Transmittance (unnormalized)
            R = R + ((self.sim_fields['Kernel_Freq']**curr_iteration)*self.sim_fields['E'][1])
            T = T + ((self.sim_fields['Kernel_Freq']**curr_iteration)*self.sim_fields['E'][-2])
            self.sim_fields['S'] = self.sim_fields['S'] + ((self.sim_fields['Kernel_Freq']**curr_iteration)*self.sim_source['Esrc'][curr_iteration])

            # Normalizing the computed frequency respones
            #self.sim_fields['R'] = (R/self.sim_fields['S']).real ** 2
            #self.sim_fields['T'] = (T/self.sim_fields['S']).real ** 2
            self.sim_fields['R'] = np.abs(R/self.sim_fields['S']) ** 2
            self.sim_fields['T'] = np.abs(T/self.sim_fields['S']) ** 2
            self.sim_fields['C'] = self.sim_fields['R'] + self.sim_fields['T']

            # Step 8: Saving the current snapshot in time in the output matrices
            if curr_iteration == 0:

                # Saving the field values for the first time
                self.output_data['E'] = self.sim_fields['E']
                self.output_data['H'] = self.sim_fields['H']
                self.output_data['m_E'] = self.sim_fields['m_E']
                self.output_data['m_H'] = self.sim_fields['m_H']
                self.output_data['R'] = self.sim_fields['R']
                self.output_data['T'] = self.sim_fields['T']
                self.output_data['S'] = self.sim_fields['S']
                self.output_data['C'] = self.sim_fields['C']

            else:
                # If this is not the first time, stack the results vertically to create 2D matrices
                # Each row is attributed to the current iteration in time and each column is each cell in the comp domain
                # 2D matrices shape: (Nt,Nz)

                self.output_data['E'] = np.vstack((self.output_data['E'],self.sim_fields['E']))
                self.output_data['H'] = np.vstack((self.output_data['H'],self.sim_fields['H']))
                self.output_data['R'] = np.vstack((self.output_data['R'],self.sim_fields['R']))
                self.output_data['T'] = np.vstack((self.output_data['T'],self.sim_fields['T']))
                self.output_data['S'] = np.vstack((self.output_data['S'],self.sim_fields['S']))
                self.output_data['C'] = np.vstack((self.output_data['C'],self.sim_fields['C']))

            print(f"\rCurrent iteration: {curr_iteration}/{self.sim_param['Nt']}",end="")



        print("\n Computing FFT using numpy.rfft()...")
        # Compute the frequency response using numpy's fft methods
        self.output_data['S_numpy'] = np.fft.rfft(self.sim_source['Esrc'])
        self.output_data['R_numpy'] = np.fft.rfft(self.output_data['E'][:,1])/self.output_data['S_numpy']
        self.output_data['T_numpy'] = np.fft.rfft(self.output_data['E'][:,-2])/self.output_data['S_numpy']
        self.output_data['C_numpy'] = self.output_data['R_numpy'] + self.output_data['T_numpy']
        self.output_data['Freq_numpy'] = np.fft.rfftfreq(int(self.sim_param['Nt']),self.sim_param['df'])

        print("End of simulation.")
        fdtd_end_time = time.perf_counter()
        fdtd_duration = fdtd_end_time - fdtd_start_time
        self.output_data['algo_time'] = fdtd_duration

        return fdtd_duration

    def simulate_fdtd_schwarz(self):
        
        fdtd_schwarz_start_time = time.perf_counter()

        if self.sim_param['multithread'] == False:
            # Loop using for loop (for the FDTD time loop)
            print_breaker('minor')
            print("Starting FDTD-SCHWARZ algorithm...")
            numLoops = 0 #variable to track the number of loops
            isConverged = False

            # Convergence loop
            while(isConverged == False):
                
                if numLoops > 0:
                    self.update_sim_param(n_wavelength=numLoops,n_dim=numLoops) # Call update_sim_parameters() here 

                numLoops += 1
                
                # FDTD Algorithm
                # FDTD Time Loop
                for curr_iteration in range(int(self.sim_param['Nt'])):
                    print(f"\rCurrent Iteration: {curr_iteration}/{self.sim_param['Nt']}" ,end='')
                    for subdom_index in range(int(self.sim_param['n_subdom'])):
                        
                        subdom_start_time = time.perf_counter()
                        # Ghost cell transfer
                        # Left side = Magnetic Field Value
                        # Right side = Electric Field Value
                        #print(self.subdomains)
                        #print(f"Index: {subdom_index}  | Current subdom: {self.subdomains[subdom_index].subdom_param['id']}")
                        if subdom_index == 0:
                            
                            self.subdomains[subdom_index].r_ghost = self.subdomains[subdom_index + 1].subdom_fields['E'][ self.sim_param['overlap'] + 1]
                        
                        elif subdom_index == self.sim_param['n_subdom'] - 1:

                            self.subdomains[subdom_index].l_ghost = self.subdomains[subdom_index - 1].subdom_fields['H'][-1]
                        
                        else:
                            
                            self.subdomains[subdom_index].r_ghost = self.subdomains[subdom_index + 1].subdom_fields['E'][ self.sim_param['overlap'] + 1]
                            self.subdomains[subdom_index].l_ghost = self.subdomains[subdom_index - 1].subdom_fields['H'][-1]
                        

                       
                       
                       
                        #FDTD Algo main
                        self.subdomains[subdom_index].simulate_subdom(curr_iter = curr_iteration, boundary_condition = self.sim_param['boundary_condition'],excitation_method = self.sim_param['excitation_method'])

                        
                        # Measure the end time
                        subdom_end_time = time.perf_counter()

                        # Measure the execution time and append it to the output 
                        subdom_duration = subdom_end_time - subdom_start_time
                        self.subdomains[subdom_index].output['subdom_time'] += subdom_duration

                print(f"E: {self.subdomains[0].output['E']}" )
                self.plot_fields()
                print()
                print("Transferring boundary data from the subdomains")
                # Transfer boundary data
                for subdom_index in range(int(self.sim_param['n_subdom'])):

                    # Transferring the boundary data of each 2D matrices of the subdomains
                    # Left side index = overlap - 1 (since index starts at 0)
                    # Right side index = Subdom size - overlap
                    # Convention: We will always GET the data from the RIGHT ADJACENT SUBDOMAIN so
                    #  the last subdomain is not included

                    subdom_start_time = time.perf_counter()

                    if subdom_index < self.sim_param['n_subdom'] - 1:

                        # Transfer for the Electric field
                        buffer_E_cur = self.subdomains[subdom_index].output['E'][:,self.sim_param['subdomain_size'] - self.sim_param['overlap']]
                        buffer_E_right = self.subdomains[subdom_index + 1].output['E'][:,self.sim_param['overlap'] - 1]

                        # Transfer the field values
                        self.subdomains[subdom_index].output['E'][:,self.sim_param['subdomain_size'] - self.sim_param['overlap']] = self.subdomains[subdom_index + 1].output['E'][:,0]

                        self.subdomains[subdom_index + 1].output['E'][:,self.sim_param['overlap'] - 1] = self.subdomains[subdom_index].output['E'][:,-1]

                        self.subdomains[subdom_index].output['E'][:,-1] = buffer_E_right

                        self.subdomains[subdom_index + 1].output['E'][:,0] = buffer_E_cur

                        # Transfer for Magnetic field
                        buffer_H_cur = self.subdomains[subdom_index].output['H'][:,self.sim_param['subdomain_size'] - self.sim_param['overlap']]
                        buffer_H_right = self.subdomains[subdom_index + 1].output['H'][:,self.sim_param['overlap'] - 1]

                        # Transfer the field values
                        self.subdomains[subdom_index].output['H'][:,self.sim_param['subdomain_size'] - self.sim_param['overlap']] = self.subdomains[subdom_index + 1].output['H'][:,0]

                        self.subdomains[subdom_index + 1].output['H'][:,self.sim_param['overlap'] - 1] = self.subdomains[subdom_index].output['H'][:,-1]

                        self.subdomains[subdom_index].output['H'][:,-1] = buffer_H_right

                        self.subdomains[subdom_index + 1].output['H'][:,0] = buffer_H_cur

                    # Measure the additional execution time
                    subdom_end_time = time.perf_counter()

                    subdom_duration = subdom_end_time - subdom_start_time

                    self.subdomains[subdom_index].output['subdom_time'] += subdom_duration


                # Check for convergence
                print_breaker('minor')
                print(f"Checking convergence | Current number of loops (Schwarz loop): {numLoops}")
                isConverged = self.check_convergence(numLoops=numLoops)

             # Construct main 2D array
            self.reconstruct_output_matrix()




        elif self.sim_param['multithread'] == True: # if we will be using OpenMP or multiprocessing module
            pass


        fdtd_schwarz_end_time = time.perf_counter()
        fdtd_schwarz_duration = fdtd_schwarz_end_time - fdtd_schwarz_start_time
        self.output_data['algo_time'] = fdtd_schwarz_duration

        
        return fdtd_schwarz_duration


    def plot_fields(self):
        rows = int(self.sim_param['n_subdom']/2)
        plt.ion()
        fig = plt.figure(1,[10,6])
        axs = []
        lineE = []
        lineH = []

        for subdom_index in range(int(self.sim_param['n_subdom'])):

            axs.append(fig.add_subplot(rows,rows,subdom_index+1))

            axs[subdom_index].set(xlabel="Cells",ylabel="Levels", title=f"Subdomain {subdom_index}",ylim=[-2,2])

            lineE.append(axs[subdom_index].plot(self.subdomains[subdom_index].comp_dom['z'],self.subdomains[subdom_index].output['E'][0,:]))
            lineH.append(axs[subdom_index].plot(self.subdomains[subdom_index].comp_dom['z'],self.subdomains[subdom_index].output['H'][0,:]))
            fig.suptitle(f"FDTD Simulation Iteration: {0}/{int(self.sim_param['Nt'])} [ boundary_cond: {self.sim_param['boundary_condition']} | excitation_method: {self.sim_param['excitation_method']} ]")

        for i in range(int(self.sim_param['Nt'])):

            for subdom_index in range(int(self.sim_param['n_subdom'])):

                lineE[subdom_index][0].set_ydata(self.subdomains[subdom_index].output['E'][i,:])
                lineH[subdom_index][0].set_ydata(self.subdomains[subdom_index].output['H'][i,:])

            fig.suptitle(f"FDTD Simulation Iteration: {i}/{int(self.sim_param['Nt'])} [ boundary_cond: {self.sim_param['boundary_condition']} | excitation_method: {self.sim_param['excitation_method']} ]")

            fig.canvas.draw()
            fig.canvas.flush_events()


    def reconstruct_output_matrix(self):
        
        # This function compiles all of the computed field values into a single 1D vector.
        # The reconstructed 1D vector will be used in plotting and FFT computations.
        # At this point, the vecotrs sim_fields.E and sim_fields.H has no size so we will use hstack()...
        
        start = int(self.sim_param['overlap'])
        stop = int(self.sim_param['non_overlap'] + self.sim_param['overlap'])

        for subdom_index in range(int(self.sim_param['n_subdom'])):

            if subdom_index == 0:

                self.output_data['E'] = self.subdomains[subdom_index].output['E'][:,start:]
                self.output_data['H'] = self.subdomains[subdom_index].output['H'][:,start:]
                

            elif subdom_index == self.sim_param['n_subdom'] - 1:
                
                self.output_data['E'] = np.hstack((self.output_data['E'],self.subdomains[subdom_index].output['E'][:,start:stop]))
                self.output_data['H'] = np.hstack((self.output_data['H'],self.subdomains[subdom_index].output['H'][:,start:stop]))

            else:

                self.output_data['E'] = np.hstack((self.output_data['E'],self.subdomains[subdom_index].output['E'][:,start:]))
                self.output_data['H'] = np.hstack((self.output_data['H'],self.subdomains[subdom_index].output['H'][:,start:]))


        return True


    def check_convergence(self,numLoops= 0):
        
        # In checking for converegence, we can also use isclose() or allclose() 
        # method where instead of L2 norm, we are just using element-wise 
        # comparison of values based on a threshold (that can be adjusted)
        # (To be implemented)

        # Layout of the error vector (Row layout):
        # | Error of sub 1 and sub 2 | Error of sub 2 and sub 3 | Error of sub 3 and 4 | ..... | Error of sub n-1 and n |
        # [              error1,                 error2,                    error3,      .....             error n-1    ]

        # Error vector length is n-1 because each error takes into account two subdomains
        # Each row corresponds to 1 iteration in the Schwarz method so 

        # # of cols = n-1 error values
        # # of rows = numLoops (number of iterations in the Schwarz method)

        isConverged = False

        # Initialize the error vectors
        self.output_data['E_error'] = np.zeros([self.sim_param['n_subdom'] - 1])
        self.output_data['H_error'] = np.zeros([self.sim_param['n_subdom'] - 1])
        

        for subdom_index in range(int(self.sim_param['n_subdom'] -1)):

            # Get the difference of the two overlapping region (in the same positions)
            E_sub = self.subdomains[subdom_index].output['E'][:,int(self.sim_param['subdomain_size'] - self.sim_param['overlap']):] - self.subdomains[subdom_index + 1].output['E'][:,:int(self.sim_param['overlap'])]
            H_sub = self.subdomains[subdom_index].output['H'][:,int(self.sim_param['subdomain_size'] - self.sim_param['overlap']):] - self.subdomains[subdom_index + 1].output['H'][:,:int(self.sim_param['overlap'])]

            # Get the errors by using the norm()
            self.output_data['E_error'][subdom_index] = np.linalg.norm(E_sub,ord=2)
            self.output_data['H_error'][subdom_index] = np.linalg.norm(H_sub,ord=2)

        # Printing error values
        print_breaker('major')
        print("Error Values:")
        print(f"E: {self.output_data['E_error']} | H: {self.output_data['H_error']}")
        # Save each error vector in the 2D matrix
        if numLoops == 1: # At the 1st call of this function...

            self.output_data['E_error_list'] = self.output_data['E_error']
            self.output_data['H_error_list'] = self.output_data['H_error']
            

        else:

            self.output_data['E_error_list'] = np.vstack((self.output_data['E_error_list'] ,self.output_data['E_error']))
            self.output_data['H_error_list'] = np.vstack((self.output_data['H_error_list'] ,self.output_data['H_error']))


        return isConverged

    
 

        

    def save_sim(self,name="",type="hdf5",output_dir="",username="",description=""):
        
        # Check if the simulation is done before calling this function
        try:
            assert self.output_data != None
      

        except AssertionError:
            print_breaker('minor')
            print("ERROR: Output data not detected.")
            return None

        # Getting the current datetime
        now = datetime.datetime.now()
        self.date_str = now.strftime("%Y-%m-%d")

        # Verify if the directory exists
        try:
            assert os.path.isdir(output_dir)
        except AssertionError:
            print_breaker('minor')
            print("ERROR: Directory is not valid.")
            print("The current working directory will be used instead.")
            print(f"Current working directory: {output_dir}")
            output_dir = os.getcwd()

        
        # Saving data in HDF5 format
        if type == "hdf5":
            print_breaker("major")
            print(f"Saving the simulation data. | File type: {type}")
            filename = output_dir + self.date_str + "_" + name + ".h5" 
            with h5py.File(filename,mode='w') as file:

                # Storing simulation metadata in the input group
                input_grp = file.create_group("input")
                print("Saving metadata -----", end="")
                file['/'].attrs['date'] = self.date_str
                file['/'].attrs['user'] = username
                file['/'].attrs['description'] = description
                print("---> Saved!")

                # Storing the input data
                print("Saving input data -----",end="")
                for key,val in self.input_data.items():
                    input_grp.create_dataset(key,data=val,compression="gzip")
                print("---> Saved!")

                print("Saving comp domain vectors -----",end="")
                comp_domain_grp = file.create_group("comp_domain")
                for key,val in self.comp_domain.items():
                    comp_domain_grp.create_dataset(key,data=val,compression="gzip")
                print("---> Saved!")

                # Storing simulation parameters as attributes of the comp_domain group
                print("Saving simulation parameters -----",end="")
                for key,val in self.sim_param.items():
                    comp_domain_grp.attrs[key] = val
                print("---> Saved!")

               
                # Storing the output data in output group
                output_grp = file.create_group("output")
                subdom_grp = output_grp.create_group("subdomain")
                print("Saving output data -----",end="")
                for key,val in self.output_data.items():
                    output_grp.create_dataset(key,data=val,compression="gzip")
                print("---> Saved!")

        # Saving data in NPY format
        elif type == "npy":
            # 3D matrix total size: (face,row,col) = (5,Nt,Nz)
            # npy file will save a 3D matrix containing everything from the simulation parameter to the simulated data.
            # 1st 2D matrix (output[0,:,:]) will contain the simulation parameters (Nz, dz, etc.) and the time and source vectors.
            # |- In the 1st matrix, the 1st row corresponds to the sim params, specifically in this order from left to right: 
            # |   |- fmax, Nt, Nz, dz, dt, n_model, l_spacer, r_spacer, inj_point, sim_time
            # |- the 2nd row corrseponds to the source data (Esrc) and 3rd row to the (Hsrc)
            # |- 4th row: t_E, 5th row: t_H, 6th row: t
            # 2nd 2D matrix will contain E field simulation data
            # 3rd 2D matrix will contain H field simulation data
            
            row = int(self.sim_param['Nt'])
            col = int(self.sim_param['Nz'])
            output_data = np.zeros((5,row,col))
            for i in range(5):

                if i == 0:
                    output_data[i,0,0] = self.sim_param['fmax']
                    output_data[i,1,0] = self.sim_param['Nt']
                    output_data[i,2,0] = self.sim_param['Nz']
                    output_data[i,3,0] = self.sim_param['dz']
                    output_data[i,4,0] = self.sim_param['dt']
                    output_data[i,5,0] = self.sim_param['n_model']
                    output_data[i,6,0] = self.sim_param['left_spacer']
                    output_data[i,7,0] = self.sim_param['right_spacer']
                    output_data[i,8,0] = self.sim_param['inj_point']
                    output_data[i,9,0] = self.sim_param['sim_time']


                    output_data[i,:,1] = self.output_data['Esrc']
                    output_data[i,:,2] = self.output_data['Hsrc']
                    output_data[i,:,3] = self.output_data['t_E']
                    output_data[i,:,4] = self.output_data['t_H']
                    output_data[i,:,5] = self.output_data['t']
                    

                elif i == 1:
                    output_data[i,:,:] = self.output_data['E']
                elif i == 2:
                    output_data[i,:,:] = self.output_data['H']

            print_breaker("major")
            print(f"Saving the simulation data. | File type: {type}",end="")
            filename = output_dir + self.date_str + "_" + name + ".npy"
            np.save(filename,output_data) 
            print("---> Saved!")

        # Saving data in CSV format
        elif type == "csv":

            # Creates several files: electric field, magnetic field, source (E and H), reflectance, and transmittance
            csv_files = ['electric', 'magnetic', 'source', 'reflectance', 'transmittance']
            print_breaker("major")
            print(f"Saving the simulation data. | File type: {type}",end="")
            for i in range(len(csv_files)):
                filename = output_dir + self.date_str + "_" + name + "_" + csv_files[i]
                if i == 0:
                    np.savetxt(filename,self.output_data['E'],delimiter=',')
                elif i == 1:
                    np.savetxt(filename,self.output_data['H'],delimiter=',')
                elif i == 2:
                    source_output = np.hstack((self.output_data['Esrc'],self.output_data['Hsrc']))
                    np.savetxt(filename,source_output,delimiter=',')

                elif i == 3:
                    np.savetxt(filename,self.output_data['R'],delimiter=',')

                elif i == 4:
                    np.savetxt(filename,self.output_data['T'],delimiter=',')
            print("---> Saved!")
