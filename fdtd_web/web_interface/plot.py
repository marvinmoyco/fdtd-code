from asyncio import new_event_loop
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
import h5py
import sys
import os

from FDTD.simulation import *
#Converter function from numpy ndarray to str object
def np_to_str(input):
    return np.array_str(np.array(input).astype(str))

class Plot:
    # This is a Plot object that generates an animated plot (either in video or HTML format)
    # for the FDTD simulations. The different input parameters are defined as follows:
    # simulation - simulation object created (simulation must be done before creating a plot object)
    # type - tells whether the output is a video ("mp4") or browser-based ("html")
    # n_frame - integer value to get samples from the data (e.g. n_frame=50, getting every 50th sample)
    # save - boolean flag to tell the object whether to save the generated plot or not
    # output_path - file path where the generated plot will be saved
    # read - boolean flag to tell whether the object uses the simulation object (False) or hdf5 file (True)
    # input - input name of the hdf5 file (used only when read argument is True)
    # output - output name of the generated plot (used only when the save flag is True)

    def __init__(self,
                simulation,
                type="html",
                n_frame=5,
                save=True,
                output_path="./",
                read=False,
                input="",
                output="output"):

        #Initialize instance attributes
        self.sim_fields = {}
        self.FFT = {}
        self.source = {}
        self.sim_param = {}
        self.comp_domain = {}
        self.read = read
        self.save = save
        
        self.new_Nt = 0

        if self.save == True:
            self.type = type
            self.output = output
            self.output_path = output_path
        
        if self.read == True:
            self.input = input



        if read == True:
            #Read Electric field
            data = h5py.File(input,mode='r')
            print(f" HDF5 file keys: {data.keys()}")
            print(f"Nz: {np.array(data['sim param']['Nz'])}")
            #z = np.arange(start=1,stop=34)
            z = data['comp domain']['z']
            algo = np_to_str(data['sim param']['algorithm'])

            n = data['comp domain']['n'][:]
            Nt = int(np.ceil(np.array(data['sim param']['Nt'])))
            freq_axis = data['/output/freq_range'][:].T
            dz = float(np.array(data['sim param']['dz']))
            spacer = int(np.array(data['sim param']['left spacer']))*dz
            dt = float(data['sim param']['dt'][()])
            Nz = data['sim param']['Nz'][()]
            t = np.linspace(0,Nt,Nt)
            #Read the main output data
            E = np.nan_to_num(data['output']['E'][:,:])
            H = np.nan_to_num(data['output']['H'][:,:])

            #t = data['source']['t'][:]
            t_E = data['source']['t_E'][:]
            Esrc = data['source']['Esrc'][:]
            t_H = data['source']['t_H'][:]
            Hsrc = data['source']['Hsrc'][:]
            fmax = data['sim param']['fmax'][()]
            #Compute for the FFT here...
            #Toggle this variable if you want to recreate the FFT of the simulation...
            recreate_fft = True

            if recreate_fft == True:
                freq_axis = np.linspace(0,fmax,Nt)
                Kernel_freq = np.exp(-1j*2*np.pi*dt*freq_axis)
                rowR = np.zeros(Kernel_freq.shape)
                rowT = np.zeros(Kernel_freq.shape)
                rowC = np.zeros(Kernel_freq.shape)
                rowS = np.zeros(Kernel_freq.shape)
                for i in range(Nt):
                    print(f"Current iteration: {i}/{Nt}")
                    rowR = rowR + np.power(Kernel_freq,i)*E[i,1]
                    rowT = rowT + np.power(Kernel_freq,i)*E[i,-2]
                    rowS = rowS +  np.power(Kernel_freq,i)*Esrc[i]

                    rR = np.power(np.abs(rowR/rowS),2)
                    rT = np.power(np.abs(rowT/rowS),2)
                    rC = rR+rT

                    if i == 0:
                        R = rR
                        T = rT
                        C = rC
                    #elif i % 10 == 0:
                    #   new_R = np.vstack((new_R,rR))
                    #  new_T = np.vstack((new_T,rT))
                    # new_C = np.vstack((new_C,rC))
                    else:
                        R = np.vstack((R,rR))
                        T = np.vstack((T,rT))
                        C = np.vstack((C,rC))
                else:
                    R = np.nan_to_num(data['output']['reflectance'][:,:])
                    T = np.nan_to_num(data['output']['transmittance'][:,:])
                    C = np.nan_to_num(data['output']['conservation_of_energy'][:,:])
        else:
            # Get the simulation parameters that will be used in the plot
            self.sim_param['Nz'] = simulation.sim_param['Nz']
            self.sim_param['Nt'] = int(simulation.sim_param['Nt'])
            self.sim_param['dz'] = simulation.sim_param['dz']
            self.sim_param['dt'] = simulation.sim_param['dt']
            
            # Get the simulation EM fields
            # self.sim_fields['E'] = simulation.output_data['E']
            # self.sim_fields['H'] = simulation.output_data['H']
            self.sim_fields['z'] = simulation.comp_domain['z']

            # Get FFT Response
            # self.FFT['R'] = simulation.output_data['reflectance']
            # self.FFT['T'] = simulation.output_data['transmittance']
            # self.FFT['C'] = simulation.output_data['conservation_of_energy']

            for i in range(self.sim_param['Nt']):
                if i == 0:

                    self.sim_fields['E'] = simulation.output_data['E'][i,:]
                    self.sim_fields['H'] = simulation.output_data['H'][i,:]
                    self.sim_fields['R'] = simulation.output_data['R'][i,:]
                    self.sim_fields['T'] = simulation.output_data['T'][i,:]
                    self.sim_fields['C'] = simulation.output_data['C'][i,:]
                    self.new_Nt += 1

                elif i % n_frame == 0:
                    
                    self.sim_fields['E'] = np.vstack((self.sim_fields['E'],simulation.output_data['E'][i,:]))
                    self.sim_fields['H'] = np.vstack((self.sim_fields['H'],simulation.output_data['H'][i,:]))
                    self.sim_fields['R'] = np.vstack((self.sim_fields['R'],simulation.output_data['R'][i,:]))
                    self.sim_fields['T'] = np.vstack((self.sim_fields['T'],simulation.output_data['T'][i,:]))
                    self.sim_fields['C'] = np.vstack((self.sim_fields['C'],simulation.output_data['C'][i,:]))

                    self.new_Nt += 1
                elif i == self.sim_param['Nt'] -1:
                    
                    self.sim_fields['E'] = np.vstack((self.sim_fields['E'],simulation.output_data['E'][i,:]))
                    self.sim_fields['H'] = np.vstack((self.sim_fields['H'],simulation.output_data['H'][i,:]))
                    self.sim_fields['R'] = np.vstack((self.sim_fields['R'],simulation.output_data['R'][i,:]))
                    self.sim_fields['T'] = np.vstack((self.sim_fields['T'],simulation.output_data['T'][i,:]))
                    self.sim_fields['C'] = np.vstack((self.sim_fields['C'],simulation.output_data['C'][i,:]))

                    self.new_Nt += 1
            total_frames,_ = self.sim_fields['E'].shape
            self.total_frames = total_frames
            self.l_spacer = simulation.sim_param['left_spacer']*self.sim_param['dz']
            self.r_spacer = simulation.sim_param['right_spacer']
            self.input_layers = simulation.input_data['layer_size']
            self.comp_domain['freq_range'] = simulation.output_data['Freq_range']


            
            self.sim_param['sim_date'] = simulation.date_str
            self.sim_param['boundary_cond'] = simulation.sim_param['boundary_condition']
            self.sim_param['excitation_method'] = simulation.sim_param['excitation_method']


            # Get EM Source
            self.source['Esrc'] = simulation.sim_source['Esrc']
            self.source['Hsrc'] = simulation.sim_source['Hsrc']
            self.source['t'] = simulation.sim_source['t']
            self.source['t_E'] = simulation.sim_source['t_E']
            self.source['t_H'] = simulation.sim_source['t_H']
            self.source['type'] = simulation.sim_param['source_type']
            
            # Getting comp domain 
            #self.comp_domain['mu'] = simulation.comp_domain['mu']
            #self.comp_domain['epsilon'] = simulation.comp_domain['epsilon']
            self.comp_domain['n'] = simulation.comp_domain['n']
            self.comp_domain['z'] = simulation.comp_domain['z']



    def plot_html(self):
        # Plot the data using Plotly.

        # Creating a figure object
        plot_fig = make_subplots(rows=3, cols=1, subplot_titles = ('FDTD Simulation','Frequency Response', f'Source: { self.source["type"] }'))

        range_vals = [np.max([np.amax(self.sim_fields['E']),np.amax(self.sim_fields['H'])]),np.min([np.amin(self.sim_fields['E']),np.amin(self.sim_fields['H'])])]
    
        # Create the layers in the simulation
        start = self.r_spacer*self.sim_param['dz']
        end = self.r_spacer*self.sim_param['dz']
        n_model =len(self.input_layers)
        for i in range(n_model):

            end += self.input_layers[i]

            plot_fig.add_trace(go.Scatter(
                            x = [start, start,end,end,start],
                            y = [range_vals[1],range_vals[0],range_vals[0],range_vals[1],range_vals[1]],
                            hoverinfo='text',
                            fill='toself',
                            opacity=0.3,
                            text=f'Layer {i+ 1} \n Refractive index: {self.comp_domain["n"][i]:.03f}',
                            showlegend=False
            ), row= 1, col = 1)

            start = end

        # Add traces in the 1st plot (FDTD Simulation)
        # For the Electric Field plot
        plot_fig.add_trace(go.Scatter(
                            x = self.comp_domain['z'],
                            y = self.sim_fields['E'][0,:],
                            mode = 'lines',
                            legendgroup = 'Electric Field',
                            hovertemplate = 'x: %{x} <br> y:%{y}',
                            name= 'Electric Field',
                            showlegend=True
            ), row= 1, col = 1)

        # For the Magnetic Field plot
        plot_fig.add_trace(go.Scatter(
                            x = self.comp_domain['z'],
                            y = self.sim_fields['H'][0,:],
                            mode = 'lines',
                            legendgroup = 'Magnetic Field',
                            hovertemplate = 'x: %{x} <br> y:%{y}',
                            name= 'Magnetic Field',
                            showlegend=True
            ), row= 1, col = 1)

        # Adding traces for the 2nd plot (Frequency response)
        
        #For the Reflectance plot
        plot_fig.add_trace(go.Scatter(
                            x = self.comp_domain['freq_range'],
                            y = self.sim_fields['R'][0,:],
                            mode = 'lines',
                            legendgroup = 'Reflectance',
                            hovertemplate = 'x: %{x} <br> y:%{y}',
                            name= 'Reflectance',
                            showlegend=True
            ), row= 2, col = 1)

        #For the Transmittance plot
        plot_fig.add_trace(go.Scatter(
                            x = self.comp_domain['freq_range'],
                            y = self.sim_fields['T'][0,:],
                            mode = 'lines',
                            legendgroup = 'Transmittance',
                            hovertemplate = 'x: %{x} <br> y:%{y}',
                            name= 'Transmittance',
                            showlegend=True
            ), row= 2, col = 1)

        #For the Conservation of energy
        plot_fig.add_trace(go.Scatter(
                            x = self.comp_domain['freq_range'],
                            y = self.sim_fields['C'][0,:],
                            mode = 'lines',
                            legendgroup = 'Conservation of Energy',
                            hovertemplate = 'x: %{x} <br> y:%{y}',
                            name= 'Conservation of Energy',
                            showlegend=True
            ), row= 2, col = 1)

        # For the 3rd subplot (Input Excitation)

        #For the Electric Field plot of the input excitation
        plot_fig.add_trace(go.Scatter(
                            x = self.source['t_E'],
                            y = self.source['Esrc'],
                            mode = 'lines',
                            legendgroup = 'Source (Electric Field)',
                            hovertemplate = 'x: %{x} <br> y:%{y}',
                            name= 'Source (Electric Field)',
                            showlegend=True
            ), row= 3, col = 1)
    
        #For the Magnetic Field plot of the input excitation
        plot_fig.add_trace(go.Scatter(
                            x = self.source['t_H'],
                            y = self.source['Hsrc'],
                            mode = 'lines',
                            legendgroup = 'Source (Magnetic Field)',
                            hovertemplate = 'x: %{x} <br> y:%{y}',
                            name= 'Source (Magnetic Field)',
                            showlegend=True
            ), row= 3, col = 1)

        # Create the necessary layout for the plot
        title = f"FDTD Simulation [Date: {self.sim_param['sim_date']} | Source Type: {self.source['type']} | Boundary Condition: {self.sim_param['boundary_cond']} | Excitation Method: {self.sim_param['excitation_method']}]"
        plot_fig.update_layout(title_text = title )
        
        # Changing the axes labels

        #For x-axis
        plot_fig.update_xaxes(title_text="Computationaldomain [m]",row=1,col=1)
        plot_fig.update_xaxes(title_text="Frequency [Hz]",row=2,col=1)
        plot_fig.update_xaxes(title_text= "Time [s]",row=3,col=1)

        # For y-axis
        plot_fig.update_yaxes(title_text="Level",row=1,col=1)
        plot_fig.update_yaxes(title_text="Magnitude",row=2,col=1)
        plot_fig.update_yaxes(title_text= "Level",row=3,col=1)


        # Adding 'frames' to each traces
        end_indices = len(self.input_layers)
        frames=[dict(name=i,
                    data=[
                    go.Scatter(y=self.sim_fields['E'][i,:]),
                    go.Scatter(y=self.sim_fields['H'][i,:]),
                    go.Scatter(y=self.sim_fields['R'][i,:]),
                    go.Scatter(y=self.sim_fields['T'][i,:]),
                    go.Scatter(y=self.sim_fields['C'][i,:]),
                    go.Scatter(y=self.source['Esrc']),
                    go.Scatter(y=self.source['Hsrc'])
                    ],
                    #Added 4 to the end indices since there 4 traces after (E,H,R,T)
                    traces = [x for x in range(end_indices,end_indices+7)]#This is done because the rectangles are traces (to have hover information)
        ) for i in range(self.new_Nt)]

        # Adding button and slider
        updatemenus = [dict(type='buttons',
                        buttons=[dict(label='Play',
                                method='animate',
                                args=[[f'{i}' for i in range(self.new_Nt)],
                                        dict(frame=dict(duration=0,redraw=False),
                                        transition=dict(duration=0),
                                        easing='linear',
                                        fromcurrent=True,
                                        mode='immediate')])],
                        direction='left',
                        pad=dict(r=10,t=85),
                        showactive= True,x=0.1,y=0,xanchor='right',yanchor='top'

        )]

        sliders = [{'yanchor': 'top',
            'xanchor': 'left', 
            'currentvalue': {'font': {'size': 16}, 'prefix': 'Iteration: ', 'visible': True, 'xanchor': 'right'},
            'transition': {'duration': 0.0, 'easing': 'linear-in-out'},
            'pad': {'b': 10, 't': 50}, 
            'len': 0.9, 'x': 0.1, 'y': 0, 
            'steps': [{'args': [[i], {'frame': {'duration': 0.0, 'easing': 'linear-in-out', 'redraw': False},
                                    'transition': {'duration': 0, 'easing': 'linear-in-out'}}], 
                    'label': f"{i}/{self.new_Nt}", 'method': 'animate'} for i in range(self.new_Nt)       
                    ]}]

        # Update the figure to include the menu and slider
        plot_fig.update(frames=frames)
        plot_fig.update_layout(updatemenus=updatemenus,sliders=sliders)

        # Saving the plot
        if self.save == True:
            numFile = 2
            while True:
                isFileExists = os.path.exists(self.output_path + self.output + '.' + self.type)
                if isFileExists ==False:
                    break
                numFile += 1
                num_str = str(numFile).zfill(3)
                self.output = self.output + '_' + num_str + '.' + self.type
            
            plot_fig.write_html(self.output_path + self.output,auto_play=False)

        else:

            plot_fig.show('browser')

        

    #Converter function from numpy ndarray to str object
    def np_to_str(input):
        return np.array_str(np.array(input).astype(str))