# Main Reference: https://chart-studio.plotly.com/~empet/15243/animating-traces-in-subplotsbr/#/


from asyncio import new_event_loop
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
import h5py
import sys

#Converter function from numpy ndarray to str object
def np_to_str(input):
    return np.array_str(np.array(input).astype(str))

filename = ''
#Check the length of input arguments
if len(sys.argv) == 3:
    filename = sys.argv[1] + sys.argv[2]

#Read Electric field
data = h5py.File(filename,mode='r')
print(f" HDF5 file keys: {data.keys()}")
print(f"Nz: {np.array(data['sim param']['Nz'])}")
#z = np.arange(start=1,stop=34)
z = data['comp domain']['z']
algo = np_to_str(data['sim param']['algorithm'])

#Read the main output data
E = np.nan_to_num(data['output']['E'][:,:])
H = np.nan_to_num(data['output']['H'][:,:])
R = np.nan_to_num(data['output']['reflectance'][:,:])
T = np.nan_to_num(data['output']['transmittance'][:,:])
C = np.nan_to_num(data['output']['conservation_of_energy'][:,:])

#Initialize number of subdomain to 1 (default since fdtd has 1 thread only)
num_subdom = 1
# Creating an empty list that will contain dictionaries (list of dictionaries where each dictionary is about the data of the subdomain)
subdomain_data = []
#Read the 
if algo == "fdtd-schwarz":
    num_subdom = data['sim param']['num subdomains'][()]
    print(f"Number of subdomains in {algo}: {num_subdom} threads/subdomains")

    #Read all of the simulation data in each subdomain
    for subdom_index in range(num_subdom):
        subdom  = {}
        subdom['Subdomain id'] = subdom_index
        subdom['E'] = np.nan_to_num(data['output']['subdomain'][str(subdom_index)]['E'][:,:])
        subdom['H'] = np.nan_to_num(data['output']['subdomain'][str(subdom_index)]['H'][:,:])
        subdomain_data.append(subdom) 

print(f"Field shapes (before adjusting): E: {E.shape} | H: {H.shape} | R: {R.shape} | T: {T.shape} | C: {C.shape}")
#Read source data
t = data['source']['t'][:]
t_E = data['source']['t_E'][:]
Esrc = data['source']['Esrc'][:]
t_H = data['source']['t_H'][:]
Hsrc = data['source']['Hsrc'][:]
#print(Esrc,Hsrc)
##print(R,T)
#print(C)
#print(E,H)

#Read the config parameters
sim_date = np_to_str((data['metadata']['date']))
source_type = np_to_str(data['source']['type'])
boundary_cond = np_to_str(data['/sim param/boundary cond'])
excitation_method = np_to_str(data['/sim param/excitation method'])

n = data['comp domain']['n'][:]
Nt = int(np.ceil(np.array(data['sim param']['Nt'])))
freq_axis = data['/output/freq_range'][:]
dz = float(np.array(data['sim param']['dz']))
spacer = int(np.array(data['sim param']['left spacer']))*dz


#print(z.shape,E.shape,Nt)
input_layer = data['input']['layer size'][:]
mu = data['input']['magnetic permeability'][:]
epsilon = data['input']['electric permittivity'][:]
#print(input_layer,spacer)

multithreading = ""
if data['sim param']['multithreading'][()] == 0:
    multithreading = "False"
elif data['sim param']['multithreading'][()] == 1:
    multithreading = "True"

new_Nt = 0
new_subdomain_data = []
#Adjust the frames here
for i in range(Nt):
    if i == 0:
        new_E = E[i,:]
        new_H = H[i,:]
        new_R = R[i,:]
        new_T = T[i,:]
        new_C = C[i,:]
        #For the subdomains
        if algo == "fdtd-schwarz":
            for subdom_index in range(num_subdom):
                new_subdom = {}
                new_subdom['Subdomain id'] = subdomain_data[subdom_index]['Subdomain id']
                new_subdom['E'] = subdomain_data[subdom_index]['E'][i,:]
                new_subdom['H'] = subdomain_data[subdom_index]['H'][i,:]
                new_subdomain_data.append(new_subdom)



        new_Nt += 1
    elif i % 10 == 0:
        new_E = np.vstack((new_E,E[i,:]))
        new_H = np.vstack((new_H,H[i,:]))
        new_R = np.vstack((new_R,R[i,:]))
        new_T = np.vstack((new_T,T[i,:]))
        new_C = np.vstack((new_C,C[i,:]))
        #For the subdomains
        if algo == "fdtd-schwarz":
            for subdom_index in range(num_subdom):
                
                new_subdomain_data[subdom_index]['E'] = np.vstack((new_subdomain_data[subdom_index]['E'],
                                                        subdomain_data[subdom_index]['E'][i,:]))
                new_subdomain_data[subdom_index]['H'] = np.vstack((new_subdomain_data[subdom_index]['H'],
                                                        subdomain_data[subdom_index]['H'][i,:]))
                
        new_Nt += 1
    elif i == Nt - 1:
        new_E = np.vstack((new_E,E[i,:]))
        new_H = np.vstack((new_H,H[i,:]))
        new_R = np.vstack((new_R,R[i,:]))
        new_T = np.vstack((new_T,T[i,:]))
        new_C = np.vstack((new_C,C[i,:]))
        #For the subdomains
        if algo == "fdtd-schwarz":
            for subdom_index in range(num_subdom):
                
                new_subdomain_data[subdom_index]['E'] = np.vstack((new_subdomain_data[subdom_index]['E'],
                                                        subdomain_data[subdom_index]['E'][i,:]))
                new_subdomain_data[subdom_index]['H'] = np.vstack((new_subdomain_data[subdom_index]['H'],
                                                        subdomain_data[subdom_index]['H'][i,:]))
                
        new_Nt += 1


print(f"Field shapes (after adjusting): E: {new_E.shape} | H: {new_H.shape} | R: {new_R.shape} | T: {new_T.shape} | C: {new_C.shape}")
print(f"New Nt (total number of frames to be rendered): {new_Nt}")
# Only print the subdomain plots (another HTML file) if the algo is fdtd-schwarz

    #Create subplots of the main simulation data (E and H), 
fig = make_subplots(rows=3, cols=1, subplot_titles = ('FDTD Simulation', 'Frequency Response (Reflectance and Transmittance)',f'Source Plot (type: {source_type})'))

#Get the total maximum for E and H fields
max_point = np.max([np.amax(new_E),np.amax(new_H)])
min_point = np.min([np.amin(new_E),np.amin(new_H)])
#Get the starting indices for the shapes
x_start = spacer
x_end = spacer
for index in range(len(input_layer)):
    x_end += input_layer[index]

    # Add shapes
    #fig.add_shape(type="rect",
    #    xref="x", yref="y",
    #    x0=x_start, y0=max_point,
    #    x1=x_end, y1=min_point,
    #    line=dict(
    #        #color="RoyalBlue",
    #        width=1,
    #        dash='dashdot'
    #    ),
    #    fillcolor="forestgreen",
    #    layer='below',
    #    opacity=0.5
    #)
    fig.add_trace(go.Scatter(
                x= [x_start, x_start, x_end, x_end, x_start],
                y= [min_point, max_point, max_point, min_point, min_point],
                hoverinfo='text',
                fill = 'toself',
                opacity=0.3,
                #line_color= 'rgb(255, 79, 38)',
                text= f'Layer {index + 1} \n Refractive index: {np.sqrt(mu[index]*epsilon[index]):.03f}',
                showlegend= False), row=1, col=1)

    x_start = x_end

fig.add_trace(go.Scatter(
                x= t_E,
                y= Esrc,
                mode = 'lines',
                hovertemplate="x: %{x} <br> y: %{y}",
                legendgroup= 'Input Source (Electric Field)',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Input Source (Electric Field)',
                showlegend= True), row=3, col=1)

fig.add_trace(go.Scatter(
                x= t_H,
                y= Hsrc,
                mode = 'lines',
                hovertemplate="x: %{x} <br> y: %{y}",
                legendgroup= 'Input Source (Magnetic Field)',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Input Source (Magnetic Field)',
                showlegend= True), row=3, col=1)


#Add traces in subplot 1
fig.add_trace(go.Scatter(
                x= z,
                y= new_E[0,:],
                mode = 'lines',
                legendgroup= 'Electric Field',
                hovertemplate="x: %{x} <br> y: %{y}",
                #line_color= 'rgb(255, 79, 38)',
                name= 'Electric Field',
                showlegend= True), row=1, col=1)
                
fig.add_trace(go.Scatter(
                x= z,
                y= new_H[0,:],
                mode = 'lines',
                hovertemplate="x: %{x} <br> y: %{y}",
                legendgroup= 'Magnetic Field',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Magnetic Field',
                showlegend= True), row=1, col=1)



#Add traces in subplot 2
fig.add_trace(go.Scatter(
                x= freq_axis,
                y= new_R[0,:],
                mode = 'lines',
                hovertemplate="x: %{x} <br> y: %{y}",
                legendgroup= 'Reflectance',
                #line_color= 'rgb(255, 79, 38)',
                name = 'Reflectance',
                showlegend= True), row=2, col=1)

fig.add_trace(go.Scatter(
                x= freq_axis,
                y= new_T[0,:],
                mode = 'lines',
                hovertemplate="x: %{x} <br> y: %{y}",
                legendgroup= 'Transmittance',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Transmittance',
                showlegend= True), row=2, col=1)

fig.add_trace(go.Scatter(
                x= freq_axis,
                y= new_C[0,:],
                mode = 'lines',
                line=dict(dash='dash'),
                hovertemplate="x: %{x} <br> y: %{y}",
                legendgroup= 'Conservation of Energy',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Conservation of Energy',
                showlegend= True), row=2, col=1)


print(len(fig.data))



#width=1024, height=720,
print("Setting the titles for the plots....")

fig.update_layout(title_text=f"FDTD Simulation [Date: {sim_date} | Algorithm: {algo} | Multithreading: {multithreading} | Source Type: {source_type} | Boundary Condition: {boundary_cond} | Excitation Method: {excitation_method}]")
#Adjust the axes of the three subplots
fig.update_xaxes(title_text="Computational Domain (m)", row=1,col=1)
fig.update_xaxes(title_text="Frequency (Hz)",row=2,col=1)
fig.update_xaxes(title_text="Time (s)",row=3,col=1)


fig.update_yaxes(title_text="Level",range=[-1,1],row=1,col=1)
fig.update_yaxes(title_text="Magnitude",range=[0,1],row=2,col=1)
fig.update_yaxes(title_text="Level",range=[-1,1],row=3,col=1)

end_indices = len(input_layer) 
#Add the frames of each traces
print("Adding the frames of the plots...")
frames=[dict(name=i,
            data=[go.Scatter(y=Esrc),
            go.Scatter(y=Hsrc),
            go.Scatter(y=new_E[i,:]),#update the E-Field
            go.Scatter(y=new_H[i,:]),#update the H-Field
            go.Scatter(y=new_R[i,:]),
            go.Scatter(y=new_T[i,:]),
            go.Scatter(y=new_C[i,:])],
            #Added 4 to the end indices since there 4 traces after (E,H,R,T)
            traces = [x for x in range(end_indices,end_indices+7)] #This is done because the rectangles are traces (to have hover information)
            ) for i in range(new_Nt)] #Iterate from the 2nd row to Nt-th row (now up until 1/3rd of Nt since it is very expensive to compute)

#Add the button and slider
print("Adding the button and sliders...")
updatemenus = [dict(type='buttons',
                    buttons=[dict(label='Play',
                                method='animate',
                                args=[[f'{i}' for i in range(new_Nt)], 
                                        dict(frame=dict(duration=0, redraw=False), 
                                            transition=dict(duration=0),
                                            easing='linear',
                                            fromcurrent=True,
                                            mode='immediate'
                                                                )])],
                    direction= 'left', 
                    pad=dict(r= 10, t=85), 
                    showactive =True, x= 0.1, y= 0, xanchor= 'right', yanchor= 'top')
            ]

sliders = [{'yanchor': 'top',
            'xanchor': 'left', 
            'currentvalue': {'font': {'size': 16}, 'prefix': 'Iteration: ', 'visible': True, 'xanchor': 'right'},
            'transition': {'duration': 0.0, 'easing': 'linear-in-out'},
            'pad': {'b': 10, 't': 50}, 
            'len': 0.9, 'x': 0.1, 'y': 0, 
            'steps': [{'args': [[i], {'frame': {'duration': 0.0, 'easing': 'linear-in-out', 'redraw': False},
                                    'transition': {'duration': 0, 'easing': 'linear-in-out'}}], 
                    'label': f"{i}/{new_Nt}", 'method': 'animate'} for i in range(new_Nt)       
                    ]}]





fig.update(frames=frames),
fig.update_layout(updatemenus=updatemenus,
                sliders=sliders)
#fig.show('chrome') #in jupyter notebook
#fig.show('browser') # in browser (the former offline.plot)
print("Writing to html....")
fig.write_html(filename[:-5]+".html")



# Only plot the subdomains if the algorithm is fdtd-schwarz...
if algo == "fdtd-schwarz":
    subplot_titles = ['Whole Computational Domain Plot']

    # Initialize x-axis vector for the subdom plots
    _,subdom_z_shape = new_subdomain_data[0]['E'].shape
    subdom_z = np.linspace(1,subdom_z_shape,subdom_z_shape)
    spec_list = [[{"colspan": 2}, None]]
    # Initialize the subplot titles
    for subdom_index in range(num_subdom):
        subplot_titles.append(f"Subdomain {subdom_index} Plot")
        if subdom_index % 2 == 0:
            spec_list.append([{},{}])
    
    #Create a new figure
    subdom_fig = make_subplots(rows=(int(num_subdom/2)+1), cols=2,specs=spec_list,subplot_titles=subplot_titles)

    # Add the trace for the computational domain...
    #Add traces in subplot 1
    subdom_fig.add_trace(go.Scatter(
                    x= z,
                    y= new_E[0,:],
                    mode = 'lines',
                    legendgroup= 'Comp Domain E-Field',
                    hovertemplate="x: %{x} <br> y: %{y}",
                    #line_color= 'rgb(255, 79, 38)',
                    name= 'Comp Domain E-Field',
                    showlegend= True), row=1, col=1)
                    
    subdom_fig.add_trace(go.Scatter(
                    x= z,
                    y= new_H[0,:],
                    mode = 'lines',
                    hovertemplate="x: %{x} <br> y: %{y}",
                    legendgroup= 'Comp Domain H-Field',
                    #line_color= 'rgb(255, 79, 38)',
                    name= 'Comp Domain H-Field',
                    showlegend= True), row=1, col=1)



    #Add the traces of each subdomain...
    row_tracker = 2

    for subdom_index in range(num_subdom):
        
        # If the subdom index is divisible by 2 (it is placed in the left column)
        if subdom_index % 2 == 0:
            # Add trace for Electric Field
            subdom_fig.add_trace(go.Scatter(
                        x= subdom_z,
                        y= new_subdomain_data[subdom_index]['E'][0,:],
                        mode = 'lines',
                        legendgroup= f'Subdomain {subdom_index} E-Field',
                        hovertemplate="x: %{x} <br> y: %{y}",
                        #line_color= 'rgb(255, 79, 38)',
                        name= f'Subdomain {subdom_index} E-Field',
                        showlegend= True), row=row_tracker, col=1)

            # Add trace for Magnetic Field
            subdom_fig.add_trace(go.Scatter(
                        x= subdom_z,
                        y= new_subdomain_data[subdom_index]['H'][0,:],
                        mode = 'lines',
                        legendgroup= f'Subdomain {subdom_index} H-Field',
                        hovertemplate="x: %{x} <br> y: %{y}",
                        #line_color= 'rgb(255, 79, 38)',
                        name= f'Subdomain {subdom_index} H-Field',
                        showlegend= True), row=row_tracker, col=1)

        # If the subdom index has a remainder 1 (we need to place it on the right column)
        elif subdom_index % 2 == 1: 
            # Add trace for Electric Field
            subdom_fig.add_trace(go.Scatter(
                        x= subdom_z,
                        y= new_subdomain_data[subdom_index]['E'][0,:],
                        mode = 'lines',
                        legendgroup= f'Subdomain {subdom_index} E-Field',
                        hovertemplate="x: %{x} <br> y: %{y}",
                        #line_color= 'rgb(255, 79, 38)',
                        name= f'Subdomain {subdom_index} E-Field',
                        showlegend= True), row=row_tracker, col=2)

            # Add trace for Magnetic Field
            subdom_fig.add_trace(go.Scatter(
                        x= subdom_z,
                        y= new_subdomain_data[subdom_index]['H'][0,:],
                        mode = 'lines',
                        legendgroup= f'Subdomain {subdom_index} H-Field',
                        hovertemplate="x: %{x} <br> y: %{y}",
                        #line_color= 'rgb(255, 79, 38)',
                        name= f'Subdomain {subdom_index} H-Field',
                        showlegend= True), row=row_tracker, col=2)
            row_tracker += 1

    subdom_fig.update_layout(title_text=f"FDTD Simulation [Date: {sim_date} | Algorithm: {algo} | Multithreading: {multithreading} | Source Type: {source_type} | Boundary Condition: {boundary_cond} | Excitation Method: {excitation_method}]")
    #Adjust the axes of the three subplots
    subdom_fig.update_xaxes(title_text="Computational Domain (m)", row=1,col=1)
    subdom_fig.update_yaxes(title_text="Level",range=[-1,1],row=1,col=1)

    row_tracker = 2
    for subdom_index in range(num_subdom):
        if subdom_index % 2 == 0:
            subdom_fig.update_xaxes(title_text="Cells",row=row_tracker,col=1)
            subdom_fig.update_yaxes(title_text="Level",range=[-1,1],row=row_tracker,col=1)
        elif subdom_index % 2 == 1:
            subdom_fig.update_xaxes(title_text="Cells",row=row_tracker,col=2)
            subdom_fig.update_yaxes(title_text="Level",range=[-1,1],row=row_tracker,col=2)
            row_tracker += 1

    print("Adding the frames for the subdomain plots...")
    # Adding the frames for each traces
    subdom_frames = []

    for i in range(new_Nt):
        new_dict = {}
        new_data = []
        # Add the data for the E and H
        new_dict['name'] = i
        new_data.append(go.Scatter(y=new_E[i,:]))
        new_data.append(go.Scatter(y=new_H[i,:]))
        
        # Add the frames for each subdomain
        for subdom_index in range(num_subdom):
            new_data.append(go.Scatter(y=new_subdomain_data[subdom_index]['E'][i,:]))
            new_data.append(go.Scatter(y=new_subdomain_data[subdom_index]['H'][i,:]))

        new_dict['data'] = new_data

        subdom_frames.append(new_dict)

    #Add the button and slider
    print("Adding the button and sliders...")
    subdom_updatemenus = [dict(type='buttons',
                        buttons=[dict(label='Play',
                                    method='animate',
                                    args=[[f'{i}' for i in range(new_Nt)], 
                                            dict(frame=dict(duration=0, redraw=False), 
                                                transition=dict(duration=0),
                                                easing='linear',
                                                fromcurrent=True,
                                                mode='immediate'
                                                                    )])],
                        direction= 'left', 
                        pad=dict(r= 10, t=85), 
                        showactive =True, x= 0.1, y= 0, xanchor= 'right', yanchor= 'top')
                ]

    subdom_sliders = [{'yanchor': 'top',
                'xanchor': 'left', 
                'currentvalue': {'font': {'size': 16}, 'prefix': 'Iteration: ', 'visible': True, 'xanchor': 'right'},
                'transition': {'duration': 0.0, 'easing': 'linear-in-out'},
                'pad': {'b': 10, 't': 50}, 
                'len': 0.9, 'x': 0.1, 'y': 0, 
                'steps': [{'args': [[i], {'frame': {'duration': 0.0, 'easing': 'linear-in-out', 'redraw': False},
                                        'transition': {'duration': 0, 'easing': 'linear-in-out'}}], 
                        'label': f"{i}/{new_Nt}", 'method': 'animate'} for i in range(new_Nt)       
                        ]}]
    subdom_fig.update(frames=subdom_frames),
    subdom_fig.update_layout(updatemenus=subdom_updatemenus,
                    sliders=subdom_sliders)
    #fig.show('chrome') #in jupyter notebook
    #fig.show('browser') # in browser (the former offline.plot)
    print("Writing to html....")
    subdom_fig.write_html(filename[:-5]+f"_{num_subdom}_subdomains"+".html")