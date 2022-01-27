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
E = np.nan_to_num(data['output']['E'][:,:])
H = np.nan_to_num(data['output']['H'][:,:])
R = np.nan_to_num(data['output']['reflectance'][:,:])
T = np.nan_to_num(data['output']['transmittance'][:,:])
C = np.nan_to_num(data['output']['conservation_of_energy'][:,:])





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
#print(input_layer,spacer)
new_Nt = 0
#Adjust the frames here
for i in range(Nt):
    if i == 0:
        new_E = E[i,:]
        new_H = H[i,:]
        new_R = R[i,:]
        new_T = T[i,:]
        new_C = C[i,:]
        new_Nt += 1
    elif i % 10 == 0:
        new_E = np.vstack((new_E,E[i,:]))
        new_H = np.vstack((new_H,H[i,:]))
        new_R = np.vstack((new_R,R[i,:]))
        new_T = np.vstack((new_T,T[i,:]))
        new_C = np.vstack((new_C,C[i,:]))
        new_Nt += 1


print(f"Field shapes (after adjusting): E: {new_E.shape} | H: {new_H.shape} | R: {new_R.shape} | T: {new_T.shape} | C: {new_C.shape}")
print(f"New Nt (total number of frames to be rendered): {new_Nt}")

#Create subplotsd
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
                text= f'Layer {index + 1} \n Refractive index: {n[index]}',
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
fig.update_layout(title_text=f"FDTD Simulation [Date: {sim_date} | Source Type: {source_type} | Boundary Condition: {boundary_cond} | Excitation Method: {excitation_method}]")
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
