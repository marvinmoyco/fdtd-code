# Main Reference: https://chart-studio.plotly.com/~empet/15243/animating-traces-in-subplotsbr/#/


from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
import h5py

#Converter function from numpy ndarray to str object
def np_to_str(input):
    return np.array_str(np.array(input).astype(str))

#Read Electric field
data = h5py.File('./outputs/2021-11-11_fft.hdf5',mode='r')
print(data.keys())
z = data['Computational domain z (vector)'][:]
E = data['E'][:,:]
H = data['H'][:,:]
R = data['Reflectance'][:,:]
T = data['Transmittance'][:,:]
sim_date = np_to_str((data['Date of Simulation']))
source_type = np_to_str(data['Source Type'])
n = data['Refractive Index (vector)'][:]
Nt = int(np.array(data['Total number of time iteration (Nt)']))

dz = float(np.array(data['Cell size (dz)']))
spacer = int(np.array(data['Amount of spacing (number of cells)']))*dz
#print(z.shape,E.shape,Nt)
input_layer = data['Input Layer size'][:]
print(input_layer,spacer)

#Create subplotsd
fig = make_subplots(rows=2, cols=1, subplot_titles = ('FDTD Simulation', 'Frequency Response (Reflectance and Transmittance)'))

#Get the total maximum for E and H fields
max_point = np.max([np.amax(E),np.amax(H)])
min_point = np.min([np.amin(E),np.amin(H)])
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

#Add traces in subplot 1
fig.add_trace(go.Scatter(
                x= z,
                y= E[0,:],
                mode = 'lines',
                hoverinfo='name',
                legendgroup= 'Electric Field',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Electric Field',
                showlegend= True), row=1, col=1)
                
fig.add_trace(go.Scatter(
                x= z,
                y= H[0,:],
                mode = 'lines',
                hoverinfo='name',
                legendgroup= 'Magnetic Field',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Magnetic Field',
                showlegend= True), row=1, col=1)



#Add traces in subplot 2
fig.add_trace(go.Scatter(
                x= z,
                y= R[0,:],
                mode = 'lines',
                hoverinfo='name',
                legendgroup= 'Reflectance',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Reflectance',
                showlegend= True), row=2, col=1)

fig.add_trace(go.Scatter(
                x= z,
                y= T[0,:],
                mode = 'lines',
                hoverinfo='name',
                legendgroup= 'Transmittance',
                #line_color= 'rgb(255, 79, 38)',
                name= 'Transmittance',
                showlegend= True), row=2, col=1)


print(len(fig.data))





fig.update_layout(width=1920, height=1080,title_text=f"FDTD Simulation [Date: {sim_date} | Source Type: {source_type}]")
#Adjust the axes of the two subplots
fig.update_xaxes(title_text="Computational Domain (m)", row=1,col=1)
fig.update_xaxes(title_text="Frequency (Hz)",row=2,col=1)

fig.update_yaxes(title_text="Level",range=[-1,1],row=1,col=1)
fig.update_yaxes(title_text="Magnitude",range=[0,1],row=2,col=1)

end_indices = len(input_layer) 
#Add the frames of each traces
frames=[dict(name=i,
             data=[go.Scatter(y=E[i,:]),#update the E-Field
             go.Scatter(y=H[i,:]),#update the H-Field
             go.Scatter(y=R[i,:]),
             go.Scatter(y=T[i,:])],
             #Added 4 to the end indices since there 4 traces after (E,H,R,T)
             traces = [x for x in range(end_indices,end_indices+4)] #This is done because the rectangles are traces (to have hover information)
             ) for i in range(Nt)] #Iterate from the 2nd row to Nt-th row

#Add the button and slider
updatemenus = [dict(type='buttons',
                    buttons=[dict(label='Play',
                                  method='animate',
                                  args=[[f'{i}' for i in range(Nt)], 
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
                       'label': f"{i}/{Nt}", 'method': 'animate'} for i in range(Nt)       
                    ]}]





fig.update(frames=frames),
fig.update_layout(updatemenus=updatemenus,
                  sliders=sliders)
fig.show() #in jupyter notebook
#fig.show('browser') # in browser (the former offline.plot)
fig.write_html('sample.html')
