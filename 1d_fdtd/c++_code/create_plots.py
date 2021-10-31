import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import csv
import h5py
from datetime import date, datetime
import sys







#Plotting functions for csv files
def plot_source_csv(fname,row,output_dir="",source_name="source"):
    plt.ion()
    fig  = plt.figure(1,[10,6])
    ax = fig.add_subplot(111)
    plt.ylabel("Level")
    x = np.zeros([3,row])




    with open(fname) as csv_f:
        csv_reader = csv.reader(csv_f,delimiter=',')
        i = 0
        for row in csv_reader:
            x[i,:] = row
            i += 1


    ax.plot(x[1,:])
    ax.plot(x[2,:])
    date_str = datetime.today().strftime('%Y-%m-%d')
    
    plt.show()
    plt.savefig(output_dir + date_str + '_' + source_name + '.png')


def plot_fields_csv(field = "", fname=["","",""],row = 0,col=0,col_fft = 0,output_dir=""):
    efield = np.zeros([row,col])
    hfield = np.zeros([row,col])
    refl = np.zeros([row,col_fft])
    trans = np.zeros([row,col_fft])
    con = np.zeros([row,col_fft])

    
    iteration = 0
    #Read efield
    with open(fname[1]) as csv_f:
        csv_reader = csv.reader(csv_f,delimiter=',')
        i = 0
        for row in csv_reader:
            efield[i,:] = row
            i += 1
            iteration += 1


    #Read hfield
    with open(fname[2]) as csv_f:
        csv_reader = csv.reader(csv_f,delimiter=',')
        i = 0
        for row in csv_reader:
            hfield[i,:] = row
            i += 1

    #Read refl
    with open(fname[3]) as csv_f:
        csv_reader = csv.reader(csv_f,delimiter=',')
        i = 0
        for row in csv_reader:
            refl[i,:] = row
            i += 1

    #Read trans
    with open(fname[4]) as csv_f:
        csv_reader = csv.reader(csv_f,delimiter=',')
        i = 0
        for row in csv_reader:
            trans[i,:] = row
            if row == 0:
                print(row.dtype)
            i += 1

    #Read con_of_energy
    with open(fname[5]) as csv_f:
        csv_reader = csv.reader(csv_f,delimiter=',')
        i = 0
        for row in csv_reader:
            con[i,:] = row
            i += 1

    plt.ion()
    fig, ax = plt.subplots(2)
    
    
    lineE, = ax[0].plot(efield[0,:])
    lineH, = ax[0].plot(hfield[0,:])

    lineR, = ax[1].plot(refl[0,:])
    lineT, = ax[1].plot(trans[0,:])
    lineC, = ax[1].plot(con[0,:])
    ax[0].set_ylim([-2,2])
    ax[1].set_ylim([-10,10])
   
    np.nan_to_num(efield)
    np.nan_to_num(hfield)
    np.nan_to_num(refl)
    np.nan_to_num(trans)
    np.nan_to_num(con)

    for k in range(0,iteration):
        ax[0].set(xlabel="Cells",ylabel="Levels",title = f"FDTD Simulation Iteration: {k}/{iteration}")
        ax[1].set(xlabel="Frequency",ylabel="Levels",title="FFT ")
        #Only for ax[0]
        #ax[0].set_ylim([min(hfield[k,:]),max(efield[k,:])])
        ax[0].legend(handles = [lineE,lineH],labels=["Electric Field","Magnetic Field"])
        lineE.set_ydata(efield[k,:])
        lineH.set_ydata(hfield[k,:])

        #For ax[1]
        #ax[1].set_ylim([min(refl[k,:]),max(con[k,:])])
        ax[1].legend(handles= [lineR,lineT,lineC],labels=["Reflectance","Transmittance","Conservation of Energy"])
        lineR.set_ydata(refl[k,:])
        lineT.set_ydata(trans[k,:])
        lineC.set_ydata(con[k,:])


        plot_name = "E_H_FFT_images_{num:07d}.png".format(num=k)
        plt.savefig(output_dir + plot_name)
        fig.canvas.draw()
        fig.canvas.flush_events()



#Plotting functions for npy files
def plot_fields_npy(fields = None,num_x = 0, num_y = 0, x_fft = 0,output_dir=""):
    plt.ion()
    fig, ax = plt.subplots(2)
    
    
    lineE, = ax[0].plot(fields[0,:,0])
    lineH, = ax[0].plot(fields[0,:,1])

    lineR, = ax[1].plot(fields[0,:,2])
    lineT, = ax[1].plot(fields[0,:,3])
    lineC, = ax[1].plot(fields[0,:,4])
    #ax[0].set_ylim([-2,2])
    #ax[1].set_ylim([-5,5])
    _,row,_ = fields.shape
    print(row)
    for k in range(0,row):
        ax[0].set(xlabel="Cells",ylabel="Levels",title = f"FDTD Simulation (Gaussian, PABC, Soft) Iteration: {k}/{row}")
        ax[1].set(xlabel="Frequency",ylabel="Levels",title="FFT ")
        #Only for ax[0]
        ax[0].set_ylim([min(fields[k,:,1]),max(fields[k,:,0])])
        ax[0].legend(handles = [lineE,lineH],labels=["Electric Field","Magnetic Field"])
        lineE.set_ydata(fields[k,:,0])
        lineH.set_ydata(fields[k,:,1])
        
        #For ax[1]
        ax[1].set_ylim([min(fields[k,:,2]),max(fields[k,:,4])])
        ax[1].legend(handles= [lineR,lineT],labels=["Reflectance","Transmittance"])
        lineR.set_ydata(fields[k,:,2])
        lineT.set_ydata(fields[k,:,3])
        lineC.set_ydata(fields[k,:,4])


        plot_name = "E_H_FFT_images_{num:07d}.png".format(num=k)
        plt.savefig(output_dir + plot_name)
        fig.canvas.draw()
        fig.canvas.flush_events()

def plot_source_npy(source_data = None,num_x = 0, num_y = 0, x_fft = 0,output_dir = "",source_name = ""):
    plt.ion()
    fig  = plt.figure(1,[10,6])
    ax = fig.add_subplot(111)
    plt.ylabel("Level")
    plt.xlabel("Time")
    plt.title("FDTD Source Excitation")
    ax.plot(source_data[:,1],source_data[:,2])  
    ax.plot(source_data[:,1],source_data[:,3])
    date_str = datetime.today().strftime('%Y-%m-%d')

    plt.savefig(output_dir + date_str + '_' + source_name + '.png')

#Converter function from numpy ndarray to str object
def np_to_str(input):
    return np.array_str(np.array(input).astype(str))

#Plotting function for hdf5 files
def plot_hdf5(data=None,save_plots=True,output_dir=""):

    #Read pre-processing data
    Nt = np.array(data['Total number of time iteration (Nt)'])
    source_type = np_to_str(data["Source Type"])
    boundary_cond = np_to_str(data['Boundary Condition'])
    excitation_method = np_to_str(data['Source Excitation Method'])
    z = np.array(data['Computational domain z (vector)'])
    source_title = f"Source Excitation ({source_type})"
    sim_title = "1D FDTD Simulation [SRC: " + source_type + " |BC: " + boundary_cond + " |EX: " + excitation_method + "]" 
    date_str = datetime.today().strftime('%Y-%m-%d')

    #Read source data
    Esrc = np.array(data['Esrc'])
    Hsrc = np.array(data['Hsrc'])
    t_src = np.array(data['t'])

    #Read simulation data
    E = np.array(data['E'])
    H = np.array(data['H'])
    Reflectance = np.array(data['Reflectance'])
    f = np.array(data['FFT Frequency Range'])
    Transmittance = np.array(data['Transmittance'])
    Con_of_Energy = np.array(data['Conservation of Energy'])

    #Plot source
    plt.figure(1)
    plt.plot(t_src,Esrc,linewidth=0.5)
    plt.plot(t_src,Hsrc,linewidth=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(source_title)
    #plt.show()
    if save_plots == True:
        plt.savefig(output_dir + date_str + '_' + source_type + '.png')

    #Find min and max of data to properly set ylim of plots
    #For simulation data
    sim_max = []
    sim_max.append(np.amax(E))
    sim_max.append(np.amax(H))

    sim_min = []
    sim_min.append(np.amin(E))
    sim_min.append(np.amin(H))

    #For FFT plot
    fft_max = []
    fft_max.append(np.amax(Reflectance))
    fft_max.append(np.amax(Transmittance))
    fft_max.append(np.amax(Con_of_Energy))

    fft_min = []
    fft_min.append(np.amin(Reflectance))
    fft_min.append(np.amin(Transmittance))
    fft_min.append(np.amin(Con_of_Energy))
    

    #Plot simulation data
    fig,ax = plt.subplots(2,figsize=(16, 9), dpi=(1920/16))
    lineE, = ax[0].plot(z,E[0,:],linewidth=0.5)
    lineH, = ax[0].plot(z,H[0,:],linewidth=0.5)
    ax[0].set(xlabel="Cells",ylabel="Levels",title = sim_title)
    ax[0].legend(handles = [lineE,lineH],labels=["Electric Field","Magnetic Field"])
    ax[0].set_ylim(np.amin(sim_min),np.amax(sim_max))

    lineR, = ax[1].plot(f,Reflectance[0,:],linewidth=0.5)
    lineT, = ax[1].plot(f,Transmittance[0,:],linewidth=0.5)
    lineC, = ax[1].plot(f,Con_of_Energy[0,:],linewidth=0.5,linestyle='dashed')
    ax[1].set(xlabel="Frequency (Hz) ",ylabel="Levels",title="FFT Response of Simulation")
    ax[1].legend(handles= [lineR,lineT,lineC],labels=["Reflectance","Transmittance","Conservation of Energy"])
    #ax[1].set_ylim(np.amin(fft_min),np.amax(fft_max))
    

    fargs = (E,H,Reflectance,Transmittance,Con_of_Energy,lineE,lineH,lineR,lineT,lineC,ax,sim_title,Nt)
    print("Start animation....")
    ani = animation.FuncAnimation(fig=fig,func=update_plot,frames=int(Nt),fargs=fargs)
    if save_plots == True:
        video_filename = output_dir + date_str + '_' + 'video_plot.mp4'
        writer = animation.FFMpegWriter(fps=200, metadata=dict(artist='Marvin',year=datetime.now().year), bitrate=8000)
        ani.save(video_filename, writer=writer)
        print("Finished saving animation....")
        
    print("Showing Animation.....s")
    plt.show()


def update_plot(i,E=0,H=0,R=0,T=0,C=0,l_E=None,l_H=None,l_R=None,l_T=None,l_C=None,ax=None,sim_title="",Nt=0):
    sim_title = sim_title + f" ({i}/{Nt})"
    fft_max = []
    fft_max.append(np.amax(R[i,:]))
    fft_max.append(np.amax(T[i,:]))
    fft_max.append(np.amax(C[i,:]))

    fft_min = []
    fft_min.append(np.amin(R[i,:]))
    fft_min.append(np.amin(T[i,:]))
    fft_min.append(np.amin(C[i,:]))
    ax[0].set(title=sim_title)
    ax[1].set_ylim(np.amin(fft_min),np.amax(fft_max))

    l_E.set_ydata(E[i,:])
    l_H.set_ydata(H[i,:])
    l_R.set_ydata(R[i,:])
    l_T.set_ydata(T[i,:])
    l_C.set_ydata(C[i,:])
    print(f"Iteration {i}/{Nt}")



def main():
    #Get the current date
    
    date_str = datetime.today().strftime('%Y-%m-%d')
    #Get the input arguments; Format: (1) input directory (2) Custom file name (3) Output file type (4)Source Name (for image) (5) Output directory (for plots)

    name = sys.argv[2]
    input_dir = sys.argv[1]
    type = sys.argv[3]
    source_name = sys.argv[4]
    output_dir = sys.argv[5]

    if type == 'csv':
        row = 5744
        col = 1336
        col_fft = 1336
        file_names = ["source.csv", "e_field.csv","h_field.csv","refl.csv","trans.csv","refl_trans.csv"]


        for i in range(len(file_names)):
            file_names[i] = input_dir + date_str + "_" + name + "_" + file_names[i]


        print(file_names)
        plot_source_csv(file_names[0],row,output_dir,source_name)
        plot_fields_csv(field="1D FDTD Simulation", fname = file_names,row=row,col=col,col_fft=col_fft,output_dir=output_dir)
    elif type == 'npy':
        npy_filename = input_dir + date_str + "_" + name + ".npy"
        data = np.load(input_dir+npy_filename)
        print(data.shape)
        plot_source_npy(data[:,:,0],num_x = data[0,2,0],output_dir=output_dir,source_name=source_name)

        plot_fields_npy(fields = data[:,:,1:], num_x = data[0,1,0], num_y = data[0,2,0], x_fft = data[0,0,0],output_dir=output_dir)

    elif type == 'hdf5':
        hdf5_filename = input_dir + date_str + "_" + name + ".hdf5"
        data = h5py.File(hdf5_filename,'r')
        plot_hdf5(data,save_plots=True,output_dir=output_dir)




if __name__ == "__main__":
    main()