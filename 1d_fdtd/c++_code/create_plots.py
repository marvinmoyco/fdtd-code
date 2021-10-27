import matplotlib.pyplot as plt
import numpy as np
import csv
from datetime import datetime








#Plot source
def plot_source_csv(fname,row):
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


    plt.savefig("./plots/source.jpeg")


def plot_fields_csv(field = "", fname=["","",""],row = 0,col=0,col_fft = 0):
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
    ax[0].set_ylim([-1,1])
    ax[1].set_ylim([-5,5])
    
    for k in range(1300,iteration):
        ax[0].set(xlabel="Cells",ylabel="Levels",title = f"FDTD Simulation Iteration: {k}/{iteration}")
        ax[1].set(xlabel="Frequency",ylabel="Levels",title="FFT ")
        #Only for ax[0]
        #ax[0].set_ylim([-3,3])
        ax[0].legend(handles = [lineE,lineH],labels=["Electric Field","Magnetic Field"])
        lineE.set_ydata(efield[k,:])
        lineH.set_ydata(hfield[k,:])

        #For ax[1]
        #ax[1].set_ylim([0,1])
        ax[1].legend(handles= [lineR,lineT,lineC],labels=["Reflectance","Transmittance","Conservation of Energy"])
        lineR.set_ydata(refl[k,:])
        lineT.set_ydata(trans[k,:])
        lineC.set_ydata(con[k,:])


        plot_name = "./plots/E_H_FFT_images_{num:07d}.jpeg".format(num=k)
        plt.savefig(plot_name)
        fig.canvas.draw()
        fig.canvas.flush_events()


def plot_fields_npy(fields = None,num_x = 0, num_y = 0, x_fft = 0):
    plt.ion()
    fig, ax = plt.subplots(2)
    
    
    lineE, = ax[0].plot(fields[0,0,:])
    lineH, = ax[0].plot(fields[1,0,:])

    lineR, = ax[1].plot(fields[2,0,:])
    lineT, = ax[1].plot(fields[3,0,:])
    #lineC, = ax[1].plot(fields[0,:])
    ax[0].set_ylim([-2,2])
    #ax[1].set_ylim([-5,5])
    _,row,_ = fields.shape
    print(row)
    print(not fields[2,:].any())
    print(not fields[3,:].any())
    print(np.amax(fields[2,:]))
    print(np.amax(fields[3,:]))
    for k in range(1350,row):
        ax[0].set(xlabel="Cells",ylabel="Levels",title = f"FDTD Simulation (Gaussian, PABC, Soft) Iteration: {k}/{row}")
        ax[1].set(xlabel="Frequency",ylabel="Levels",title="FFT ")
        #Only for ax[0]
        #ax[0].set_ylim([-3,3])
        ax[0].legend(handles = [lineE,lineH],labels=["Electric Field","Magnetic Field"])
        lineE.set_ydata(fields[0,k,:])
        lineH.set_ydata(fields[1,k,:])
        
        #For ax[1]
        #ax[1].set_ylim([0,1])
        ax[1].legend(handles= [lineR,lineT],labels=["Reflectance","Transmittance"])
        lineR.set_ydata(fields[2,k,:])
        lineT.set_ydata(fields[3,k,:])
        #lineC.set_ydata(con[k,:])


        plot_name = "./plots/E_H_FFT_images_{num:07d}.jpeg".format(num=k)
        plt.savefig(plot_name)
        fig.canvas.draw()
        fig.canvas.flush_events()

def plot_source_npy(source_data = None,num_x = 0, num_y = 0, x_fft = 0):
    plt.ion()
    fig  = plt.figure(1,[10,6])
    ax = fig.add_subplot(111)
    plt.ylabel("Level")
    plt.xlabel("Time")
    plt.title("FDTD Source Excitation")
    ax.plot(source_data[:,1],source_data[:,2])  
    ax.plot(source_data[:,1],source_data[:,3])

    plt.savefig("./plots/source.jpeg")

def main():
    #Get the current date
    
    date_str = datetime.today().strftime('%Y-%m-%d')
    name = "marvin"
    curr_dir = "./input_output/"
    type = 'npy'

    if type == 'csv':
        row = 5744
        col = 1336
        col_fft = 1000
        file_names = ["source.csv", "e_field.csv","h_field.csv","refl.csv","trans.csv","refl_trans.csv"]


        for i in range(len(file_names)):
            file_names[i] = curr_dir + date_str + "_" + name + "_" + file_names[i]


        print(file_names)
        plot_source_csv(file_names[0],row)
        plot_fields_csv(field="1D FDTD (PABC, Hard Source - Gaussian)", fname = file_names,row=row,col=col,col_fft=col_fft)
    elif type == 'npy':
        npy_filename = 'output.npy'
        data = np.load(curr_dir+npy_filename)
        plot_source_npy(data[0,:,:],num_x = data[0,0,2])

        plot_fields_npy(fields = data[1:,:,:], num_x = data[0,0,1], num_y = data[0,0,2], x_fft = data[0,0,0])





if __name__ == "__main__":
    main()