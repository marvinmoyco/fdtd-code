import matplotlib.pyplot as plt
import numpy as np
import csv
from datetime import datetime




#Plot source
def plot_source(fname,row):
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

    ax.plot(x[0,:],x[1,:])
    ax.plot(x[0,:],x[2,:])

    plt.savefig("./plots/source.jpeg")

def plot_field(field = "", fname=["","",""],row = 0,col=0):
    plt.ion()
    fig  = plt.figure(1,[10,6])
    ax = fig.add_subplot(111)
    plt.ylabel("Level")
    efield = np.zeros([row,col])
    hfield = np.zeros([row,col])
    plt.title(field)
    
    
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

    lineE, = ax.plot(efield[0,:])
    lineH, = ax.plot(hfield[0,:])

    for k in range(1,iteration):
        plt.ylim(-3,3)
        plt.legend(handles = [lineE,lineH],labels=["Electric Field","Magnetic Field"])
        lineE.set_ydata(efield[k,:])
        lineH.set_ydata(hfield[k,:])
        plot_name = "./plots/E_H_images_{num:07d}.jpeg".format(num=k)
        plt.savefig(plot_name)
        fig.canvas.draw()
        fig.canvas.flush_events()



def main():
    #Get the current date
    row = 2100
    col = 400
    date_str = datetime.today().strftime('%Y-%m-%d')
    name = "marvin"
    curr_dir = "./csv/"
    file_names = ["source.csv", "e_field.csv","h_field.csv"]

    for i in range(len(file_names)):
        file_names[i] = curr_dir + date_str + "_" + name + "_" + file_names[i]

    print(file_names)
    plot_source(file_names[0],row)
    plot_field(field="1D FDTD (PABC, Hard Source - Gaussian)", fname = file_names,row=row,col=col)


if __name__ == "__main__":
    main()