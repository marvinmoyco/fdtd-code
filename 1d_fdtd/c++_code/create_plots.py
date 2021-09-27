import matplotlib.pyplot as plt
import numpy as np
import csv

plt.ion()
fig  = plt.figure(1,[10,6])
ax = fig.add_subplot(111)
plt.ylabel("Level")
x = np.zeros([3,2113])


with open("./csv/source.csv") as csv_f:
    csv_reader = csv.reader(csv_f,delimiter=',')
    i = 0
    for row in csv_reader:
        x[i,:] = row
        i += 1

ax.plot(x[0,:],x[1,:])
ax.plot(x[0,:],x[2,:])

plt.savefig("source.jpeg")