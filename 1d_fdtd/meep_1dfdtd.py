import numpy as np
from numpy import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import animation


k_total = 500

ex = []
hy = []

def initialize():
	for k in range(0,k_total):
		ex.append(0)
		hy.append(0)
	return ex, hy

ex, hy = initialize()

k_center = k_total/2
t0 = 40.0
spread = 12
T = 0
nsteps = 1000

exLow2 = 0
exLow1 = 0
exHigh1 = 0
exHigh2 = 0

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, k_total), ylim=(-1, 2))
# ax.autoscale(True,'both')
point, = ax.plot([], [], lw=2)


# initialization function: plot the background of each frame
def init():
    point.set_data([], [])
    return point,

def animate(i):
	global k_center, t0, spread, T, nsteps, ex, hy, exLow2,exLow1,exHigh2,exHigh1
	for i in arange(i):

		T = T + 1

		for k in range(1,k_total):
			ex[k] = ex[k] + 0.5*(hy[k-1] - hy[k])

		pulse = np.exp(-0.5*(np.power((t0-T)/spread,2.)))
		ex[k_center] = pulse
		# ex[k_center-20] = pulse
		# ex[k_center+20] = pulse

		# absorbing boundary condition
		ex[0] = exLow2
		exLow2 = exLow1
		exLow1 = ex[1]

		ex[k_total-1] = exHigh2
		exHigh2 = exHigh1
		exHigh1 = ex[k_total-2]

		for k in range(0,k_total-1):
			hy[k] = hy[k] + 0.5*(ex[k] - ex[k+1])


	point.set_data(arange(k_total),ex)
	return point,


# for i in arange(nsteps):
# 	T = T + 1

# 	for k in range(1,k_total):
# 		ex[k] = ex[k] + 0.5*(hy[k-1] - hy[k])

# 	pulse = np.exp(-0.5*(np.power((t0-T)/spread,2.)))
# 	ex[k_center] = pulse

# 	for k in range(0,k_total-1):
# 		hy[k] = hy[k] + 0.5*(ex[k] - ex[k+1])

# fig = plt.figure()
# ax = fig.add_subplot(212)
# ax2 = fig.add_subplot(211)
# ax.autoscale(True,'both')
# ax2.autoscale(True,'both')
# ax.plot(arange(k_total),ex)
# ax2.plot(arange(k_total),hy)

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=100, blit=True)

# anim.save('1d_fdtd3.gif', writer='imagemagick', fps=30)

plt.show()


# material_thickness = 0.0625

# # x_step = 0.001
# theta_step = 0.001
# # x = arange(0.001,.25,x_step)
# # ion()
# l1 = 0.5
# l2 = 2.5
# # theta1 = np.arcsin(x/l1)
# theta1 = arange(0,pi/4,theta_step)
# theta2 = np.arcsin(l1*np.sin(theta1)/l2)
# y = l1*cos(theta1)+l2*cos(theta2)
# print(y)
# dthdy = theta_step/np.gradient(y)
# print(dthdy)
# # x = 0.235
# # y_start = y[len(y)-1] + material_thickness
# y = -(y-y[0])
# th_min = np.interp(material_thickness,y,theta1)
# dthdy_min = np.interp(material_thickness,y,dthdy)
# print(dthdy_min)

# fig = plt.figure()
# ax = fig.add_subplot(212, xlim=(0, np.rad2deg(0.25)), ylim=(-40, 10))
# ax2 = fig.add_subplot(211, xlim=(0, np.rad2deg(0.25)), ylim=(0, 0.125))
# # ax.autoscale(True,'y')
# ax.plot(np.rad2deg(theta1),dthdy,'green')
# # ax.plot(x,y*100-min(y)*100,'red')
# ax.plot(np.rad2deg(th_min),dthdy_min,'ko',lw=10)
# ax.plot([np.rad2deg(th_min), np.rad2deg(th_min)],[10,-40],'k--')
# s = str("Mechanical Advantage = %.3f\n@ theta = %.1f degrees (%.3f radians)" % (dthdy_min,rad2deg(th_min),th_min,))
# ax.annotate(s,xy=(np.rad2deg(th_min), dthdy_min), arrowprops=dict(arrowstyle='->'), xytext=(np.rad2deg(0.0125), 0))
# ax.set_xlabel('Input linkage angle (degrees)')
# ax.set_ylabel('Mechanical Advantage')

# ax2.plot(np.rad2deg(theta1),y,'red')
# ax2.plot(np.rad2deg(th_min),material_thickness,'ko',lw=10)
# ax2.plot([0,np.rad2deg(th_min),np.rad2deg(th_min)],[material_thickness, material_thickness,0],'k--')
# ax2.set_ylabel('Material Thickness (inches)')
# plt.show()


# xx = []
# yy = []

# plt.xkcd()

# First set up the figure, the axis, and the plot element we want to animate
# fig = plt.figure()
# ax = fig.add_subplot(111, xlim=(0, 1.25), ylim=(-10, 100))
# point, = ax.plot([], [], lw=5)


 
# # initialization function: plot the background of each frame
# def init():
#     point.set_data([], [])
#     return point,
 
# # animation function.  This is called sequentially
# def animate(i):
# 	global xx, yy, y, x, cor, vx, vy, tmax, g
# 	if (y > 0):
# 		vy -= g
# 	else:
# 		vy = -vy*cor
# 	y += vy
# 	x += vx
# 	if (y<0):
# 		y=0
# 	xx.append(x)
# 	yy.append(y)
# 	print "(",x,", ",y,")"
# 	point.set_data(xx, yy)
# 	return point,
 
# # call the animator.  blit=True means only re-draw the parts that have changed.
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=125, interval=20, blit=True)
 
# # this is how you save your animation to file:
# #anim.save('animation.gif', writer='imagemagick_file', fps=30)
# anim.save('animation.gif', writer='imagemagick', fps=30)
 
# plt.show()