#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
from mpl_toolkits.mplot3d import Axes3D                      
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time
import sys
import os.path
from matplotlib import cm


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_limits[0], z_limits[0] + 2.0*plot_radius])

if len(sys.argv)==2: 
 
    column_file = sys.argv[1]

else:

    print('Please provide the following argument:\n')
    print('1) File name (.col)\n')
    sys.exit()

umbrella_file = column_file.replace('.col','_std.p_2d')

print('umbrella file',umbrella_file)
data = np.loadtxt(umbrella_file,skiprows=0)

x = data[:,0]
y = data[:,1]
h = data[:,2]
u = data[:,3]
v = data[:,4]

mag_vel = np.sqrt(u**2+v**2)

x0_idx = np.asarray((np.where(data[:,0]==data[0,0])))

ny_cells = x0_idx[0,1]
nx_cells = data.shape[0] / ny_cells

nx_cells = nx_cells.astype(int)
ny_cells = ny_cells.astype(int)

X_cent = x.reshape((nx_cells,ny_cells))
Y_cent = y.reshape((nx_cells,ny_cells))
H_cent = h.reshape((nx_cells,ny_cells))
MAG_VEL_cent = mag_vel.reshape((nx_cells,ny_cells))


dx = X_cent[0,1]-X_cent[0,0]
dy = Y_cent[1,0]-Y_cent[0,0]

X_uni = np.unique(X_cent)
X_uni = X_uni - 0.5*dx
X_uni = np.append(X_uni,X_uni[-1]+dx)


Y_uni = np.unique(Y_cent)
Y_uni = Y_uni - 0.5*dy
Y_uni = np.append(Y_uni,Y_uni[-1]+dy)

X_stag, Y_stag = np.meshgrid(X_uni,Y_uni)

# idx = np.ma.masked_where(H_cent<=1.e-2,H_cent)
# H_cent[np.where(np.ma.getmask(idx)==True)] = np.nan

# create a figure for the plot
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')


col = np.loadtxt(column_file,skiprows=1)
zc = col[:,0]
rc = col[:,1]
xc = col[:,2]
yc = col[:,3]
rhom = col[:,4]

f = open(column_file)
line = f.readline()
f.close()
line_str = line.split()

i=0
for str in line_str:
    if "atm.rho" in str:
       idx = i
    i+=1 

rhoa = col[:,idx]

delta_rho = rhom - rhoa
below_nbl = np.argwhere(delta_rho<0.0)
ncol = below_nbl[-1]
z_nbl = zc[ncol]

v = np.linspace(zc[0],z_nbl,20)
vidx = np.searchsorted(zc,v)

theta = np.linspace(0, 2 * np.pi, 201)

for i in vidx:

    x = xc[i] + rc[i]*np.cos(theta)
    y = yc[i] + rc[i]*np.sin(theta)
    z = zc[i] + 0.0*np.sin(theta)

    ax2.plot(x,y,z,'b')


levels = zc[ncol]+np.linspace(10.0,H_cent.max(),20)
cset = ax2.contour(X_cent, Y_cent, H_cent+zc[ncol],levels=levels, cmap=cm.coolwarm)

dat0 = cset.allsegs[0][0]
plt.xlim(np.min(dat0[:,0]),np.max(dat0[:,0]))
plt.ylim(np.min(dat0[:,1]),np.max(dat0[:,1]))


ax2.clabel(cset, fontsize=9, inline=1)

"""
max_range = np.array([X_cent.max()-X_cent.min(), Y_cent.max()-Y_cent.min(), H_cent.max()+zc[ncol]]).max() / 2.0

mid_x = (X_cent.max()+X_cent.min()) * 0.5
mid_y = (Y_cent.max()+Y_cent.min()) * 0.5
mid_z = (H_cent.max()+zc[ncol]) * 0.5
ax2.set_xlim(mid_x - max_range, mid_x + max_range)
ax2.set_ylim(mid_y - max_range, mid_y + max_range)
ax2.set_zlim(0, 2.0*max_range)
"""

set_axes_equal(ax2)

plt.show()    


