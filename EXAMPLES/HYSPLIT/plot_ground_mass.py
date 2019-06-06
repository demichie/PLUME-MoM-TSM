from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")

import Tkinter, tkFileDialog
import numpy as np
import sys
from haversine import haversine
import os

from input_file import *
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as ticker

import pandas as pd
import salem
from salem import get_demo_file, DataLevels, GoogleVisibleMap, Map

import easygui


fname = 'CON2ASC.GROUND'
with open(fname) as f:
    lines = f.read().splitlines()

for filename in lines:
    exists = os.path.isfile(filename.strip())
    if exists:
        
        os.rename(filename.strip(),filename.strip()+'.gnd')

filename = easygui.fileopenbox( filetypes=['*.gnd'])

# filename = 'cdumpcum_part_sunset_may_1e8_131_1000'


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

GROUND=[]
AIR=[]

#print npart, n_levels, H_LEVELS
print ' '
print '*** MASS ON THE GROUND ***'
print ' '
# Check mass deposited on the ground
print filename

f = open(filename)

und_where = ( [pos for pos, char in enumerate(filename) if char == '_'])
dot_where = ( [pos for pos, char in enumerate(filename) if char == '.'])

time = filename[und_where[-1]+1:dot_where[0]]
day = filename[und_where[-2]+1:und_where[-1]]

# time = filename.strip()[-8:-4]
# day = filename.strip()[-12:-9]
           
print ' ---> day and time ',day,' ',time,' '

count = 0
header = []
with open(filename) as f:
    lines = f.readlines()

    for line in lines:
        if line[0].isdigit() == False:
            header = header + line.split()
            count = count + 1
        else:
            break

f.close()    
         
m = []        

for j in range(99):

    h_new = []

    to_find = 'CL'+str(int(j)).zfill(2)

    occurrence = 0
             
    for i in range(len(header)):

        if to_find in header[i]:

            occurrence = occurrence + 1

            h_new.append(int(header[i][4:]))


    if occurrence > 0 :

        h_old=np.asarray(h_new)

        h_old=h_old.reshape((-1,1))

        m.append([j,occurrence])

        n_height = h_old.shape[0]

    if occurrence == 0 :
 
        h = h_old
                                 
        break                 


m = np.asarray(m)
m=m.reshape((-1,2))
        
npart = m.shape[0]
n_levels = h.shape[0]
H_LEVELS = h

total_mass = 0

#print 'Number of particle classes :'npart
#print 'Heights :',H_LEVELS

a = np.loadtxt(filename, skiprows = int(count))
         
if a.shape[0] == 0 :
         
    print 'No mass deposited at ',time
        
else:

    a = np.asarray(a) 

    a = a.reshape((-1,(npart * n_levels + 4)))

    lat = a[:,2]

    lat = lat.reshape((-1,1))

    lon = a[:,3]

    lon = lon.reshape((-1,1))

    lon_unique = np.unique(lon)

    lat_unique = np.unique(lat)
    
    Lon,Lat = np.meshgrid(lon_unique,lat_unique) 

    lon_unique = lon_unique.reshape((-1,1))
    lat_unique = lat_unique.reshape((-1,1))

    Lon,Lat = np.meshgrid(lon_unique,lat_unique) 

    spacing_lat = np.absolute(lat_unique[0,0] - lat_unique[1,0])

    spacing_lon = np.absolute(lon_unique[0,0] - lon_unique[1,0])

    lon_stag = lon_unique-0.5*spacing_lon
    lon_stag = np.append(lon_stag,lon_unique[-1]+0.5*spacing_lon)
    
    lat_stag = lat_unique-0.5*spacing_lat
    lat_stag = np.append(lat_stag,lat_unique[-1]+0.5*spacing_lat)
    Lon_stag,Lat_stag = np.meshgrid(lon_stag,lat_stag) 

    dist = []

    loading2D = np.zeros((lat_unique.shape[0],lon_unique.shape[0],npart * n_levels))
    loading2D_sum = np.zeros((lat_unique.shape[0],lon_unique.shape[0]))
    
    for i in range(a.shape[0]):

        idx_lon = np.where(lon_unique==lon[i,0])
        idx_lat = np.where(lat_unique==lat[i,0])

        loading2D[idx_lat,idx_lon,:] = a[i,4:]
        
        point1 = ((lat[i,0]-spacing_lat/float(2)),lon[i,0])
        point2 = ((lat[i,0]+spacing_lat/float(2)),lon[i,0])
                  
        d_lat = haversine(point1, point2) 

        point3 = (lat[i,0],(lon[i,0]-spacing_lon/float(2)))
        point4 = (lat[i,0],(lon[i,0]+spacing_lon/float(2)))
                  
        d_lon = haversine(point3, point4)       
                 
        dist.append([d_lat*1000,d_lon*1000])

    dist = np.asarray(dist)
    dist = dist.reshape((-1,2))

    a = a[:,4:]

    a = a.reshape((-1,(npart * n_levels)))

    min_x = np.amin(lon_stag)
    max_x = np.amax(lon_stag)

    min_y = np.amin(lat_stag)
    max_y = np.amax(lat_stag)

    g = GoogleVisibleMap(x=[min_x, max_x], y=[min_y, max_y],
                         scale=2,  # scale is for more details
                         maptype='terrain')  # try out also: 'terrain'

    ggl_img = g.get_vardata()

    column = 0

    
    for i in range(npart):

            conc = a[:, column]

            conc = conc.reshape((-1,1))
            min_conc = np.amin(conc)
            max_conc = np.amax(conc)

            half_conc = 0.5 * ( min_conc + max_conc )

            half_int = np.floor(np.log10(half_conc))
            
            log_scale = np.logspace(half_int-6,half_int, num=7)
            lin_scale = np.linspace(half_conc, max_conc, num=8)

            levels = np.append(log_scale,lin_scale)
            
            mass_on_the_ground = np.sum(conc[:,0] * dist[:,0] * dist[:,1])
            total_mass = total_mass +  mass_on_the_ground 

            loading_i = np.maximum(loading2D[:,:,column],1.e-20)
            loading2D_sum += loading_i

            print 'Deposit class CL',str(i+1).zfill(2),' mass ','%.1e'%mass_on_the_ground,' kg'
            
            f = plt.figure(i)
            plt.rcParams["font.size"] = 8.0
            cmap = plt.get_cmap('Spectral')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            
            sm = Map(g.grid, factor=1, countries=False)
            sm.set_rgb(ggl_img)  # add the background rgb image
            sm.visualize()

            loading_i[np.log10(loading_i)<half_int-6] = np.nan
            x, y = sm.grid.transform(Lon_stag, Lat_stag)
    
            plt.pcolormesh(x,y,loading_i, cmap=cmap, norm=norm,alpha=0.50,zorder=1)
            plt.pcolormesh(Lon_stag, Lat_stag,loading_i, cmap=cmap, norm=norm,alpha=1.0)

            xvent, yvent = sm.grid.transform(vent_lon, vent_lat)
            plt.plot(xvent, yvent,"^m",markersize=1,zorder=2)

            plt.xlim(left=np.amin(x))
            plt.xlim(right=np.amax(x))
            plt.ylim(bottom=np.amax(y))
            plt.ylim(top=np.amin(y))
            plt.grid()
            plt.title('Class CL'+str(i+1).zfill(2)+' mass '+'%.1e'%mass_on_the_ground+' kg' )
            clb = plt.colorbar(format=ticker.FuncFormatter(fmt))
            clb.set_label('Loading (kg/m^2)', labelpad=-40, y=1.05, rotation=0)

            f.savefig('CL'+str(i+1)+'_'+day+'_'+time+'.pdf', bbox_inches='tight')
            column = column + n_levels

f = plt.figure(i+1)
plt.rcParams["font.size"] = 8.0

min_conc = np.amin(loading2D_sum)
max_conc = np.amax(loading2D_sum)

half_conc = 0.5 * ( min_conc + max_conc )

half_int = np.floor(np.log10(half_conc))
            
log_scale = np.logspace(half_int-6,half_int, num=7)
lin_scale = np.linspace(half_conc, max_conc, num=8)

levels = np.append(log_scale,lin_scale)
cmap = plt.get_cmap('Spectral')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                     
sm = Map(g.grid, factor=1, countries=False)
sm.set_rgb(ggl_img)  # add the background rgb image
sm.visualize()

loading2D_sum[np.log10(loading2D_sum)<half_int-6] = np.nan

x, y = sm.grid.transform(Lon_stag, Lat_stag)
plt.pcolormesh(x,y,loading2D_sum, cmap=cmap, norm=norm,alpha=0.50,zorder=1)
plt.pcolormesh(Lon_stag, Lat_stag,loading2D_sum, cmap=cmap, norm=norm,alpha=1.0)

xvent, yvent = sm.grid.transform(vent_lon, vent_lat)
plt.plot(xvent, yvent,"^m",markersize=1,zorder=2)

plt.grid()
plt.xlim(left=np.amin(x))
plt.xlim(right=np.amax(x))
plt.ylim(bottom=np.amax(y))
plt.ylim(top=np.amin(y))
plt.title('Total deposit')
clb = plt.colorbar(format=ticker.FuncFormatter(fmt))
clb.set_label('Loading (kg/m^2)', labelpad=-40, y=1.05, rotation=0)
f.savefig('CL_sum'+'_'+day+'_'+time+'.pdf', bbox_inches='tight')

# plt.show()        

            

