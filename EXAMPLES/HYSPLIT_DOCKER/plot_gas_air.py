import numpy as np
import sys
from haversine import haversine
import os
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as ticker
import pandas as pd
import salem
from salem import get_demo_file, DataLevels, GoogleVisibleMap, Map
import warnings
warnings.simplefilter(action='ignore')

import streamlit as st

from input_file import *

def file_selector(folder_path='.', ext='gas'):
    filenames = os.listdir(folder_path)
    filelist = []
    for file in filenames:
        if file.endswith(".gas"):
            filelist.append(file)

    filelist.insert(0, '')    

    selected_filename = st.selectbox('Select a file', filelist, format_func=lambda x: 'Select a file' if x == '' else x)

    if selected_filename !='':
        st.write('You selected `%s`' % selected_filename)
    elif selected_filename =='':
        st.stop()
        st.warning('No option is selected')

    return os.path.join(folder_path,selected_filename)


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


# Read the list of air concentration files created by Hysplit
# utility con2asc
fname = 'CON2ASC.GAS'
with open(fname) as f:
    lines = f.read().splitlines()

for filename in lines:
    exists = os.path.isfile(filename.strip())
    if exists:
        
        # Add .air extension to distinguish the files from other types 
        os.rename(filename.strip(),filename.strip()+'.gas')

f.close()

filename = file_selector(ext='gas')

AIR=[]

#print ( npart, n_levels, H_LEVELS
print ( ' ' )
print ( '*** VOLCGAS MASS IN THE AIR ***' )
print ( ' ' )
print ( "Run name: ",runname )
print ( ' ' )

f = open(filename)

# find day and time from filename
und_where = ( [pos for pos, char in enumerate(filename) if char == '_'])
dot_where = ( [pos for pos, char in enumerate(filename) if char == '.'])

time = filename[und_where[-1]+1:dot_where[1]]
day = filename[und_where[-2]+1:und_where[-1]]
           
print ( ' ---> day and time ',day,' ',time,' ' )
print ( ' ' )

data = f.read()
first_line = data.split('\n', 1)[0]
        
line_split = first_line.split()

m = []        

total_mass = 0.0

for j in range(99):

    h_new = []

    to_find = 'GS'+str(int(j)).zfill(2)

    occurrence = 0
             
    for i in range(len(line_split)):

        if to_find in  line_split[i]:

            occurrence = occurrence + 1

            h_new.append(int(line_split[i][4:]))


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
        
ngas = m.shape[0]
n_levels = h.shape[0]
H_LEVELS = h

print ( 'Number of volcanic gasses : ',ngas )
print ( 'Number of atmospheric levels :',n_levels - 1 )

a = np.loadtxt(filename, skiprows = 1)
         
if a.shape[0] == 0 :
         
    print ( 'No volcanic gasses in the air at time ',time )
        
else:

    a = np.asarray(a) 

    a = a.reshape((-1,(ngas * n_levels + 4)))

    # extract and reshape grid latitude values (at pixel centers)
    lat = a[:,2]
    lat = lat.reshape((-1,1))
    lat_unique = np.unique(lat)
    lat_unique = lat_unique.reshape((-1,1))

    # extract and reshape grid longitude values (at pixel centers) 
    lon = a[:,3]
    lon = lon.reshape((-1,1))
    lon_unique = np.unique(lon)
    lon_unique = lon_unique.reshape((-1,1))

    # create two meshgrid 2D arrays 
    Lon,Lat = np.meshgrid(lon_unique,lat_unique) 

    # compute the grid spacing
    spacing_lat = np.absolute(lat_unique[0,0] - lat_unique[1,0])
    spacing_lon = np.absolute(lon_unique[0,0] - lon_unique[1,0])

    # create staggered (at pixel edges) longitude values
    lon_stag = lon_unique-0.5*spacing_lon
    lon_stag = np.append(lon_stag,lon_unique[-1]+0.5*spacing_lon)
    
    # create staggered (at pixel edges) latitude values
    lat_stag = lat_unique-0.5*spacing_lat
    lat_stag = np.append(lat_stag,lat_unique[-1]+0.5*spacing_lat)

    # create two meshgrid staggered 2D arrays 
    Lon_stag,Lat_stag = np.meshgrid(lon_stag,lat_stag) 

    # ll_utm = utm.from_latlon(lat_stag[0],lon_stag[0])

    dist = []

    # allocate leading array for all the classes
    loading3D = np.zeros((lat_unique.shape[0],lon_unique.shape[0],ngas * n_levels))

    # allocate leading array for total loading
    loading3D_sum_kg = np.zeros((lat_unique.shape[0],lon_unique.shape[0]))
    
    for i in range(a.shape[0]):

        idx_lon = np.where(lon_unique==lon[i,0])
        idx_lat = np.where(lat_unique==lat[i,0])

        loading3D[idx_lat,idx_lon,:] = a[i,4:]
        
        point1 = ((lat[i,0]-spacing_lat/float(2)),lon[i,0])
        point2 = ((lat[i,0]+spacing_lat/float(2)),lon[i,0])
                  
        d_lat = haversine(point1, point2) 

        point3 = (lat[i,0],(lon[i,0]-spacing_lon/float(2)))
        point4 = (lat[i,0],(lon[i,0]+spacing_lon/float(2)))
                  
        d_lon = haversine(point3, point4)       
                 
        dist.append([d_lat*1000,d_lon*1000])

    dist = np.asarray(dist)
    dist = dist.reshape((-1,2))

    area = dist[:,0] * dist[:,1]

    # array for pixel area
    pixel_area = area.reshape((lat_unique.shape[0],lon_unique.shape[0]))

    a = a[:,4:]

    a = a.reshape((-1,(ngas * n_levels)))

    min_x = np.amin(lon_stag)
    max_x = np.amax(lon_stag)

    min_y = np.amin(lat_stag)
    max_y = np.amax(lat_stag)

    g = GoogleVisibleMap(x=[min_x, max_x], y=[min_y, max_y],
                         scale=2,  # scale is for more details
                         maptype='terrain')  # try out also: 'terrain,hybrid'

    ggl_img = g.get_vardata()

    column = 0

    header = "ncols     %s\n" % loading3D.shape[1]
    header += "nrows    %s\n" % loading3D.shape[0]
    header += "xllcorner " + str(lon_stag[0]) +"\n"
    header += "yllcorner " + str(lat_stag[1]) +"\n"
    header += "cellsize " + str(spacing_lat) +"\n"
    header += "NODATA_value 0"
    
    for i in range(ngas):

        for j in range(n_levels-1):

            conc = a[:, column + j + 1]

            conc = conc.reshape((-1,1))

            if np.sum(conc) == 0:

                print ( 'No volcanic gases of VG_'+str(i+1).zfill(2)+' H: '+str(H_LEVELS[j,0])+' m - '+str(H_LEVELS[j+1,0])+' m' )
                pass

            else:       

                # compute the range of values to plot
                min_conc = np.amin(conc)
                max_conc = np.amax(conc)


                half_conc = 0.5 * ( min_conc + max_conc )

                half_int = np.floor(np.log10(half_conc))
                base_scale = np.logspace(half_int-10,half_int-6, num=2)
                log_scale = np.logspace(half_int-4,half_int, num=5)
                lin_scale = np.linspace(half_conc, max_conc, num=8)

                # Assembre discrete levels to plot
                levels = np.append(log_scale,lin_scale)
                levels = np.append(base_scale,levels)
            
                mass_in_the_air = np.sum(conc[:,0] * dist[:,0] * dist[:,1] * (int(H_LEVELS[j+1,0])-int(H_LEVELS[j,0])))
                total_mass = total_mass + mass_in_the_air

                loading_i = loading3D[:,:,column + j + 1]
                loading3D_sum_kg += loading_i * pixel_area * (int(H_LEVELS[j+1,0])-int(H_LEVELS[j,0]))

                print ( 'VG ',str(i+1).zfill(2),': ','%.1e'%mass_in_the_air,' kg, H: '+str(H_LEVELS[j,0])+' m - '+str(H_LEVELS[j+1,0])+' m' )

            

                # Save the deposit for the single classes on a ESRI rater ascii file
                output_file = runname+'_VG'+str(i+1).zfill(2)+'_H_'+str(H_LEVELS[j,0])+'_'+str(H_LEVELS[j+1,0])+'_'+day+'_'+time+'.asc'

                np.savetxt(output_file, np.flipud(loading_i), \
                       header=header, fmt='%.3E',comments='')

                # Create a new figure
                f = plt.figure(i)
                plt.rcParams["font.size"] = 8.0
                cmap = plt.get_cmap('YlGnBu')
                norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            
                sm = Map(g.grid, factor=1, countries=False)
                sm.set_rgb(ggl_img)  # add the background rgb image
                sm.visualize()

                # Zero-loading pixels should not be plotted
                loading_i[loading_i==0.0] = np.nan
                Zm = np.ma.masked_where(np.isnan(loading_i),loading_i)

                x, y = sm.grid.transform(Lon_stag, Lat_stag)
                plt.pcolormesh(x,y,Zm, cmap=cmap, norm=norm,alpha=0.70,zorder=1)
                plt.pcolormesh(Lon_stag, Lat_stag,loading_i, cmap=cmap, norm=norm,alpha=1.0)

                xvent, yvent = sm.grid.transform(vent_lon, vent_lat)
                plt.plot(xvent, yvent,"^m",markersize=1,zorder=2)

                plt.xlim(left=np.amin(x))
                plt.xlim(right=np.amax(x))
                plt.ylim(bottom=np.amax(y))
                plt.ylim(top=np.amin(y))
                plt.grid()
                plt.title('VG'+str(i+1).zfill(2)+'\nH: '+str(H_LEVELS[j,0])+' - '+str(H_LEVELS[j+1,0])+' m\n'+'%.1e'%mass_in_the_air+' kg')
                clb = plt.colorbar(format=ticker.FuncFormatter(fmt))
                clb.set_label('conc. (kg/m$^3$)', labelpad=-40, y=1.05, rotation=0)

                f.savefig(runname+'_VG'+str(i+1).zfill(2)+'_H_'+str(H_LEVELS[j,0])+'_'+str(H_LEVELS[j+1,0])+'_'+day+'_'+time+'_AIR.pdf', bbox_inches='tight')
                plt.close()
            

        column = column + n_levels



f = plt.figure(i+1)
plt.rcParams["font.size"] = 8.0

loading3D_sum = loading3D_sum_kg / (pixel_area * (int(H_LEVELS[j+1,0])-int(H_LEVELS[j,0])))

min_conc = np.amin(loading3D_sum)
max_conc = np.amax(loading3D_sum)

half_conc = 0.5 * ( min_conc + max_conc )

half_int = np.floor(np.log10(half_conc))
  
base_scale = np.logspace(half_int-10,half_int-6, num=2)
log_scale = np.logspace(half_int-4,half_int, num=5)
lin_scale = np.linspace(half_conc, max_conc, num=8)

levels = np.append(log_scale,lin_scale)
levels = np.append(base_scale,levels)

cmap = plt.get_cmap('YlGnBu')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                     
sm = Map(g.grid, factor=1, countries=False)
sm.set_rgb(ggl_img)  # add the background rgb image
sm.visualize()

output_file = runname+'_'+'VG_sum'+'_'+day+'_'+time+'.asc'

np.savetxt(output_file, np.flipud(loading3D_sum), \
                       header=header, fmt='%.3E',comments='')


loading3D_sum[loading3D_sum<10**(-10)] = np.nan

x, y = sm.grid.transform(Lon_stag, Lat_stag)

Zm = np.ma.masked_where(np.isnan(loading3D_sum),loading3D_sum)
    
plt.pcolormesh(x,y,Zm, cmap=cmap, norm=norm,alpha=0.70,zorder=1)
plt.pcolormesh(Lon_stag, Lat_stag,loading3D_sum, cmap=cmap, norm=norm,alpha=1.0)

xvent, yvent = sm.grid.transform(vent_lon, vent_lat)
plt.plot(xvent, yvent,"^m",markersize=1,zorder=2)

plt.grid()
plt.xlim(left=np.amin(x))
plt.xlim(right=np.amax(x))
plt.ylim(bottom=np.amax(y))
plt.ylim(top=np.amin(y))
plt.title('Total atmospheric loading')
clb = plt.colorbar(format=ticker.FuncFormatter(fmt))
clb.set_label('conc. (kg/m$^3$)', labelpad=-40, y=1.05, rotation=0)
f.savefig(runname+'_'+'VG_sum'+'_'+day+'_'+time+'_AIR.pdf', bbox_inches='tight')
plt.close()
# plt.show()        
