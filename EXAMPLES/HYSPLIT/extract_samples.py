import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from extract_wind import write_atm
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pylab as plot

from input_file import *

f1=open('./con2stn.inp', 'w+')


sampled_deposit = np.zeros((100,1))

for i in range(0,100):
    checkP = 'P'+'{:02d}'.format(i)
    if (checkP in vars()):
        print ( '{:03d}'.format(i),'{:15.10f}'.format(eval(checkP)[0]),'{:15.10f}'.format(eval(checkP)[1]) )

        f1.write('{:03d}'.format(i)+' {:15.10f} '.format(eval(checkP)[0])+' {:15.10f} '.format(eval(checkP)[1]))
        f1.write("\n")

        if len(eval(checkP))==3:

            sampled_deposit[i-1] =eval(checkP)[2]

f1.close()

 
#'a' in vars() or 'a' in globals()

INDEX, LAT, LON = np.loadtxt('con2stn.inp',skiprows=0, unpack=True)

if isinstance(INDEX, (np.ndarray) ):

    nsampl = len(INDEX)

else:

    nsampl = 1
    INDEX = np.array([INDEX])
    LAT = np.array([LAT])
    LON = np.array([LON])

block_length = nsampl * n_sections

print ( 'nsampl',nsampl,'block_length',block_length )

sampling_cmd = hysplit_dir+'/exec/con2stn -i'+'cdumpcum_part_'+runname+' -scon2stn.inp -d0 -p0 -xi -z1 -r0 -ocon2stn.txt'

os.system(sampling_cmd)

con2stn = np.loadtxt('con2stn.txt',skiprows=1)

# print 'size',con2stn.shape

JDAY = con2stn[:,0]

Value = con2stn[:,-1]

# print 'JDAY',JDAY
# print 'Value',Value

con2std_len = len(JDAY)

nblocks = int(con2std_len/block_length)

loading = np.zeros((nblocks,nsampl,n_sections))

for i in range(nblocks):
    for j in range(nsampl):
        loading[i,j,:] = Value[i*block_length+(j+1)+nsampl*np.arange(n_sections)-1]


fig = plt.figure()
ax = plt.subplot(111)

params = {'legend.fontsize': 6,
	  'legend.handlelength': 2}
plot.rcParams.update(params)

width = 1.5/nsampl

legend_strings = []
color=iter(cm.rainbow(np.linspace(0,1,nsampl)))

deposit_file = 'sample_deposit'


deposit_file_total = deposit_file+'_total.txt'
ftot = open(deposit_file_total,'w')

chi_squared = 0.0
chi_squared_j = 0.0

diam_phi = np.linspace(phi_min+(delta_phi),phi_min+(delta_phi)+((n_sections-1)*delta_phi),n_sections, endpoint=True)

for j in range(nsampl):

    checkP = 'P'+'{:02d}'.format(j+1)

    sim_deposit = np.sum(loading[nblocks-1,j,:]) 

    if ( sampled_deposit[j]>0 ):
    
        chi_squared_j = ( np.sum(loading[nblocks-1,j,:]) - sampled_deposit[j] )**2 / sampled_deposit[j]
        ftot.write('{:15.10f}'.format(eval(checkP)[0])+'{:15.10f}'.format(eval(checkP)[1])  \
                   +'{:18.10f}'.format(sim_deposit)+'{:20.10f}'.format(chi_squared_j[0])+'\n')

        chi_squared += chi_squared_j

    else:


        ftot.write('{:15.10f}'.format(eval(checkP)[0])+'{:15.10f}'.format(eval(checkP)[1])  \
                   +'{:15.10f}'.format(sim_deposit)+'\n')


    if ( np.sum(loading[nblocks-1,j,:]) > 0.0 ):


        c=next(color)
        ax.bar(np.array(diam_phi)+(j-0.5*nsampl)*width, loading[nblocks-1,j,:]/np.sum(loading[nblocks-1,j,:])*100,width,color=c)
        stringj = "Loc %s (%s,%s), Loading=%.2e [kg/m2]" % (str(j+1), str(LAT[j]), str(LON[j]),np.sum(loading[nblocks-1,j,:]) )
        legend_strings.append(stringj)

        fig2 = plt.figure()
        ax2 = plt.subplot(111)
        ax2.bar(np.array(diam_phi), loading[nblocks-1,j,:]/np.sum(loading[nblocks-1,j,:])*100,color=c)
        ax2.legend(["Loc %s (%s,%s), Loading=%.2e [kg/m2]" % (str(j+1), str(LAT[j]), str(LON[j]),np.sum(loading[nblocks-1,j,:]))],loc='upper left')
        ax2.set_ylabel('Loading [mass wt%]')
        ax2.set_title('Loading at sampling loc '+str(j+1).zfill(3))
        ax2.set_xlabel('Diameter [phi]')
        ax2.set_xticks(np.array(diam_phi))
        fig2.savefig('gsd_Loc'+str(j+1).zfill(3)+'.pdf')
        plt.close()

        deposit_file_j = deposit_file+'{:03d}'.format(j)+'.txt'
        f = open(deposit_file_j,'w')

        out_data = np.vstack((np.array(diam_phi), loading[nblocks-1,j,:]/np.sum(loading[nblocks-1,j,:])*100))
      	

        f.write(str(out_data.T))

        f.close()


ftot.close()

deposit_file_chi = deposit_file+'_chi_squared.txt'

fchi = open(deposit_file_chi,'w')
fchi.write('{:20.10f}'.format(chi_squared[0]))
fchi.close()

ax.legend(legend_strings)
ax.set_ylabel('Loading [mass wt%]')
ax.set_title('Loading at sampling locs')
ax.set_xlabel('Diameter [phi]')
ax.set_xticks(np.array(diam_phi))

fig.savefig('gsd.pdf')

"""
# plot sampling points on a map using a *.tif image as background
import georaster
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

llcrnrlon0=np.min([np.min(LON),vent_lon]) -0.005
urcrnrlon0=np.max([np.max(LON),vent_lon]) +0.005

llcrnrlat0=np.min([np.min(LAT),vent_lat]) -0.005
urcrnrlat0=np.max([np.max(LAT),vent_lat]) +0.005

fig = plt.figure(figsize=(8,8))

# full path to the geotiff file
fpath = "ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif"  # Tif image

# read extent of image without loading
# good for values in degrees lat/long
# geotiff may use other coordinates and projection
my_image = georaster.SingleBandRaster(fpath, load_data=False)

# grab limits of image's extent
minx = llcrnrlon0
maxx = urcrnrlon0
miny = llcrnrlat0
maxy = urcrnrlat0

# set Basemap with slightly larger extents
# set resolution at intermediate level "i"
m = Basemap( projection='cyl', \
            llcrnrlon=minx, \
            llcrnrlat=miny, \
            urcrnrlon=maxx, \
            urcrnrlat=maxy, \
            resolution='i')

m.drawcoastlines(color="gray")
m.fillcontinents(color='beige')

# load the geotiff image, assign it a variable
image = georaster.SingleBandRaster( fpath, \
                        load_data=(minx, maxx, miny, maxy), \
                        latlon=True)

image2 = np.flip(image.r, axis=0)

Xa = np.linspace(minx,maxx,image.r.shape[1])
Ya = np.linspace(miny,maxy,image.r.shape[0])
X,Y = np.meshgrid(Xa,Ya)

# plot the image on matplotlib active axes
# set zorder to put the image on top of coastlines and continent areas
# set alpha to let the hidden graphics show through
m.imshow(image2, extent=(minx, maxx, miny, maxy), zorder=5, alpha=0.6,cmap='winter')
m.contour(X,Y,image2, colors="white", levels = list(range(0, 1000, 10)),zorder=10, alpha=0.8, linewidths=0.5)

x_vent, y_vent = m(vent_lon, vent_lat)
m.scatter(x_vent, y_vent,zorder=15,s=7,c="k", marker = "^")

x, y = m(LON, LAT)
for j in range(nsampl):

    m.scatter(x[j], y[j],zorder=15,s=2,c="k")
    plt.annotate("P"+str(j+1),(x[j],y[j]),zorder=15)

fig.savefig('sample_locs_map.pdf')
"""
