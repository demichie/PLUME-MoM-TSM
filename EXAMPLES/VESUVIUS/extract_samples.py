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
        print '{:03d}'.format(i),'{:15.10f}'.format(eval(checkP)[0]),'{:15.10f}'.format(eval(checkP)[1])

        f1.write('{:03d}'.format(i)+'{:15.10f}'.format(eval(checkP)[0])+'{:15.10f}'.format(eval(checkP)[1]))
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

block_length = nsampl * npart

print 'nsampl',nsampl,'block_length',block_length

sampling_cmd = hysplit_dir+'/exec/con2stn -i'+'cdumpcum_part_'+runname+' -scon2stn.inp -d0 -p0 -xi -z1 -r0 -ocon2stn.txt'

os.system(sampling_cmd)

con2stn = np.loadtxt('con2stn.txt',skiprows=1)

# print 'size',con2stn.shape

JDAY = con2stn[:,0]

Value = con2stn[:,-1]

# print 'JDAY',JDAY
# print 'Value',Value

con2std_len = len(JDAY)

nblocks = con2std_len/block_length

loading = np.zeros((nblocks,nsampl,npart))

for i in range(nblocks):
    for j in range(nsampl):
        loading[i,j,:] = Value[i*block_length+(j+1)+nsampl*np.arange(npart)-1]


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
ftot = open(deposit_file_total,'wb')

chi_squared = 0.0
chi_squared_j = 0.0

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

        deposit_file_j = deposit_file+'{:03d}'.format(j)+'.txt'
        f = open(deposit_file_j,'wb')

        out_data = np.vstack((np.array(diam_phi), loading[nblocks-1,j,:]/np.sum(loading[nblocks-1,j,:])*100))
      	

        f.write(str(out_data.T))

        f.close()


ftot.close()

deposit_file_chi = deposit_file+'_chi_squared.txt'

fchi = open(deposit_file_chi,'wb')
fchi.write('{:20.10f}'.format(chi_squared[0]))
fchi.close()


ax.legend(legend_strings)

ax.set_ylabel('Loading [mass wt%]')
ax.set_title('Loading at sampling locs')

ax.set_xlabel('Diameter [phi]')
ax.set_xticks(np.array(diam_phi))

# ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))

fig.savefig('gsd.pdf')

