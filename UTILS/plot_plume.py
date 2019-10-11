import matplotlib
matplotlib.use("TkAgg")

#import Tkinter, tkFileDialog
import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
from mpl_toolkits.mplot3d import Axes3D
import easygui

# increase the default widht of figures
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] *= 1.5 
plt.rcParams["figure.figsize"] = fig_size


#option 1
#filename = easygui.fileopenbox( filetypes=['*.col'])

#option 2 (in case option 1 doesn't work)
from tkFileDialog import askopenfilename
filename = askopenfilename(filetypes=[("col files", "*.col")])

bakfile = filename.replace('col','bak')
print(filename)
print(bakfile)

with open(bakfile) as fp:  
   for cnt, line in enumerate(fp):
       if "N_SECTIONS" in line:
           n_sections_str= line.replace('N_SECTIONS=','')
           n_sections_str= n_sections_str.replace(',','')
           n_sections = np.int(n_sections_str)
           print("n_sections",n_sections)
       if "PHI_MIN" in line:
           phi_min_str = line.replace('PHI_MIN=','')
           phi_min_str = phi_min_str.replace(',','')
           phi_min = np.float(phi_min_str)
           print("phi_min",phi_min)
       if "DELTA_PHI" in line:
           delta_phi_str = line.replace('DELTA_PHI=','')
           delta_phi_str = delta_phi_str.replace(',','')
           delta_phi = np.float(delta_phi_str)
           print("delta_phi",delta_phi)

       if "AGGREGATION_FLAG" in line:
           aggregation_flag_str = line.replace('AGGREGATION_FLAG=','')
           aggregation_flag = ('T' in aggregation_flag_str)
           print("aggregation_flag",aggregation_flag,aggregation_flag_str)

phi_max = phi_min + (n_sections-1)*delta_phi


filename = filename.split('/')[-1]
filename = re.sub('\.col$', '', filename)


class smf_array:
  pass


d1={}

with open("%s.col" % filename, "r") as file1:
    line=file1.readlines()
    header_line=line[0]
    header_split = header_line.split()
    n_part = 0
    n_part_sect = 0
    n_bin = 0
    n_gas = 0


    for i in range(len(header_split)):

        d1["smf{0}".format(i)] = smf_array() 

        if header_split[i][0:4] == "rhoB":

            last_part = header_split[i][4:]

            n_part_sect = n_part_sect +1
                         
            d1["smf{0}".format(n_part_sect)].org = n_part_sect 

            d1["smf{0}".format(n_part_sect)].column_org = i
       
        elif header_split[i] == "volgas.massf":
             n_gas = n_gas + 1

    n_part = int(last_part[0:2])
    n_bin = int(last_part[3:])

       
print( 'number of particles phases', n_part)
print( 'number of bins ', n_bin)
print( 'number of volcanic gases ',n_gas)

labels =[]
for i in range(n_part):
    for j in range(n_bin):
        labels.append(r'$\phi=$'+"{}".format(phi_min+j*delta_phi))
        

results = np.loadtxt("%s.col" % filename, skiprows = 1)


z_levels = results.shape[0]


results=results.reshape((z_levels,-1))


z = results[:,0]/float(1000)
r_mt = results[:,1]
r = results[:,1]/float(1000)
x = results[:,2]/float(1000)
y = results[:,3]/float(1000)
rho_mix = results[:,4]
temp = results[:,5]
w = results[:,6]
mag_u = results[:,7]
dry_air_mass_fraction = results[:,8]
wvapour_mass_fraction = results[:,9]
liquid_water_mass_fraction = results[:,10]
ice_mass_fraction = results[:,11]

z=z.reshape((-1,1))
r_mt = r_mt.reshape((-1,1))
r = r.reshape((-1,1))
x = x.reshape((-1,1))
y = y.reshape((-1,1))
rho_mix = rho_mix.reshape((-1,1))
temp = temp.reshape((-1,1))
w = w.reshape((-1,1))
mag_u = mag_u.reshape((-1,1))
dry_air_mass_fraction = dry_air_mass_fraction.reshape((-1,1))
wvapour_mass_fraction = wvapour_mass_fraction.reshape((-1,1))
liquid_water_mass_fraction = liquid_water_mass_fraction.reshape((-1,1))
ice_mass_fraction = ice_mass_fraction.reshape((-1,1))

n_levels = results.shape[0]

solid_partial_mass_fraction = np.zeros((results.shape[0],n_part_sect))
rhoBsolid = np.zeros((results.shape[0],n_part_sect))

for i in range(n_part_sect):

    rhoBsolid[:,i] = results[:,d1["smf"+str(i+1)].column_org]

rhoBsolidTot = np.sum(rhoBsolid, axis=1)
rhoBsolidTot = rhoBsolidTot.reshape((-1,1))

for i in range(n_levels):

    solid_partial_mass_fraction[i,:] = rhoBsolid[i,:] / rhoBsolidTot[i]


if n_gas == 0:

    volcgas_mass_fraction = np.zeros((results.shape[0],1))
    volcgas_mix_mass_fraction = np.zeros((results.shape[0],1))

else:

    volcgas_mass_fraction = np.zeros((results.shape[0],n_gas))

    for i in range(n_gas):

        volcgas_mass_fraction[:,i] = results[:,12+n_part+i]

    volcgas_mix_mass_fraction = results[:,12+n_part+n_gas]

#print volcgas_mass_fraction

volcgas_mass_fraction_tot = np.sum(volcgas_mass_fraction, axis = 1)
volcgas_mass_fraction_tot=volcgas_mass_fraction_tot.reshape((-1,1))

gas_mass_fraction = np.zeros((results.shape[0],1))

for i in range(gas_mass_fraction.shape[0]):
    gas_mass_fraction[i,0] = dry_air_mass_fraction[i,0] + wvapour_mass_fraction[i,0] + volcgas_mass_fraction_tot[i,0] 


solid_mass_fraction = np.zeros((results.shape[0],n_part_sect))

for i in range(n_part_sect):
    solid_mass_fraction[:,i] = solid_partial_mass_fraction[:,i] * ( 1 - gas_mass_fraction[:,0] - ice_mass_fraction[:,0] - liquid_water_mass_fraction[:,0])
       



solid_tot_mass_fraction = np.zeros((results.shape[0],1))
solid_tot_mass_fraction[:,0] = np.sum(solid_mass_fraction,axis=1)



rho_atm = results[:,12+n_part_sect+n_gas+1]
rho_atm = rho_atm.reshape((-1,1))

mfr = results[:,12+n_part_sect+n_gas+2]
mfr = mfr.reshape((-1,1))

temp_atm = results[:,12+n_part_sect+n_gas+3]
temp_atm = temp_atm.reshape((-1,1))

p_atm = results[:,12+n_part_sect+n_gas+4]
p_atm = p_atm.reshape((-1,1))


n_z = z.shape[0]

rho_rel = rho_mix - rho_atm
rho_rel = rho_rel.reshape((-1,1))


# PLOT FIGURES

cm_subsection = np.linspace(0.0,1.0,n_bin) 

colors = [ cm.jet(xj) for xj in cm_subsection ]

colors = np.tile(colors,(n_part,1))

linestyle_str = [ 'solid','dotted','dashed','dashdot'] 

linestyle_str = np.repeat(linestyle_str,n_bin)


# MASS FRACTION 

fig = plt.figure()

plt.subplot(2, 2, 1)

lines = plt.plot(dry_air_mass_fraction,z, volcgas_mix_mass_fraction,z,wvapour_mass_fraction,z,gas_mass_fraction,z)

names = ['dry air','volcgas','wv','totalgas']

plt.legend(lines, [names[j] for j in range(len(names))])
plt.xlabel('Gas mass fraction')
plt.ylabel('Height (km)')

plt.subplot(2, 2, 2)

water = wvapour_mass_fraction + liquid_water_mass_fraction + ice_mass_fraction

lines = plt.plot(wvapour_mass_fraction/water,z,'-',liquid_water_mass_fraction/water,z,'-', ice_mass_fraction/water, z,'-')

plt.xlabel('Water mass fraction')
plt.ylabel('Height (km)')
names = ['wv','lq','ice']

plt.legend(lines, [names[j] for j in range(len(names))])

plt.subplot(2, 2, 3)

for i in range(n_part_sect):

    plt.plot(solid_mass_fraction[:,i],z, color=colors[i],linestyle=linestyle_str[i])
    
plt.legend(labels,ncol=n_part,fontsize = 'x-small')


plt.xlabel('Particles mass fraction')
plt.ylabel('Height (km)')

plt.subplot(2, 2, 4)

lines = plt.plot(solid_tot_mass_fraction ,z, gas_mass_fraction,z,liquid_water_mass_fraction,z, ice_mass_fraction, z,'--')

plt.xlabel('Phases mass fraction')
plt.ylabel('Height (km)')
names = ['part','gas','lq','ice']
plt.legend(lines, [names[j] for j in range(len(names))])
fig.tight_layout()
fig.savefig(str(filename)+'_mass_fraction.pdf')   # save the figure to file
#plt.close()


# temperature

fig = plt.figure()


plt.plot(temp+273,z)
plt.axvline(273, c = 'r')
plt.axvline(233, c = 'r')
#plt.plot(temp_atm,z,'.r')
plt.xlabel('Temp [K]')
plt.ylabel('Height (km)')
fig.savefig(str(filename)+'_temp.pdf')   # save the figure to file
#plt.close()

# PARTICLE LOSS FRACTION 
fig = plt.figure()

solid_mass_flux = np.zeros((results.shape[0],n_part_sect))

solid_mass_loss_cum = np.zeros((results.shape[0],n_part_sect))

for i in range(n_part_sect):

    solid_mass_flux[:,i] = rhoBsolid[:,i] * np.pi * r_mt[:,0]**2 * mag_u[:,0] 

    solid_mass_loss_cum[:,i] =  1.0 - solid_mass_flux[:,i]/solid_mass_flux[0,i]

    plt.plot(solid_mass_loss_cum[:,i],z, color=colors[i],linestyle=linestyle_str[i])

plt.legend(labels,ncol=n_part,fontsize = 'x-small')

plt.xlabel('Particles mass loss fraction')
plt.ylabel('Height (km)')



fig.savefig(str(filename)+'_particles_fraction.pdf')   # save the figure to file
#plt.close()

solid_mass_flux_tot = np.sum(solid_mass_flux,axis=-1)
solid_mass_tot_loss_cum =  1.0 - solid_mass_flux_tot/solid_mass_flux_tot[0]



# VARIABLES

fig = plt.figure()

plt.subplot(2, 2, 1)

plt.plot(r,z)

plt.xlabel('Radius (km)')
plt.ylabel('Height (km)')

plt.subplot(2, 2, 2)

plt.plot(w,z)

plt.xlabel('Velocity (m/s)')
plt.ylabel('Height (km)')

plt.subplot(2, 2, 3)

plt.plot(rho_mix,z)

plt.xlabel('Mixture density (kg/m$^3$)')
plt.ylabel('Height (km)')

plt.subplot(2, 2, 4)

plt.plot(rho_rel,z)
#plt.plot(rho_atm,z,'.r')

plt.xlabel('Relative density (kg/m$^3$)')
plt.ylabel('Height (km)')

fig.tight_layout()
fig.savefig(str(filename)+'_profiles.pdf')   # save the figure to file
#plt.close()

# plot plume 3d

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y,z)

angle = np.linspace(0, 2*np.pi, num=50)
angle = angle.reshape((-1,1))

x_plume = np.cos(angle)
y_plume = np.sin(angle)

z_max = max(z)
z_min = min(z)

n_sect = 50

zeta_grid = np.linspace(z_min,z_max*0.99,num = n_sect)

l_seg = []

for i in range(1,x.shape[0],1):
    l_seg.append(((x[i,0]-x[i-1,0])**2 + (y[i,0]-y[i-1,0])**2 + (z[i,0]-z[i-1,0])**2)**0.5)

l_seg = np.asarray(l_seg)
l_seg = l_seg.reshape((-1,1))
    
s_axis = np.zeros((x.shape[0],1))

s_axis[0,0] = 0
s_axis[1:,0] =  np.cumsum (l_seg[:,0])

s_grid = np.linspace(0,max(s_axis)*0.99,num = n_sect)

for i in range(n_sect):
  
    ind0 =  np.where(s_axis[:,0]-s_grid[i]>0)
    ind0 = np.asarray(ind0)
    ind0 = ind0.reshape((-1,1))
    ind = min(ind0[:,0])
    
    vect = np.zeros((1,3))  
   
    vect[0,0] = x[ind,0] - x[ind-1,0]
    vect[0,1] = y[ind,0] - y[ind-1,0]
    vect[0,2] = z[ind,0] - z[ind-1,0]

    vect = vect / float(np.linalg.norm(vect, ord=2)) 
    
    vect0 = np.zeros((1,3))  

    vect0[0,0] = 0
    vect0[0,1] = 0
    vect0[0,2] = 1

    v = np.cross(vect0,vect)

    s = np.linalg.norm(v, ord=2)

    c = np.vdot(vect0,vect)

   
    mat_v = np.zeros((3,3))

    mat_v[1,0] = v[0,2]
    mat_v[0,1] = -v[0,2]
   
    mat_v[2,0] = -v[0,1]
    mat_v[0,2] = v[0,1]
   
    mat_v[1,2] = -v[0,0]
    mat_v[2,1] = v[0,0]

    R = np.eye(3) + mat_v + np.dot(mat_v,mat_v) * (1 - c) / s**2

    plume = np.zeros((3,x_plume.shape[0]))

    plume[0,:] = r[ind,0]*x_plume[:,0] 
    plume[1,:] = r[ind,0]*y_plume[:,0] 

    plume_rotated = np.dot(R, plume)

    ax.scatter(x[ind,0]+plume_rotated[0,:], y[ind,0]+plume_rotated[1,:],z[ind,0]+plume_rotated[2,:])

ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_zlabel('z (km)')
fig.tight_layout()   
fig.savefig(str(filename)+'_plume.pdf')   # save the figure to file

#--------- Animation of evolution of GSD

from matplotlib import animation

cm_subsection = np.linspace(0.0,1.0,n_part) 

bar_colors = [ cm.jet(xj) for xj in cm_subsection ]

bar_colors = np.tile(bar_colors,n_bin)
bar_colors = bar_colors.reshape((n_part_sect,4))

fig_size[0] *= 2.5 
plt.rcParams["figure.figsize"] = fig_size

fig=plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax4 = fig.add_subplot(133, frame_on=False)

time = np.zeros((x.shape[0],1))
time[0] = 0.0
for i in range(1,n_levels):

    dz = 1000.0*(z[i,0] - z[i-1,0])
    dt = dz/(0.5*(w[i]+w[i-1]))
    time[i] = time[i-1]+dt
    # print dt,dz,0.5*(w[i]+w[i-1])

n_frames = 100

time_steps = np.linspace(0,time[-1],n_frames)

idx_steps = []
for value in time_steps:
    idx_steps.append((np.abs(time - value)).argmin())

plt.subplot(1, 3, 1)

x_bin = np.flip(phi_min + delta_phi*np.arange(n_bin),axis=0)
x_bin = np.tile(x_bin,n_part)

solid_pmf_bot = np.zeros((n_part_sect))

for i in range(1,n_part):

    solid_pmf_bot[i*n_bin:(i+1)*n_bin] += solid_partial_mass_fraction[0,(i-1)*n_bin:i*n_bin]


barcollection1 = plt.bar(x_bin, solid_partial_mass_fraction[0,:], width=0.9*delta_phi, bottom=solid_pmf_bot)

solid_pmf = np.sum(solid_partial_mass_fraction.reshape((n_levels,n_part,n_bin)),axis=1)
max_solid_pmf = np.max(solid_pmf)
plt.ylim(0.0, 1.1*max_solid_pmf)
plt.xlim(phi_min-1,phi_max+1)
plt.xlabel('phi')
ax1.title.set_text('Plume GSD')



plt.subplot(1, 3, 2)

sed1={}

with open("%s.sed" % filename, "r") as file1:
    line=file1.readlines()
    header_line=line[0]
    header_split = header_line.split()
    n_part_sect = 0

    for i in range(len(header_split)):

        sed1["smf{0}".format(i)] = smf_array() 

        if header_split[i][0:4] == "rhoB":

            last_part = header_split[i][4:]

            n_part_sect = n_part_sect +1
             
            sed1["smf{0}".format(n_part_sect)].org = n_part_sect 

            sed1["smf{0}".format(n_part_sect)].column_org = i

sed_results = np.loadtxt("%s.sed" % filename, skiprows = 1)

sed_results = sed_results.reshape((z_levels,-1))

# the values read from file are cumulative. We want the diff
sed_results[1:,:] = np.diff(sed_results,axis=0)

sed_solid_partial_mass_fraction = np.zeros((sed_results.shape[0],n_part_sect))
sed_solid = np.zeros((sed_results.shape[0],n_part_sect))


for i in range(n_part_sect):

    sed_solid[:,i] = sed_results[:,sed1["smf"+str(i+1)].column_org]


sed_solidTot = np.sum(sed_solid, axis=1)
sed_solidTot = sed_solidTot.reshape((-1,1))

sed_solid_partial_mass_fraction[0,:] = sed_solid[0,:] 

for i in range(1,n_levels):

    sed_solid_partial_mass_fraction[i,:] = sed_solid[i,:] / sed_solidTot[i]

sed_solid_pmf_bot = np.zeros((n_part_sect))

for i in range(1,n_part):

    sed_solid_pmf_bot[i*n_bin:(i+1)*n_bin] += sed_solid_partial_mass_fraction[0,(i-1)*n_bin:i*n_bin]


barcollection2 = plt.bar(x_bin, sed_solid_partial_mass_fraction[0,:], width=0.9*delta_phi, bottom=sed_solid_pmf_bot)

sed_solid_pmf = np.sum(sed_solid_partial_mass_fraction.reshape((n_levels,n_part,n_bin)),axis=1)
max_sed_solid_pmf = np.max(sed_solid_pmf)
plt.ylim(0.0, 1.1*max_sed_solid_pmf)
plt.xlim(phi_min-1,phi_max+1)
ax2.title.set_text('Sedimentation GSD')
plt.xlabel('phi')


ax = plt.subplot(1, 3, 3)

plt.plot(solid_mass_tot_loss_cum,z)


ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

mark_pos, = plt.plot(solid_mass_tot_loss_cum[0],z[0],'o')
plt.ylabel('Height [km]')
# plt.xlabel('Fraction of solid flux lost')
ax3.title.set_text('Fraction of solid flux lost')

title = ax.text(0.70,0.05,"t="+"{:6.1f}".format(time_steps[0])+'s', bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")



""" 

mark_pos, = ax3.plot(np.sqrt(x[0]**2+y[0]**2),z[0],'o')

ax3.plot(np.sqrt(x**2+y**2),z)

ax4.set_ylabel('Height [km]')
ax4.set_xlabel('[km]')
#ax3.title.set_text('Position on plume axis')
ax4.set_xlim(0,np.amax(np.sqrt(x**2+y**2)))

ax4.xaxis.tick_top()
ax4.tick_params(axis='x', colors="r")
ax4.xaxis.set_label_position('top') 



ax3.plot(solid_mass_tot_loss_cum,z,color='k')
mark_pos, = ax3.plot(solid_mass_tot_loss_cum[0],z[0],'o')
ax3.set_xlim(-0.00001,np.amax(solid_mass_tot_loss_cum)+0.00001)
ax3.set_xlabel('fraction of solid flux lost',color='k') 

ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()
ax3.xaxis.tick_bottom()

ax3.set_ylabel('Height [km]')



title = ax3.text(1.7,0.15,"t="+"{:6.1f}".format(time_steps[0])+'s', bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")

"""

def animate(i):

    y1 = solid_partial_mass_fraction[idx_steps[i],:]

    solid_pmf_bot[0:n_part_sect] = 0.0

    for i_part in range(1,n_part):

        solid_pmf_bot[i_part*n_bin:(i_part+1)*n_bin] += solid_partial_mass_fraction[idx_steps[i],(i_part-1)*n_bin:i_part*n_bin]

   
    for j, b in enumerate(barcollection1):
        b.set_height(y1[j])
        b.set_y(solid_pmf_bot[j])
        b.set_color(bar_colors[j])


    y2 = sed_solid_partial_mass_fraction[idx_steps[i],:]

    sed_solid_pmf_bot[0:n_part_sect] = 0.0

    for i_part in range(1,n_part):

        sed_solid_pmf_bot[i_part*n_bin:(i_part+1)*n_bin] += sed_solid_partial_mass_fraction[idx_steps[i],(i_part-1)*n_bin:i_part*n_bin]
   
    for j, b in enumerate(barcollection2):
        b.set_height(y2[j])
        b.set_y(sed_solid_pmf_bot[j])
        b.set_color(bar_colors[j])

    mark_pos.set_xdata(solid_mass_tot_loss_cum[idx_steps[i]])  
    mark_pos.set_ydata(z[idx_steps[i]])  

    title.set_text("t="+"{:6.1f}".format(time_steps[i])+'s')

anim=animation.FuncAnimation(fig,animate,repeat=False,blit=False,frames=n_frames,
                             interval=100)

#anim=animation.FuncAnimation(fig,animate,repeat=False,blit=False,frames=3,
#                             interval=100)

anim.save('mymovie.mp4',writer=animation.FFMpegWriter(fps=10))



plt.show()


