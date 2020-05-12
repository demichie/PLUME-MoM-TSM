import matplotlib
matplotlib.use("TkAgg")

import sys
import os
import re
import pkg_resources
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
from mpl_toolkits.mplot3d import Axes3D

# increase the default widht of figures
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] *= 1.5 
plt.rcParams["figure.figsize"] = fig_size

if len(sys.argv)==2:

    filename = sys.argv[1] 
    file_name, file_extension = os.path.splitext(filename)
    print(file_extension)
    if ( file_extension != '.col' ):
        sys.exit()

else:

    required = {'easygui','tkinter'}
    installed = {pkg.key for pkg in pkg_resources.working_set}

    if ( 'easygui' in installed ):
        import easygui
        filename = easygui.fileopenbox( filetypes=['*.col'])
    elif ( 'tkinter' in installed ):

        from  tkinter import *
        root = Tk()
        root.filename =  filedialog.askopenfilename(title = "choose your file",filetypes=[("col files", "*.col")])
        filename = root.filename
        root.destroy()


bakfile = filename.replace('col','bak')


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
           print("aggregation_flag",aggregation_flag)

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
        labels.append(r'$\phi=$'+"{}".format(phi_max-j*delta_phi))
        

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


change_sign = np.argwhere(rho_rel[:-1]*rho_rel[1:]<0)

if ( len(change_sign ) > 0 ):
    last_change = change_sign[-1][0]
    change = True
else:
    last_change = -1
    change = False


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


plt.subplot(2, 2, 4)

lines = plt.plot(solid_tot_mass_fraction ,z, gas_mass_fraction,z,liquid_water_mass_fraction,z, ice_mass_fraction, z,'--')

plt.xlabel('Phases mass fraction')
plt.ylabel('Height (km)')
names = ['part','gas','lq','ice']
fig.tight_layout()
plt.legend(lines, [names[j] for j in range(len(names))])


plt.subplot(2, 2, 3)

n_bin_sample = np.linspace(0, n_bin-1, num=7, endpoint=True, dtype=int)

for i in range(n_part_sect):

    if ( np.remainder(i,n_bin) in n_bin_sample):
        plt.plot(solid_mass_fraction[:,i],z, color=colors[i],linestyle=linestyle_str[i],label=labels[i])
    else:
        plt.plot(solid_mass_fraction[:,i],z, color=colors[i],linestyle=linestyle_str[i],label='_nolegend_')

plt.legend(ncol=n_part,fontsize = 'x-small')


plt.xlabel('Particles mass fraction')
plt.ylabel('Height (km)')


fig.savefig(str(filename)+'_mass_fraction.pdf')   # save the figure to file
#plt.close()


# temperature

fig = plt.figure()

plt.plot(temp[:last_change]+273,z[:last_change])
if change:
    plt.plot(temp[last_change:]+273,z[last_change:])

plt.axvline(273, c = 'r',ls='--')
plt.axvline(233, c = 'r',ls='--')
#plt.plot(temp_atm,z,'.r')
plt.xlabel('Temp [K]')
plt.ylabel('Height (km)')
fig.savefig(str(filename)+'_temp.pdf')   # save the figure to file
#plt.close()

# PARTICLE LOSS FRACTION 

solid_mass_flux = np.zeros((results.shape[0],n_part_sect))

solid_mass_loss_cum = np.zeros((results.shape[0],n_part_sect))

for i in range(n_part_sect):

    solid_mass_flux[:,i] = rhoBsolid[:,i] * np.pi * r_mt[:,0]**2 * w[:,0] 

    solid_mass_loss_cum[:,i] =  1.0 - solid_mass_flux[:,i]/solid_mass_flux[0,i]

if ( not aggregation_flag ):

    fig = plt.figure()

    for i in range(n_part_sect):

        plt.plot(solid_mass_loss_cum[:,i],z, color=colors[i],linestyle=linestyle_str[i])

    plt.legend(labels,ncol=n_part,fontsize = 'x-small')
    plt.xlabel('Particles mass loss fraction')
    plt.ylabel('Height (km)')

    fig.savefig(str(filename)+'_particles_fraction.pdf')   # save the figure to file

    #plt.close()

solid_mass_flux_tot = np.sum(solid_mass_flux,axis=-1)
solid_mass_tot_loss_cum =  1.0 - solid_mass_flux_tot/solid_mass_flux_tot[0]



# VARIABLES

# two different colors are used below and above neutral buoyancy level
fig = plt.figure()

plt.subplot(2, 2, 1)

plt.plot(r[:last_change],z[:last_change])
if change:
    plt.plot(r[last_change:],z[last_change:])

plt.xlabel('Radius (km)')
plt.ylabel('Height (km)')

plt.subplot(2, 2, 2)

plt.plot(w[:last_change],z[:last_change])
if change:
    plt.plot(w[last_change:],z[last_change:])

plt.xlabel('Velocity (m/s)')
plt.ylabel('Height (km)')

plt.subplot(2, 2, 3)

plt.plot(rho_mix[:last_change],z[:last_change])
if change:
    plt.plot(rho_mix[last_change:],z[last_change:])

plt.xlabel('Mixture density (kg/m$^3$)')
plt.ylabel('Height (km)')

plt.subplot(2, 2, 4)

plt.plot(rho_rel[:last_change],z[:last_change])
if change:
    plt.plot(rho_rel[last_change:],z[last_change:])

#plt.plot(rho_atm,z,'.r')

plt.xlabel('Relative density (kg/m$^3$)')
plt.ylabel('Height (km)')

fig.tight_layout()
fig.savefig(str(filename)+'_profiles.pdf')   # save the figure to file
#plt.close()

# plot plume 3d

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho')
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

ax.plot(x[:last_change,0],y[:last_change,0],z[:last_change,0])
if change:
    ax.plot(x[last_change:,0],y[last_change:,0],z[last_change:,0])

# ax.scatter(x, y,z)

angle = np.linspace(0, 2*np.pi, num=50)
angle = angle.reshape((-1,1))

x_plume = np.cos(angle)
y_plume = np.sin(angle)

z_max = max(z)
z_min = min(z)

# n_sect = 50
n_sect = 20

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
    
    plume = np.zeros((3,x_plume.shape[0]))

    plume[0,:] = r[ind,0]*x_plume[:,0] 
    plume[1,:] = r[ind,0]*y_plume[:,0] 
    plume[2,:] = 0.0

    # ax.scatter(x[ind,0]+plume[0,:], y[ind,0]+plume[1,:],z[ind,0]+plume[2,:])
    if (ind<=last_change):
        ax.plot(x[ind,0]+plume[0,:], y[ind,0]+plume[1,:],z[ind,0]+plume[2,:],color=colors[0])
    else:
        ax.plot(x[ind,0]+plume[0,:], y[ind,0]+plume[1,:],z[ind,0]+plume[2,:],color=colors[1])


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
#ax4 = fig.add_subplot(133, frame_on=False)

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

# plt.subplot(1, 3, 1)

x_bin = np.flip(phi_min + delta_phi*np.arange(n_bin),axis=0)
x_bin = np.tile(x_bin,n_part)

solid_pmf_bot = np.zeros((n_part_sect))

for i in range(1,n_part):

    solid_pmf_bot[i*n_bin:(i+1)*n_bin] += solid_partial_mass_fraction[0,(i-1)*n_bin:i*n_bin]


barcollection1 = ax1.bar(x_bin, solid_partial_mass_fraction[0,:], width=0.9*delta_phi, bottom=solid_pmf_bot)

solid_pmf = np.sum(solid_partial_mass_fraction.reshape((n_levels,n_part,n_bin)),axis=1)
max_solid_pmf = np.nanmax(solid_pmf)
ax1.set_ylim(0.0, 1.1*max_solid_pmf)
ax1.set_xlim(phi_min-1,phi_max+1)
ax1.set_xlabel('phi')
ax1.title.set_text('Plume GSD')



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

# read cumulative (over z) values of sedimentation rates
sed_cum_results = np.loadtxt("%s.sed" % filename, skiprows = 1)

sed_cum_results = sed_cum_results.reshape((z_levels,-1))

# the values read from file are cumulative. 
# We want the rate of release between z[i-1] and z[i]
sed_results = np.zeros_like(sed_cum_results)
sed_results[1:,:] = np.diff(sed_cum_results,axis=0)

# sed_solid contains only sed rates (without first columns for z,x,y,r) 
sed_solid = np.zeros((sed_results.shape[0],n_part_sect))

# copy from sed_results to sed_solid
for i in range(n_part_sect):

    sed_solid[:,i] = sed_results[:,sed1["smf"+str(i+1)].column_org]

# normalize to have sum of fractions=1
for i in range(1,z_levels):

    sed_solid[i,:] /= np.sum(sed_solid[i,:])

sed_solid_bot = np.zeros((n_part_sect))

for i in range(1,n_part):

    sed_solid_bot[i*n_bin:(i+1)*n_bin] += sed_solid[0,(i-1)*n_bin:i*n_bin]


barcollection2 = ax2.bar(x_bin, sed_solid[0,:], width=0.9*delta_phi, bottom=sed_solid_bot)

# sum over different part belonging to the same bin (because of stacked bar)
sed_solid_pmf = np.sum(sed_solid.reshape((n_levels,n_part,n_bin)),axis=1)
# find the maximum over bins and levels of stacked bars
max_sed_solid_pmf = np.nanmax(sed_solid_pmf)

ax2.set_ylim(0.0, 1.1*max_sed_solid_pmf)
ax2.set_xlim(phi_min-1,phi_max+1)
ax2.title.set_text('Sedimentation GSD')
ax2.set_xlabel('phi')


# ax = plt.subplot(1, 3, 3)

ax3.plot(solid_mass_tot_loss_cum[:last_change],z[:last_change])
if change:
    ax3.plot(solid_mass_tot_loss_cum[last_change-1:],z[last_change-1:])


ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()

mark_pos, = ax3.plot(solid_mass_tot_loss_cum[0],z[0],'o')
ax3.set_ylabel('Height [km]')
# plt.xlabel('Fraction of solid flux lost')
ax3.title.set_text('Fraction of solid flux lost')

title = ax3.text(0.70,0.05,"t="+"{:6.1f}".format(float(time_steps[0]))+'s', bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax3.transAxes, ha="center")



def animate(i):

    y1 = solid_partial_mass_fraction[idx_steps[i],:]

    solid_pmf_bot[0:n_part_sect] = 0.0

    for i_part in range(1,n_part):

        solid_pmf_bot[i_part*n_bin:(i_part+1)*n_bin] += solid_partial_mass_fraction[idx_steps[i],(i_part-1)*n_bin:i_part*n_bin]

   
    for j, b in enumerate(barcollection1):
        b.set_height(y1[j])
        b.set_y(solid_pmf_bot[j])
        b.set_color(bar_colors[j])


    y2 = sed_solid[idx_steps[i],:]

    sed_solid_bot[0:n_part_sect] = 0.0

    for i_part in range(1,n_part):

        sed_solid_bot[i_part*n_bin:(i_part+1)*n_bin] += sed_solid[idx_steps[i],(i_part-1)*n_bin:i_part*n_bin]
   
    for j, b in enumerate(barcollection2):
        b.set_height(y2[j])
        b.set_y(sed_solid_bot[j])
        b.set_color(bar_colors[j])

    mark_pos.set_xdata(solid_mass_tot_loss_cum[idx_steps[i]])  
    mark_pos.set_ydata(z[idx_steps[i]])  

    title.set_text("t="+"{:6.1f}".format(float(time_steps[i]))+'s')

anim=animation.FuncAnimation(fig,animate,repeat=False,blit=False,frames=n_frames,
                             interval=100)


anim.save(str(filename)+'_anim.mp4',writer=animation.FFMpegWriter(fps=10))



plt.show()


