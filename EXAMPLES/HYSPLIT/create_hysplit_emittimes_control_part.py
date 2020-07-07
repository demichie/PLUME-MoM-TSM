import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
import math as math
from collections import Counter

from part_density import calc_density
from input_file import *

def round_minutes(dt, direction, resolution):

    if ( dt.minute%resolution == 0 ):

        rounded_time = dt

    else: 

        new_minute = (dt.minute // resolution + (1 if direction == 'up' else 0)) * resolution

        rounded_time = dt + datetime.timedelta(minutes=new_minute - dt.minute)

    return rounded_time

time_format = "%y %m %d %H %M"

starttime_hhmm = datetime.datetime.strptime(starttime,time_format)
starttime_round = round_minutes(starttime_hhmm, 'down', 60) # arrotonda per difetto starttime

endruntime_hhmm = datetime.datetime.strptime(endruntime,time_format)
endruntime_round = round_minutes(endruntime_hhmm, 'up', 60) # arrotonda per eccesso endemittime

runtime=endruntime_round-starttime_round

d = datetime.datetime(2000,1,1) + runtime
runtime_hh = '{0:02}'.format( int(runtime.total_seconds()//3600) )

print ( 'Start time:',starttime )
print ( 'End run time:',endruntime )
print ( 'Total runtime',runtime_hh,'hrs' )

# to work with deltat_plumemom < 3600 

if deltat_plumemom >=3600:

    starttime_hhmm = datetime.datetime.strptime(starttime,time_format) 
    starttime_round = round_minutes(starttime_hhmm, 'down', 60) # arrotonda per difetto starttime

    endemittime_hhmm = datetime.datetime.strptime(endemittime,time_format)
    endemittime_round = round_minutes(endemittime_hhmm, 'up', 60) # arrotonda per eccesso endemittime

else:

    # for deltat_plumemom = 900 or 1800

    starttime_hhmm = datetime.datetime.strptime(starttime,time_format) 
    starttime_round = round_minutes(starttime_hhmm, 'down', 15)  # arrotonda per difettos endemittime

    starttime_round_c = round_minutes(starttime_hhmm, 'down', 60)  # arrotonda per difettos endemittime

    endemittime_hhmm = datetime.datetime.strptime(endemittime,time_format)
    endemittime_round = round_minutes(endemittime_hhmm, 'down', 15)  # arrotonda per difettos endemittime



runtime=endemittime_round-starttime_round

n_runs = np.int(np.ceil( runtime.total_seconds() / deltat_plumemom ) ) # total number of plumemom runs

if deltat_plumemom >= 3600: # deltat_hhhh cannot be lesss than 3600
    
    deltat_hhhh = deltat_plumemom

else:

    deltat_hhhh = 3600


d = datetime.datetime(2000,1,1) + datetime.timedelta(seconds=deltat_plumemom)
duration_hhmm = str(d.strftime("%H%M"))

duration_h=(int(d.strftime("%H%M")[0:2])+(int(d.strftime("%H%M")[2:4])/float(60)))

d2 = datetime.datetime(2000,1,1) + datetime.timedelta(seconds=deltat_hhhh)

duration_hhhh = '{0:04}'.format(int(str(d2.strftime("%H"))))

#number of lines of each emittimes block

hours = []

startime_count_0 = starttime_hhmm
print ( starttime_hhmm )
startime_count_hour0 = starttime_hhmm.replace(minute=0, second=0, microsecond=0) 
hours.append(startime_count_hour0)

for i in range(n_runs-1):

    startime_count = startime_count_0 + datetime.timedelta(seconds=deltat_plumemom)
    startime_count_hour = startime_count.replace(minute=0, second=0, microsecond=0) 
    hours.append(startime_count_hour)
    startime_count_0 = startime_count

dict_h = Counter(hours)
maximum = max(dict_h, key=dict_h.get)  
num_occurrence = int(dict_h[maximum]) # count the number of plumemom runs within each hour

# particle diameters phi scale
diam_phi = np.linspace(phi_min+(delta_phi),phi_min+(delta_phi)+((n_sections-1)*delta_phi),n_sections, endpoint=True)

# diameter in millimeters [mm]
diam = 2.0**(-np.asarray(diam_phi))

# density in g/cc (calc density compute it in kg/m^3)
# the function calc_density convert internally from mm to m
density = calc_density(diam_phi)/1000

# we assume that all the particles have the same shapefactor
shapefactor = np.ones(npart)*shapefactor

particles_settling_velocity = []

gi = 9.81
visc_atm = 1.8e-5
rho_atm = 1.2

for i in range(npart):

    for j in range(n_sections):

        diam_mt = diam[j]/1000.0

        k1 = shapefactor[i]**(-0.828)
        k2 = 2.0 * np.sqrt( 1.07 - shapefactor[i] )

        # print ( 'k1,k2',k1,k2 )

        mass = density[i,j]*1000.0 * 4.0/3.0 * np.pi * ( 0.5*diam_mt )**3

        # print ( 'mass',i,diam_mt,mass )

        A_cs = np.pi * ( 0.5*diam_mt )**2

        c0 = -2.0 * diam_mt * mass * gi
        c1 = 24.0 * visc_atm * k1 * A_cs
        c2 = rho_atm * diam_mt * k2 * A_cs

        sqrt_delta = np.sqrt( c1**2 - 4.0*c0*c2 )

        Us_1 = ( - c1 + sqrt_delta ) / ( 2 * c2 )
        Us_2 = ( - c1 - sqrt_delta ) / ( 2 * c2 )


        Cd_100 = 24.0/100.0 * k1 + k2
        Us_100 = np.sqrt( 2.0 * mass * gi / ( Cd_100*rho_atm * A_cs ) )

        Cd_1000 = 1.0
        Us_1000 = np.sqrt( 2.0 * mass * gi / ( Cd_1000*rho_atm * A_cs ) )

        Rey1 = rho_atm * diam_mt * Us_1 / visc_atm
        Rey2 = rho_atm * diam_mt * Us_2 / visc_atm

        #print ( 'rho_atm , diam_mt , Us_1 , visc_atm',rho_atm , diam_mt , Us_1 , visc_atm )
        #print ( 'Rey1,Rey2',Rey1,Rey2 )
        
        # Initialization only
        Us = Us_1000

        if ( ( Rey1 >= 0.0 ) and ( Rey1 <= 100.0 ) ):

           # For small Reynolds numbers the drag coefficient is given by Eq.8
           # of Pfeiffer et al. 2005 and the settling velocity is Us_1

           Us = Us_1  

        elif ( ( Rey1 > 100.0 ) and ( Rey1 <= 1000.0 ) ):

           # For intermediate Reyonlds numbers, 100<Re<1000, the drag coefficient 
           # is linearly interpolated between Cd_100 and Cd_1000

           Cd_interp = Cd_100 + ( Rey1 - 100.0 ) / ( 1000.0 - 100.0 ) * ( Cd_1000 - Cd_100)
           Us = np.sqrt( 2.0 * mass * gi / ( Cd_interp * rho_atm * A_cs ) )

        elif ( Rey1 > 1000.0 ):

            # For large Reynolds numbers the drag coefficient is taken as Cd=1,
            # as in Pfeiffer et al. 2005 with the settling velocity is Us_1000

            Us = Us_1000


        if ( ( Rey2 >= 0.0 ) and ( Rey2 <= 100.0 ) ): 

            Us = Us_2

        elif ( ( Rey2 > 100.0 ) and ( Rey2 <= 1000.0 ) ): 

            Cd_interp = Cd_100 + ( Rey2 - 100 ) / ( 1000 - 100 ) * ( Cd_1000 - Cd_100)
            Us = np.sqrt( 2 * mass * gi / ( Cd_interp * rho_atm * A_cs ) )

        elif ( Rey2 > 1000.0 ):

            Us = Us_1000


        particles_settling_velocity.append(Us)
        print ( diam_mt,Us )

released_mass=0


with open('meteo_ground_elev.txt', 'r') as f:
    for line in f:
        z_ground = float(line)     

 
print ( 'EMITTIMES: z_ground',z_ground )

#***************CREATE EMITTIMES and CONTROL FOR PARTICLES DISPERSION************************

with open('EMITTIMES.part','w') as emittimes:    
	emittimes.write('YYYY MM DD HH    DURATION(hhhh) #RECORDS \nYYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w) \n')
emittimes.close()

# search for the maximum number of lines in the .hy files
max_lines = 0

for i in range(n_runs):

    plume_hy = runname + '_{0:03}'.format(i+1)+'.hy'
  
    if os.path.isfile(str(plume_hy)):

        with open(plume_hy) as f:
            max_lines = max(max_lines,sum(1 for _ in f)-1)



"""

First EMITTIMES.part Block

"""

# name of the .hy file
plume_hy = runname + '_{0:03}'.format(1)+'.hy'


# time of the block
#timei =  datetime.datetime.strptime(starttime,time_format)
#timei_end =  starttime_round+datetime.timedelta(seconds=deltat_plumemom)

if n_runs == 1:

    timei =  datetime.datetime.strptime(starttime,time_format)
    timei_end =  datetime.datetime.strptime(endemittime,time_format)

    d = datetime.datetime(2000,1,1) + (timei_end-timei)
    duration_hhmm = str(d.strftime("%H%M"))

else:

    timei =  datetime.datetime.strptime(starttime,time_format)
    timei_end =  starttime_round+datetime.timedelta(seconds=deltat_plumemom)

    d = datetime.datetime(2000,1,1) + (timei_end-timei)
    duration_hhmm = str(d.strftime("%H%M"))

timei_old = timei

print ( 'Block 1',duration_hhmm )

timei_str = timei.strftime("%Y %m %d %H")
timei_str_mm = timei.strftime("%Y %m %d %H %M")

if os.path.isfile(str(plume_hy)):

    data=np.loadtxt(plume_hy,skiprows=1)
    data=data.reshape((-1,int(3+(npart*n_sections))))    

    # data1: array containing data from .hy file, without x,z,h
    data1=np.delete(data, [0,1,2], 1)
    data1 = np.flip(data1, axis=1)

    # array containing lat,lon, height and emission area for time i
    b=[]

    for i0 in range(len(data)):
        x=data[i0,0] #[m]
        y=data[i0,1] #[m] 
        # height=data[i0,2]-vent_height #[m] 	
        height=data[i0,2]-z_ground #[m]	
        emission_area = 0 #[m2]

        # convert from m to lat lon	  
        lon_col = vent_lon + ((x*10**-3)/float(100))
        lat_col = vent_lat + ((y*10**-3)/float(100))

        b.append([lat_col, lon_col, height, emission_area])

    b = np.asarray(b)
    b = b.reshape((-1,4))

    if umbrella_flag == str("T") :

        umbrella_file = runname + '_{0:03}'.format(1)+'.swu'
 
        a = np.loadtxt(umbrella_file,skiprows=1)
        a = np.asarray(a)

        lat_new = vent_lat + ((a[1]*10**-3)/float(100))
        lon_new = vent_lon + ((a[0]*10**-3)/float(100))
        emission_area_new = np.pi * a[2]**(2)
        h_avg = a[3]
        height_new = b[-1,2] + h_avg/float(2)

        b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]	# at nbl replace values for umbrella cloud 	

    # add lines in order to have all the blocks with the same lenght

    for i in range(max_lines-len(b)):

        b = np.vstack(( b , b[len(b)-1,:] + [0.01,0.01,100,0] ))

        data1 = np.vstack((data1,np.zeros((npart*n_sections))))

    # b1 is an array containing lat,lon and height for time i repeated npart*n_sections times
    b1=[]

    for i0 in range(len(b)):    
        for i1 in range(npart):
            for i2 in range(n_sections):
                b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3]])

    b1=np.asarray(b1)
    b1=b1.reshape((-1,4))
	

    # data3 is the array to be written in EMITTIMES for every time interval
    data3 = np.zeros((max_lines*npart*n_sections,5))

    for i0 in range(max_lines):

        i01 = i0*npart*n_sections

        for i1 in range(npart):

            for i2 in range(n_sections):

               data3[i01+(i1*n_sections)+i2,0:4] = b1[i01+(i1*n_sections)+i2,0:4]

               data3[i01+(i1*n_sections)+i2,4] = data1[i0,(i1*n_sections)+i2]

    # mass released in one hour [kg]
    emission_rate = data3[:,4]*3600

    # released_mass_i: mass [kg] released during the simulation at i run time
    released_mass_i=np.sum(emission_rate*duration_h)

    released_mass=released_mass+released_mass_i

    with open('EMITTIMES.part','a') as emittimes:	

        emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence)+'\n')

        print(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence))	

        for h in range(len(data3)):
            emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                   str(data3[h,0]) + ' ' + str(data3[h,1]) + ' ' +
                   str(data3[h,2]) + ' ' + str(emission_rate[h]) +
                   ' '+str(data3[h,3])+' 0.0\n')


else:

    data3 = np.zeros((max_lines*npart*n_sections,5))

    with open('EMITTIMES.part','a') as emittimes:

        emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence)+'\n')

        for h in range(len(data3)):
            emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                   str(vent_lat) + ' ' + str(vent_lon) + ' ' +
                   str(vent_height) + ' ' + str("0") +
                   ' 0.0 0.0\n')

"""

Central EMITTIMES.part Blocks

"""

# loop over the .hy files to write the blocks in EMITTIMES
for i in range(2,n_runs,1):

    

    # name of the .hy file
    plume_hy = runname + '_{0:03}'.format(i)+'.hy'

    # time of the block
    timei =  starttime_round+datetime.timedelta(seconds=(i-1)*deltat_plumemom)    	
    timei_end =  starttime_round+datetime.timedelta(seconds=(i)*deltat_plumemom)



    d = datetime.datetime(2000,1,1) + ( min(endemittime_hhmm,timei_end) - timei )

    duration_hhmm = str(d.strftime("%H%M"))

    print ( 'Block',i,duration_hhmm )

    timei_str = timei.strftime("%Y %m %d %H")
    timei_str_mm = timei.strftime("%Y %m %d %H %M")

    timei_end_str = timei_end.strftime("%Y %m %d %H")

    if os.path.isfile(str(plume_hy)):

        # read the whole plumemom .hy file
        with open(plume_hy, 'r') as fin:
            data = fin.read().splitlines(True)
        fin.close()

        # delete the header line and save to temp.hy
        with open('temp.hy', 'w') as fout:
            fout.writelines(data[1:])
        fout.close()

        # load the data from temp.hy
        with open('temp.hy', 'r') as fin:
            data = np.loadtxt(fin)
        fin.close()
        subprocess.call('rm temp.hy', shell=True)
        # put the data in a numpy array
        data=np.asarray(data)
        data=np.loadtxt(plume_hy,skiprows=1)
        data=data.reshape((-1,int(3+(npart*n_sections))))    

        # data1: array containing data from .hy file, without x,z,h
        data1=np.delete(data, [0,1,2], 1)
        data1 = np.flip(data1, axis=1)

        # array containing lat,lon and height for time i
        b=[]

        for i0 in range(len(data)):
            x=data[i0,0] #[m]
            y=data[i0,1] #[m] 
            # height=data[i0,2]-vent_height #[m] 	
            height=data[i0,2]-z_ground #[m]	

            emission_area = 0

            # convert from m to lat lon	  
            lon_col = vent_lon + ((x*10**-3)/float(100))
            lat_col = vent_lat + ((y*10**-3)/float(100))

            b.append([lat_col, lon_col, height,emission_area])

        b = np.asarray(b)
        b = b.reshape((-1,4))	

        if umbrella_flag == str("T") :

            umbrella_file = runname + '_{0:03}'.format(i)+'.swu' 
 
            a = np.loadtxt(umbrella_file,skiprows=1)
            a = np.asarray(a)

            lat_new = vent_lat + ((a[1]*10**-3)/float(100))
            lon_new = vent_lon + ((a[0]*10**-3)/float(100))
            emission_area_new = np.pi * a[2]**(2)
            h_avg = a[3]
            height_new = b[-1,2] + h_avg/float(2)

            b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]	# at nbl replace values for umbrella cloud 	
  
        # add lines in order to have all the blocks with the same lenght

        for i in range(max_lines-len(b)):

            b = np.vstack(( b , b[len(b)-1,:] + [0.01,0.01,100,0] ))

            data1 = np.vstack((data1,np.zeros((npart*n_sections))))

        # b1 is an array containing lat,lon and height for time i repeated npart*n_sections times
        b1=[]

        for i0 in range(len(b)):    
            for i1 in range(npart):
                for i2 in range(n_sections):
                    b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3]])

        b1=np.asarray(b1)
        b1=b1.reshape((-1,4))	

        # data3 is the array to be written in EMITTIMES for every time interval
        data3 = np.zeros((max_lines*npart*n_sections,5))

        for i0 in range(max_lines):

            i01 = i0*npart*n_sections

            for i1 in range(npart):

                for i2 in range(n_sections):

                    data3[i01+(i1*n_sections)+i2,0:4] = b1[i01+(i1*n_sections)+i2,0:4]

                    data3[i01+(i1*n_sections)+i2,4] = data1[i0,(i1*n_sections)+i2]

	
        # mass released in one hour [kg]
        emission_rate = data3[:,4]*3600
    
        # released_mass_i: mass [kg] released during the simulation at i run time
        released_mass_i=np.sum(emission_rate*duration_h)
    
        released_mass=released_mass+released_mass_i
   
   

        with open('EMITTIMES.part','a') as emittimes:


            if deltat_plumemom >= 3600:

                emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3))+'\n')


            else:


                if timei.hour != timei_old.hour:


                    for hs in dict_h: 
                        if hs == timei_old.hour:
                            num_occurrence_done = int(dict_h[hs])
                            num_occurrence_to_append = num_occurrence - num_occurrence_done
                            data3_to_append = np.zeros((max_lines*npart*n_sections*num_occurrence_to_append,4))

                            timei_str_old = timei_old.strftime("%Y %m %d %H %M")
                        

                            for h1 in range(len(data3_to_append)):

                                emittimes.write(timei_str_old+' '+duration_hhmm+' '+
                                           str(vent_lat) + ' ' + str(vent_lon) +'  '+ 
                                           str(vent_height)+' 0.0 0.0 0.0\n')


                    emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence)+'\n')	


                else:
  
                    pass

            for h in range(len(data3)):
    
                emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                           str(data3[h,0]) + ' ' + str(data3[h,1]) + ' ' +
                           str(data3[h,2]) + ' ' + str(emission_rate[h]) +
                           ' '+str(data3[h,3])+' 0.0\n')

        timei_old = timei


    else:


        data3 = np.zeros((max_lines*npart*n_sections,5)) 

        with open('EMITTIMES.part','a') as emittimes:


            if deltat_plumemom >= 3600:

                emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3))+'\n')


            else:


                if timei.hour != timei_old.hour:
                    for h1 in range(len(data3_to_append)):

                        emittimes.write(timei_str_old+' '+duration_hhmm+' '+
                                           str(vent_lat) + ' ' + str(vent_lon) +'  '+ 
                                           str(vent_height)+' 0.0 0.0 0.0\n')


                    emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence)+'\n')	


                else:
  
                    pass

            for h in range(len(data3)):
    
                emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                           str(vent_lat) + ' ' + str(vent_lon) + ' ' +
                           str(vent_height) + ' ' + str(0) +
                           ' 0.0 0.0\n')

        timei_old = timei


"""

Final EMITTIMES.part Block

"""

if ( n_runs > 1):

    # time of the block
    timei =  starttime_round+datetime.timedelta(seconds=(n_runs-1)*deltat_plumemom)

    #timei =  endemittime_round_down

  
    endemittime_round = round_minutes(endemittime_hhmm, 'up', 60)
    endemittime_round_down = round_minutes(endemittime_hhmm, 'down', 60)

    timei_end = endemittime_round

    d = datetime.datetime(2000,1,1) + (endemittime_hhmm-timei)
    duration_hhmm = str(d.strftime("%H%M"))

    timei_str = timei.strftime("%Y %m %d %H")
    timei_str_mm = timei.strftime("%Y %m %d %H %M")

    print ( 'Block',n_runs,duration_hhmm )

    # name of the .hy file
    plume_hy = runname + '_{0:03}'.format(n_runs)+'.hy'
    if os.path.isfile(str(plume_hy)):

        data=np.loadtxt(plume_hy,skiprows=1)
        data=data.reshape((-1,int(3+(npart*n_sections))))    

        # data1: array containing data from .hy file, without x,z,h
        data1=np.delete(data, [0,1,2], 1)
        data1 = np.flip(data1, axis=1)

        # array containing lat,lon and height for time i
        b=[]

        for i0 in range(len(data)):
            x=data[i0,0] #[m]
            y=data[i0,1] #[m] 
            # height=data[i0,2]-vent_height #[m] 	
            height=data[i0,2]-z_ground #[m]

            emission_area=0	

            # convert from m to lat lon	  
            lon_col = vent_lon + ((x*10**-3)/float(100))
            lat_col = vent_lat + ((y*10**-3)/float(100))

            b.append([lat_col, lon_col, height,emission_area])

        b = np.asarray(b)
        b = b.reshape((-1,4))

        if umbrella_flag == str("T") :

            umbrella_file = runname + '_{0:03}'.format(n_runs)+'.swu'
 
            a = np.loadtxt(umbrella_file,skiprows=1)
            a = np.asarray(a)

            lat_new = vent_lat + ((a[1]*10**-3)/float(100))
            lon_new = vent_lon + ((a[0]*10**-3)/float(100))
            emission_area_new = np.pi * a[2]**(2)
            h_avg = a[3]
            height_new = b[-1,2] + h_avg/float(2)

            b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]	# at nbl replace values for umbrella cloud 	

	
        # add lines in order to have all the blocks with the same lenght

        for i in range(max_lines-len(b)):

            b = np.vstack(( b , b[len(b)-1,:] + [0.01,0.01,100,0] ))

            data1 = np.vstack((data1,np.zeros((npart*n_sections))))

        # b1 is an array containing lat,lon and height for time i repeated npart*n_sections times
        b1=[]

        for i0 in range(len(b)):    
            for i1 in range(npart):
                for i2 in range(n_sections):
                    b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3]])

        b1=np.asarray(b1)
        b1=b1.reshape((-1,4))	

        # data3 is the array to be written in EMITTIMES for every time interval
        data3 = np.zeros((max_lines*npart*n_sections,5))

        for i0 in range(max_lines):

            i01 = i0*npart*n_sections

            for i1 in range(npart):

                for i2 in range(n_sections):

                    data3[i01+(i1*n_sections)+i2,0:4] = b1[i01+(i1*n_sections)+i2,0:4]

                    data3[i01+(i1*n_sections)+i2,4] = data1[i0,(i1*n_sections)+i2]


        # mass released in one hour [kg]
        emission_rate = data3[:,4]*3600

        # released_mass_i: mass [kg] released during the simulation at i run time
        released_mass_i=np.sum(emission_rate*duration_h)

        released_mass=released_mass+released_mass_i

        with open('EMITTIMES.part','a') as emittimes:


            if deltat_plumemom >= 3600:

                emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3))+'\n')


            else:


                if timei.hour != timei_old.hour:

                    for hs in dict_h: 
                        if hs == timei_old.hour:

                            num_occurrence_done = int(dict_h[hs])
                            num_occurrence_to_append = num_occurrence - num_occurrence_done

                            data3_to_append = np.zeros((max_lines*npart*n_sections*num_occurrence_to_append,4))
                            timei_str_old = timei_old.strftime("%Y %m %d %H %M")

                            for h1 in range(len(data3_to_append)):
    
                                emittimes.write(timei_str_old+' '+duration_hhmm+' '+
                                           str(vent_lat) + ' ' + str(vent_lon) +'  '+ 
                                           str(vent_height)+' 0.0 0.0 0.0\n')

                    emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence)+'\n')	


                else:
  
                    pass

            for h in range(len(data3)):

                emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                    str(data3[h,0]) + ' ' + str(data3[h,1]) + ' ' +
                    str(data3[h,2]) + ' ' + str(emission_rate[h]) +
                    ' '+str(data3[h,3])+' 0.0\n')

        timei_old = timei

    else:


        data3 = np.zeros((max_lines*npart*n_sections,4)) 

        with open('EMITTIMES.part','a') as emittimes:


            if deltat_plumemom >= 3600:

                emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3))+'\n')


            else:


                if timei.hour != timei_old.hour:
                    for h1 in range(len(data3_to_append)):

                        emittimes.write(timei_str_old+' '+duration_hhmm+' '+
                                           str(vent_lat) + ' ' + str(vent_lon) +'  '+ 
                                           str(vent_height)+' 0.0 0.0 0.0\n')


                    emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence)+'\n')	


                else:
  
                    pass

            for h in range(len(data3)):
    
                emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                           str(vent_lat) + ' ' + str(vent_lon) + ' ' +
                           str(vent_height) + ' ' + str(0) +
                           ' 0.0 0.0\n')

        timei_old = timei

for hs in dict_h: 
    if hs == timei_old.hour:

        num_occurrence_done = int(dict_h[hs])
        num_occurrence_to_append = num_occurrence - num_occurrence_done

        data3_to_append = np.zeros((max_lines*npart*n_sections*num_occurrence_to_append,4))

        with open('EMITTIMES.part','a') as emittimes:

            for h1 in range(len(data3_to_append)):

                emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                          str(vent_lat) + ' ' + str(vent_lon) + ' ' +
                          str(vent_height)+' 0.0 0.0 0.0\n')



emittimes.close()

# write CONTROL file

if deltat_plumemom >= 3600:

    starttime_round_control = starttime_round.strftime("%Y %m %d %H %M")

else:

    starttime_round_control = starttime_round_c.strftime("%Y %m %d %H %M")

file_control=open('CONTROL.part','w')

file_control.writelines(starttime_round_control+'\n')
file_control.writelines('%d\n'%(max_lines*num_occurrence))
for i in range(max_lines*num_occurrence):
    file_control.writelines("%f %f %f\n"%(vent_lat,vent_lon,vent_height))
file_control.writelines(str(runtime_hh)+'\n')
file_control.writelines('0\n')
file_control.writelines(str(model_top)+'\n')
file_control.writelines('1\n')
file_control.writelines('./\n')
file_control.writelines(meteo_file+'\n')
file_control.writelines('%d\n'%(npart*n_sections))
for i in range(npart):
    i0 = i*n_sections
    for j in range(n_sections):
        file_control.writelines('CL%02d\n'%(i0+j))
        file_control.writelines('0.0\n')
        file_control.writelines('0\n')
        file_control.writelines('00 00 00 00 00\n')
file_control.writelines('1\n')
#file_control.writelines('0.0 0.0\n')
file_control.writelines(str(lat)+' '+str(lon)+'\n')
file_control.writelines(str(spacing_lat)+' '+str(spacing_lon)+'\n')
file_control.writelines(str(span_lat)+' '+str(span_lon)+'\n')
file_control.writelines('./\n')
file_control.writelines('cdump_part_'+runname+'\n')


n_levels = len(H_LEVELS.split())
file_control.writelines(str(n_levels)+'\n')
file_control.writelines(H_LEVELS+'\n')

file_control.writelines(starttime+'\n')
file_control.writelines('00 00 00 00 00\n')
file_control.writelines(str(SI_TYPE)+' '+str(SI_HOUR)+' '+str(SI_MINUTE)+' '+'\n')
file_control.writelines('%d\n'%npart)
for i in range(npart):
    for j in range(n_sections):
        # the diameter should be converted to microns (as required by hysplit) from millimeters
        if ( shapefactor[i] < 1.0 ):
            file_control.writelines('%f %f %f \n'%(1000.0*diam[j],density[i,j],-shapefactor[i]))#50.0 6.0 1.0
        #file_control.writelines('%f %f %f \n'%(0,0,0))#50.0 6.0 1.0
        else:
        #file_control.writelines('%f %f %f \n'%(0,0,0))#50.0 6.0 1.0
            file_control.writelines('%f %f %f \n'%(1000.0*diam[j],density[i,j],shapefactor[i]))#50.0 6.0 1.0        
        # Deposition velocity (m/s), Pollutant molecular weight (Gram/Mole), Surface Reactivity Ratio, Diffusivity  Ratio, Effective Henry's Constant
        file_control.writelines('0.0 0.0 0.0 0.0 0.0 \n')#0 0.0 0.0 0.0 0.0
        # file_control.writelines('1.0 0.0 0.0 0.0 0.0 \n')#0 0.0 0.0 0.0 0.0
        # file_control.writelines(str(particles_settling_velocity[i])+' 0.0 0.0 0.0 0.0 \n')#0 0.0 0.0 0.0 0.0
        file_control.writelines('0.0 0.0 0.0 \n')#0.0 1.0E+06 1.0E-06
        file_control.writelines('0\n')#0
        file_control.writelines('0.0\n')#0.0
file_control.close()



