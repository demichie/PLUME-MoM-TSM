import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
import math as math
from collections import Counter
import math as math

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

released_mass=0


with open('meteo_ground_elev.txt', 'r') as f:
    for line in f:
        z_ground = float(line)     

 
print ( 'EMITTIMES: z_ground',z_ground )

#***************CREATE EMITTIMES and CONTROL FOR GAS DISPERSION************************

with open('EMITTIMES.gas','w') as emittimes:    
	emittimes.write('YYYY MM DD HH    DURATION(hhhh) #RECORDS \nYYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w) \n')
emittimes.close()

# search for the maximum number of lines in the .hy files
max_lines = 0

for i in range(n_runs):

    plume_hy = runname + '_{0:03}'.format(i+1)+'_volcgas.hy'
  
    if os.path.isfile(str(plume_hy)):

        with open(plume_hy) as f:
            max_lines = max(max_lines,sum(1 for _ in f)-1)

"""

First EMITTIMES.gas Block

"""

# name of the .hy file
plume_hy = runname + '_{0:03}'.format(1)+'_volcgas.hy'

# time of the block
#timei =  datetime.datetime.strptime(starttime,time_format)
#timei_end =  starttime_round+datetime.timedelta(seconds=deltat_plumemom)

if n_runs == 1:

    timei =  datetime.datetime.strptime(starttime,time_format)
    timei_end =  datetime.datetime.strptime(endemittime,time_format)

    d = datetime.datetime(2000,1,1) + (timei_end-timei)
    t_sec = (d - datetime.datetime(2000,1,1)).total_seconds()
    duration_hhmm = str(d.strftime("%H%M"))

else:

    timei =  datetime.datetime.strptime(starttime,time_format)
    timei_end =  starttime_round+datetime.timedelta(seconds=deltat_plumemom)

    d = datetime.datetime(2000,1,1) + (timei_end-timei)
    t_sec = (d - datetime.datetime(2000,1,1)).total_seconds()
    duration_hhmm = str(d.strftime("%H%M"))

timei_old = timei

T_sec = 0 # seconds from the beginnig of the emission

print ( 'Block 1',duration_hhmm , t_sec )

timei_str = timei.strftime("%Y %m %d %H")
timei_str_mm = timei.strftime("%Y %m %d %H %M")

if os.path.isfile(str(plume_hy)):

    data=np.loadtxt(plume_hy,skiprows=1)
    data=data.reshape((-1,int(8+(ngas))))    

    # data1: array containing data from .hy file, without x,z,h
    data1=np.delete(data, [0,1,2,3,4,5,6,7], 1)
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

    if umbrella_flag == str("Model") :

        umbrella_file = runname + '_{0:03}'.format(1)+'.swu'
 
        a = np.loadtxt(umbrella_file,skiprows=1)
        a = np.asarray(a)

        lat_new = vent_lat + ((a[1]*10**-3)/float(100))
        lon_new = vent_lon + ((a[0]*10**-3)/float(100))
        emission_area_new = np.pi * a[2]**(2)
        h_avg = a[3]
        height_new = b[-1,2] + h_avg/float(2)

        b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]	# at nbl replace values for umbrella cloud

    elif umbrella_flag == str("Fit") :

        print("** Umbrella Fitting**")
        r_old_nbl = data[-1,3] # radius at NBL m
        x_old_nbl = data[-1,0]
        y_old_nbl = data[-1,1]

        d_old = np.sqrt((x_old_nbl)**2+(y_old_nbl)**2)

        u_atm_nbl = data[-1,4]
        v_atm_nbl = data[-1,5]
        vel_atm_nbl = ( u_atm_nbl**2 + v_atm_nbl**2 )**(0.5)
        rho_mix_nbl = data[-1,6]
        mfr_nbl = data[-1,7]
        vfr_nbl = mfr_nbl /float(rho_mix_nbl)

        a_fit = 4.02*10**(-9)
        b_fit = 14.07
        c_fit = -0.7457

        Delta_d = a_fit * (vel_atm_nbl)**(c_fit)*(np.log10(vfr_nbl))**(b_fit)
        d_new_nbl = d_old + Delta_d

        A_fit = 1.89*10**(-8)
        B_fit = 13.9
        C_fit = -0.8856

        r_new_nbl = A_fit*(vel_atm_nbl)**(C_fit)*(np.log10(vfr_nbl))**(B_fit)

        alpha_fit = 1.9*10**(-5)
        beta_fit = 1.09*10
        gamma_fit = -1.57

        t_steady_fit =  alpha_fit*(vel_atm_nbl)**(gamma_fit)*(np.log10(vfr_nbl))**(beta_fit)

        T_sec += t_sec

        T_mean = int(t_sec/float(2))

        if r_new_nbl==0 and t_steady_fit==0 and Delta_d==0:
            vfr_nbl = 0
            d_new_nbl = 0
            d_old = 0
            f_t = 0
            r_t = 0
            d_t = 0
            lat_new = vent_lat
            lon_new = vent_lon
            height_new = 0
            emission_area_new = 0
        else:
            f_t =  2.0/float(np.pi)*np.arctan(13.6*T_mean/float(t_steady_fit))

            r_t = (1 - f_t ) * r_old_nbl + f_t * r_new_nbl 
            d_t = (1 - f_t ) * d_old + f_t * d_new_nbl

            alpha = math.acos(u_atm_nbl / float(vel_atm_nbl)) #rad
            x_new_nbl = d_t * math.cos(alpha) * np.sign(u_atm_nbl)
            y_new_nbl = d_t * math.sin(alpha) * np.sign(v_atm_nbl)

            lat_new = vent_lat + ((y_new_nbl*10**-3)/float(100))
            lon_new = vent_lon + ((x_new_nbl*10**-3)/float(100))
            height_new = b[-1,2]
            emission_area_new = np.pi * r_t**(2)
       
        b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]

    elif umbrella_flag == str("False") :

        print("No umbrella")	 	

    # add lines in order to have all the blocks with the same lenght

    for i in range(max_lines-len(b)):

        b = np.vstack(( b , b[len(b)-1,:] + [0.01,0.01,100,0] ))

        data1 = np.vstack((data1,np.zeros((ngas))))

    # b1 is an array containing lat,lon and height for time i repeated ngas times
    b1=[]

    for i0 in range(len(b)):    
        for i1 in range(ngas):
            b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3]])

    b1=np.asarray(b1)
    b1=b1.reshape((-1,4))
	

    # data3 is the array to be written in EMITTIMES for every time interval
    data3 = np.zeros((max_lines*ngas,5))

    for i0 in range(max_lines):

        i01 = i0*ngas

        for i1 in range(ngas):

           data3[i01+i1,0:4] = b1[i01+i1,0:4]

           data3[i01+i1,4] = data1[i0,i1]

    # mass released in one hour [kg]
    emission_rate = data3[:,4]*3600

    # released_mass_i: mass [kg] released during the simulation at i run time
    released_mass_i=np.sum(emission_rate*duration_h)

    released_mass=released_mass+released_mass_i

    with open('EMITTIMES.gas','a') as emittimes:	

        emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence)+'\n')	

        for h in range(len(data3)):
            emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                   str(data3[h,0]) + ' ' + str(data3[h,1]) + ' ' +
                   str(data3[h,2]) + ' ' + str(emission_rate[h]) +
                   ' '+str(data3[h,3])+' 0.0\n')


else:

    data3 = np.zeros((max_lines*ngas,5))

    with open('EMITTIMES.gas','a') as emittimes:

        emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3)*num_occurrence)+'\n')

        for h in range(len(data3)):
            emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                   str(vent_lat) + ' ' + str(vent_lon) + ' ' +
                   str(vent_height) + ' ' + str("0") +
                   ' 0.0 0.0\n')

"""

Central EMITTIMES.gas Blocks

"""

# loop over the .hy files to write the blocks in EMITTIMES
for i in range(2,n_runs,1):

    

    # name of the .hy file
    plume_hy = runname + '_{0:03}'.format(i)+'_volcgas.hy'

    # time of the block
    timei =  starttime_round+datetime.timedelta(seconds=(i-1)*deltat_plumemom)    	
    timei_end =  starttime_round+datetime.timedelta(seconds=(i)*deltat_plumemom)

    d = datetime.datetime(2000,1,1) + ( min(endemittime_hhmm,timei_end) - timei )

    duration_hhmm = str(d.strftime("%H%M"))

    print ( 'Block',i,duration_hhmm,t_sec )

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
        data=data.reshape((-1,int(8+(ngas))))    

        # data1: array containing data from .hy file, without x,z,h
        data1=np.delete(data, [0,1,2,3,4,5,6,7], 1)
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

        if umbrella_flag == str("Model") :

            umbrella_file = runname + '_{0:03}'.format(i)+'.swu'
 
            print ( umbrella_file )
 
            a = np.loadtxt(umbrella_file,skiprows=1)
            a = np.asarray(a)

            lat_new = vent_lat + ((a[1]*10**-3)/float(100))
            lon_new = vent_lon + ((a[0]*10**-3)/float(100))
            emission_area_new = np.pi * a[2]**(2)
            h_avg = a[3]
            height_new = b[-1,2] + h_avg/float(2)

            b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]	# at nbl replace values for umbrella cloud 

        elif umbrella_flag == str("Fit") :

            print("** Umbrella Fitting**")

            r_old_nbl = data[-1,3] # radius at NBL m
            x_old_nbl = data[-1,0]
            y_old_nbl = data[-1,1]
 
            u_atm_nbl = data[-1,4]
            v_atm_nbl = data[-1,5]
            vel_atm_nbl = ( u_atm_nbl**2 + v_atm_nbl**2 )**(0.5)
            rho_mix_nbl = data[-1,6]
            mfr_nbl = data[-1,7]
            vfr_nbl = mfr_nbl /float(rho_mix_nbl)

            a_fit = 4.02*10**(-9)
            b_fit = 14.07
            c_fit = -0.7457

            Delta_d = a_fit * (vel_atm_nbl)**(c_fit)*(np.log10(vfr_nbl))**(b_fit)
            d_new_nbl = d_old + Delta_d

            A_fit = 1.89*10**(-8)
            B_fit = 13.9
            C_fit = -0.8856

            r_new_nbl = A_fit*(vel_atm_nbl)**(C_fit)*(np.log10(vfr_nbl))**(B_fit)

            alpha_fit = 1.9*10**(-5)
            beta_fit = 1.09*10
            gamma_fit = -1.57

            t_steady_fit =  alpha_fit*(vel_atm_nbl)**(gamma_fit)*(np.log10(vfr_nbl))**(beta_fit)

            T_mean = T_sec + int(t_sec/float(2))

            T_sec += t_sec

            if r_new_nbl==0 and t_steady_fit==0 and Delta_d==0:
                vfr_nbl = 0
                d_new_nbl = 0
                d_old = 0
                f_t = 0
                r_t = 0
                d_t = 0
                lat_new = vent_lat
                lon_new = vent_lon
                height_new = 0
                emission_area_new = 0

            else:

                f_t =  2.0/float(np.pi)*np.arctan(13.6*T_mean/float(t_steady_fit))

                r_t = (1 - f_t ) * r_old_nbl + f_t * r_new_nbl 
                d_t = (1 - f_t ) * d_old + f_t * d_new_nbl

                alpha = math.acos(u_atm_nbl / float(vel_atm_nbl)) #rad
                x_new_nbl = d_t * math.cos(alpha) * np.sign(u_atm_nbl)
                y_new_nbl = d_t * math.sin(alpha) * np.sign(v_atm_nbl)

                lat_new = vent_lat + ((y_new_nbl*10**-3)/float(100))
                lon_new = vent_lon + ((x_new_nbl*10**-3)/float(100))
                height_new = b[-1,2]
                emission_area_new = np.pi * r_t**(2)

            b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]	

        elif umbrella_flag == str("False") :
 
            print("No umbrella")
  
        # add lines in order to have all the blocks with the same lenght

        for i in range(max_lines-len(b)):

            b = np.vstack(( b , b[len(b)-1,:] + [0.01,0.01,100,0] ))

            data1 = np.vstack((data1,np.zeros((ngas))))

        # b1 is an array containing lat,lon and height for time i repeated ngas times
        b1=[]

        for i0 in range(len(b)):    
            for i1 in range(ngas):
                b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3]])

        b1=np.asarray(b1)
        b1=b1.reshape((-1,4))	

        # data3 is the array to be written in EMITTIMES for every time interval
        data3 = np.zeros((max_lines*ngas,5))

        for i0 in range(max_lines):

            i01 = i0*ngas

            for i1 in range(ngas):

                data3[i01+i1,0:4] = b1[i01+i1,0:4]

                data3[i01+i1,4] = data1[i0,i1]
	
        # mass released in one hour [kg]
        emission_rate = data3[:,4]*3600
    
        # released_mass_i: mass [kg] released during the simulation at i run time
        released_mass_i=np.sum(emission_rate*duration_h)
    
        released_mass=released_mass+released_mass_i 

        with open('EMITTIMES.gas','a') as emittimes:


            if deltat_plumemom >= 3600:

                emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3))+'\n')


            else:


                if timei.hour != timei_old.hour:


                    for hs in dict_h: 
                        if hs.hour == timei_old.hour:
                            num_occurrence_done = int(dict_h[hs])
                            num_occurrence_to_append = num_occurrence - num_occurrence_done
                            data3_to_append = np.zeros((max_lines*ngas*num_occurrence_to_append,4))

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
    
        pass



"""

Final EMITTIMES.gas Block

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

    print ( 'Block',n_runs,duration_hhmm,t_sec )

    # name of the .hy file
    plume_hy = runname + '_{0:03}'.format(n_runs)+'_volcgas.hy'
    if os.path.isfile(str(plume_hy)):

        data=np.loadtxt(plume_hy,skiprows=1)
        data=data.reshape((-1,int(8+(ngas))))    

        # data1: array containing data from .hy file, without x,z,h
        data1=np.delete(data, [0,1,2,3,4,5,6,7], 1)
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

        if umbrella_flag == str("Model") :

            umbrella_file = runname + '_{0:03}'.format(n_runs)+'.swu'
 
            print ( umbrella_file )
 
            a = np.loadtxt(umbrella_file,skiprows=1)
            a = np.asarray(a)

            lat_new = vent_lat + ((a[1]*10**-3)/float(100))
            lon_new = vent_lon + ((a[0]*10**-3)/float(100))
            emission_area_new = np.pi * a[2]**(2)
            h_avg = a[3]
            height_new = b[-1,2] + h_avg/float(2)
  
            b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]	# at nbl replace values for umbrella cloud 	

        elif umbrella_flag == str("Fit") :

            print("** Umbrella Fitting**")
            r_old_nbl = data[-1,3] # radius at NBL m
            x_old_nbl = data[-1,0]
            y_old_nbl = data[-1,1]

            u_atm_nbl = data[-1,4]
            v_atm_nbl = data[-1,5]
            vel_atm_nbl = ( u_atm_nbl**2 + v_atm_nbl**2 )**(0.5)
            rho_mix_nbl = data[-1,6]
            mfr_nbl = data[-1,7]
            vfr_nbl = mfr_nbl /float(rho_mix_nbl)

            a_fit = 4.02*10**(-9)
            b_fit = 14.07
            c_fit = -0.7457

            Delta_d = a_fit * (vel_atm_nbl)**(c_fit)*(np.log10(vfr_nbl))**(b_fit)
            d_new_nbl = d_old + Delta_d

            A_fit = 1.89*10**(-8)
            B_fit = 13.9
            C_fit = -0.8856

            r_new_nbl = A_fit*(vel_atm_nbl)**(C_fit)*(np.log10(vfr_nbl))**(B_fit)

            alpha_fit = 1.9*10**(-5)
            beta_fit = 1.09*10
            gamma_fit = -1.57

            t_steady_fit =  alpha_fit*(vel_atm_nbl)**(gamma_fit)*(np.log10(vfr_nbl))**(beta_fit)

            T_mean = T_sec + int(t_sec/float(2))

            T_sec += t_sec

            if r_new_nbl==0 and t_steady_fit==0 and Delta_d==0:
                vfr_nbl = 0
                d_new_nbl = 0
                d_old = 0
                f_t = 0
                r_t = 0
                d_t = 0
                lat_new = vent_lat
                lon_new = vent_lon
                height_new = 0
                emission_area_new = 0

            else:

                f_t =  2.0/float(np.pi)*np.arctan(13.6*T_mean/float(t_steady_fit))

                r_t = (1 - f_t ) * r_old_nbl + f_t * r_new_nbl 
                d_t = (1 - f_t ) * d_old + f_t * d_new_nbl

                alpha = math.acos(u_atm_nbl / float(vel_atm_nbl)) #rad
                x_new_nbl = d_t * math.cos(alpha) * np.sign(u_atm_nbl)
                y_new_nbl = d_t * math.sin(alpha) * np.sign(v_atm_nbl)

                lat_new = vent_lat + ((y_new_nbl*10**-3)/float(100))
                lon_new = vent_lon + ((x_new_nbl*10**-3)/float(100))
                height_new = b[-1,2]
                emission_area_new = np.pi * r_t**(2)

            b[-1,[0,1,2,3]] = [lat_new,lon_new,height_new,emission_area_new]	

        elif umbrella_flag == str("False") :

            print("No umbrella")

        # add lines in order to have all the blocks with the same lenght

        for i in range(max_lines-len(b)):

            b = np.vstack(( b , b[len(b)-1,:] + [0.01,0.01,100,0] ))

            data1 = np.vstack((data1,np.zeros((ngas))))

        # b1 is an array containing lat,lon and height for time i repeated ngas times
        b1=[]

        for i0 in range(len(b)):    
            for i1 in range(ngas):
                b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3]])

        b1=np.asarray(b1)
        b1=b1.reshape((-1,4))	

        # data3 is the array to be written in EMITTIMES for every time interval
        data3 = np.zeros((max_lines*ngas,5))

        for i0 in range(max_lines):

            i01 = i0*ngas

            for i1 in range(ngas):

                data3[i01+i1,0:4] = b1[i01+i1,0:4]

                data3[i01+i1,4] = data1[i0,i1]


        # mass released in one hour [kg]
        emission_rate = data3[:,4]*3600

        # released_mass_i: mass [kg] released during the simulation at i run time
        released_mass_i=np.sum(emission_rate*duration_h)

        released_mass=released_mass+released_mass_i

        with open('EMITTIMES.gas','a') as emittimes:


            if deltat_plumemom >= 3600:

                emittimes.write(timei_str+' '+duration_hhhh+' '+str(len(data3))+'\n')


            else:


                if timei.hour != timei_old.hour:

                    for hs in dict_h: 
                        if hs.hour == timei_old.hour:

                            num_occurrence_done = int(dict_h[hs])
                            num_occurrence_to_append = num_occurrence - num_occurrence_done

                            data3_to_append = np.zeros((max_lines*ngas*num_occurrence_to_append,4))
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

        pass

for hs in dict_h: 
    if hs.hour == timei_old.hour:

        num_occurrence_done = int(dict_h[hs])
        num_occurrence_to_append = num_occurrence - num_occurrence_done

        data3_to_append = np.zeros((max_lines*ngas*num_occurrence_to_append,4))

        with open('EMITTIMES.gas','a') as emittimes:

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

file_control=open('CONTROL.gas','w')

file_control.writelines(starttime_round_control+'\n')
file_control.writelines('%d\n'%(max_lines*num_occurrence))
for i in range(max_lines*num_occurrence):
    file_control.writelines("%f %f %f\n"%(vent_lat,vent_lon,vent_height))
file_control.writelines(str(runtime_hh)+'\n')
file_control.writelines('0\n')
file_control.writelines(str(model_top)+'\n')

file_control.writelines(str(len(meteo_file))+'\n')
for meteo_file_i in meteo_file:
    file_control.writelines(str(meteo_file_dir)+'\n')
    file_control.writelines(meteo_file_i+'\n')

file_control.writelines('%d\n'%(ngas))
for i in range(ngas):
    file_control.writelines('GS%02d\n'%(i))
    file_control.writelines('0.0\n')
    file_control.writelines('0\n')
    file_control.writelines('00 00 00 00 00\n')
file_control.writelines('1\n')
#file_control.writelines('0.0 0.0\n')
file_control.writelines(str(lat)+' '+str(lon)+'\n')
file_control.writelines(str(spacing_lat)+' '+str(spacing_lon)+'\n')
file_control.writelines(str(span_lat)+' '+str(span_lon)+'\n')
file_control.writelines('./\n')
file_control.writelines('cdump_gas_'+runname+'\n')


n_levels = len(H_LEVELS.split())
file_control.writelines(str(n_levels)+'\n')
file_control.writelines(H_LEVELS+'\n')

file_control.writelines(starttime+'\n')
file_control.writelines('00 00 00 00 00\n')
file_control.writelines(str(SI_TYPE)+' '+str(SI_HOUR)+' '+str(SI_MINUTE)+' '+'\n')
file_control.writelines('%d\n'%ngas)
for i in range(ngas):
    file_control.writelines('0.0 0.0 0.0\n')# 0.0 0.0 0.0 
    # Deposition velocity (m/s), Pollutant molecular weight (Gram/Mole), Surface Reactivity Ratio, Diffusivity  Ratio, Effective Henry's Constant
    file_control.writelines('0.0 0.0 0.0 0.0 0.0 \n')#0 0.0 0.0 0.0 0.0
    file_control.writelines('0.0 0.0 0.0 \n')#0.0 1.0E+06 1.0E-06
    file_control.writelines('0.0\n')#0
    file_control.writelines('0.0\n')#0.0
file_control.close()


