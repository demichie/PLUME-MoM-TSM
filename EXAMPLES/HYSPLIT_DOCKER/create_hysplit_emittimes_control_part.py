import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
import math as math
from collections import Counter
import pyproj
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import csv

from part_density import calc_density
from input_file import *

from WriteEmittimesFunc import *
from UmbrellaFunc import *

def read_csv_file(file_path):

    data_array = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    data_array=np.delete(data_array,-1,axis=1)

    return data_array



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

n_runs = int(np.ceil( runtime.total_seconds() / deltat_plumemom ) ) # total number of plumemom runs

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
print("Num PM runs for each hour ",num_occurrence)

# particle diameters phi scale
diam_phi = np.linspace(phi_min+(delta_phi),phi_min+(delta_phi)+((n_sections-1)*delta_phi),n_sections, endpoint=True) - 0.5 * delta_phi
#print("diam_phi ",diam_phi)

# diameter in millimeters [mm]
diam = 2.0**(-np.asarray(diam_phi))

# density in g/cc (calc density compute it in kg/m^3)
# the function calc_density convert internally from mm to m
density = calc_density(diam_phi)/1000

if 'shapefactor' in locals():
    print("Constant shape factor")
    shapefactor = np.ones((npart,n_sections))*shapefactor

elif 'shape1' in locals() and 'shape2' in locals():
    print("Linear shape factor (shape1 and shape2)")
    # diam is in millimiters while diam1 and diam2 are in meters

    phi1=np.ones(npart)*phi1
    phi2=np.ones(npart)*phi2

    shape1 = np.ones(npart)*shape1
    shape2 = np.ones(npart)*shape2

    shapefactor = np.zeros((npart,n_sections))

    for i in range(npart):

        for j in range(n_sections):
  
            if ( diam_phi[j] <= phi1[i] ):

                shapefactor[i,j] = shape1[i]

            elif ( phi1[i] < diam_phi[j] < phi2[i] ):
 
                shapefactor[i,j] = shape1[i] + ( diam_phi[j] - phi1[i] ) / ( phi2[i] - phi1[i] ) * ( shape2[i] - shape1[i] )
       
            elif ( diam_phi[j] >= phi2[i] ):

                shapefactor[i,j] = shape2[i]

	#print 'diam',diam[i],diam1[i],diam2[i],diam_phi[i],density[i]
else:
    print("Variable shape factor")
    shapefactor = np.asarray(shape_factor_bin)
    shapefactor = shapefactor.reshape((npart,n_sections))


# shapefactor can be constant or variable with bins
#shapefactor = calc_shapefactor(diam_phi)


particles_settling_velocity = []

gi = 9.81
visc_atm = 1.8e-5
rho_atm = 1.2

for i in range(npart):

    for j in range(n_sections):

        diam_mt = diam[j]/1000.0

        k1 = shapefactor[i,j]**(-0.828)
        k2 = 2.0 * np.sqrt( 1.07 - shapefactor[i,j] )

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
        print ("diam %E mm - SF %5.3f - U_set %6.3f m/s" %(diam[j], shapefactor[i,j], Us) )

released_mass=0


with open('meteo_ground_elev.txt', 'r') as f:
    for line in f:
        z_ground = float(line)     

 
print ( 'EMITTIMES: z_ground',z_ground )

#initialize umbrella radius
r_t = 0

#***************CREATE EMITTIMES and CONTROL FOR PARTICLES DISPERSION************************

#with open('EMITTIMES.temp','w') as emittimes:    
#	emittimes.write('YYYY MM DD HH    DURATION(hhhh) #RECORDS \nYYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w) \n')
#emittimes.close()

num_release_pnts=[]


"""

First EMITTIMES.part Block

"""

# name of the .hy file
plume_hy = runname + '_{0:03}'.format(1)+'_hy.csv'


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

print ( 'Block 1',duration_hhmm, t_sec )

timei_str = timei.strftime("%Y %m %d %H")
timei_str_mm = timei.strftime("%Y %m %d %H %M")

if os.path.isfile(str(plume_hy)):

    data=read_csv_file(plume_hy)
    #data=np.loadtxt(plume_hy,skiprows=1)
    data=data.reshape((-1,int(8+(npart*n_sections))))    

    # data1: array containing data from .hy file, without x,z,h
    data1=np.delete(data, [0,1,2,3,4,5,6,7], 1)
    data1 = np.flip(data1, axis=1)

    # array containing lat,lon, height, emission area and fraction of emission for time i
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

        b.append([lat_col, lon_col, height, emission_area,1])

    b = np.asarray(b)
    b = b.reshape((-1,5))

    if umbrella_flag == str("Model") :

        umbrella_file = runname + '_{0:03}'.format(1)+'.swu'
 
        a = np.loadtxt(umbrella_file,skiprows=1)
        a = np.asarray(a)

        lat_new = vent_lat + ((a[1]*10**-3)/float(100))
        lon_new = vent_lon + ((a[0]*10**-3)/float(100))
        emission_area_new = np.pi * a[2]**(2)
        h_avg = a[3]
        height_new = b[-1,2] + h_avg/float(2)

        b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]	# at nbl replace values for umbrella cloud 	

    elif umbrella_flag == str("Fit"):   
    
        lat_new,lon_new,height_new,emission_area_new,r_t,t_sec,T_sec=UmbrellaFitting(runname,data,b,t_sec,
        T_sec,vent_lat,vent_lon,1)
        
        r_limit = (spacing_met_grid * 1000)/2
        print("Distance [m] computed from the meteo grid ",r_limit)
                
        if small_sources_flag:
        
            if  r_t < r_limit:
                b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]
                
                
            else: 
                  
                r_sm = r_limit   
                utm_x, utm_y, utm_zone = convert_lat_lon_to_utm(lat_new, lon_new)
                x_c = utm_x
                y_c = utm_y           
                emission_area_sm=np.pi*r_sm**2
                x_coords, y_coords, source = generate_small_sources(r_sm,x_c,y_c,r_t,gaussian_source)
                circle_data=[]
                circle_data.append([lat_new, lon_new,r_t,"red"])
                bnew=[]
                for index in range(len(x_coords)):
                    lat_pnt, lon_pnt = convert_utm_to_lat_lon(x_coords[index], y_coords[index], utm_zone)
                    circle_data.append([lat_pnt, lon_pnt, r_sm*(1/3),"blue"])
                    sourceTot=sum(source)
                    sourceNorm=source[index]/sourceTot
                    bnew.append([lat_pnt, lon_pnt,height_new,emission_area_sm,sourceNorm])
                
                if plot_fig:              
                    fig = create_circles_map(circle_data,lat_new, lon_new,r_t,runname + '_{0:03}'.format(1))
                    
                b = np.delete(b, -1, 0)
                b = np.vstack((b, bnew))

                additional_data1_rows = np.tile(data1[-1], (len(x_coords)-1, 1))
                data1 = np.vstack((data1, additional_data1_rows))                                                       
        else:

            b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]	
          
    elif umbrella_flag == str("False"):
        print("No umbrella")

    # b1 is an array containing lat,lon and height for time i repeated npart*n_sections times
    b1=[]
    num_release_pnts.append([len(b)])  

    for i0 in range(len(b)):    
        for i1 in range(npart):
            for i2 in range(n_sections):
                b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3],b[i0,4]])

    b1=np.asarray(b1)
    b1=b1.reshape((-1,5))
    	

    # data3 is the array to be written in EMITTIMES for every time interval
    data3 = np.zeros((len(b)*npart*n_sections,5))
    for i0 in range(len(b)):

        i01 = i0*npart*n_sections

        for i1 in range(npart):

            for i2 in range(n_sections):

               data3[i01+(i1*n_sections)+i2,0:4] = b1[i01+(i1*n_sections)+i2,0:4]

               data3[i01+(i1*n_sections)+i2,4] = data1[i0,(i1*n_sections)+i2] * b1[i01+(i1*n_sections)+i2,4]

    # mass released in one hour [kg]
    emission_rate = data3[:,4]*3600

    # released_mass_i: mass [kg] released during the simulation at i run time
    released_mass_i=np.sum(emission_rate*duration_h)

    released_mass=released_mass+released_mass_i

   
    with open('EMITTIMES.temp','w') as emittimes:	

        emittimes.write(timei_str+' '+duration_hhhh+' XXX\n')

        for h in range(len(data3)):
            emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                   str(data3[h,0]) + ' ' + str(data3[h,1]) + ' ' +
                   str(data3[h,2]) + ' ' + str(emission_rate[h]) +
                   ' '+str(data3[h,3])+' 0.0\n')

else:

    pass


"""

Central EMITTIMES.part Blocks

"""

# loop over the .hy files to write the blocks in EMITTIMES
for i in range(2,n_runs,1):

    

    # name of the .hy file
    plume_hy = runname + '_{0:03}'.format(i)+'_hy.csv'

    # time of the block
    timei =  starttime_round+datetime.timedelta(seconds=(i-1)*deltat_plumemom)    	
    timei_end =  starttime_round+datetime.timedelta(seconds=(i)*deltat_plumemom)



    d = datetime.datetime(2000,1,1) + ( min(endemittime_hhmm,timei_end) - timei )
    t_sec = (d - datetime.datetime(2000,1,1)).total_seconds()

    duration_hhmm = str(d.strftime("%H%M"))

    print ( 'Block',i,duration_hhmm,t_sec)

    timei_str = timei.strftime("%Y %m %d %H")
    timei_str_mm = timei.strftime("%Y %m %d %H %M")

    timei_end_str = timei_end.strftime("%Y %m %d %H")

    if os.path.isfile(str(plume_hy)):

        data=read_csv_file(plume_hy)        
        #data=np.loadtxt(plume_hy,skiprows=1)
        data=data.reshape((-1,int(8+(npart*n_sections))))    

        # data1: array containing data from .hy file, without x,z,h
        data1=np.delete(data, [0,1,2,3,4,5,6,7], 1)
        data1 = np.flip(data1, axis=1)

        # array containing lat,lon and height and fraction for time i
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

            b.append([lat_col, lon_col, height,emission_area,1])

        b = np.asarray(b)
        b = b.reshape((-1,5))	

        if umbrella_flag == str("Model") :

            umbrella_file = runname + '_{0:03}'.format(i)+'.swu' 
 
            a = np.loadtxt(umbrella_file,skiprows=1)
            a = np.asarray(a)

            lat_new = vent_lat + ((a[1]*10**-3)/float(100))
            lon_new = vent_lon + ((a[0]*10**-3)/float(100))
            emission_area_new = np.pi * a[2]**(2)
            h_avg = a[3]
            height_new = b[-1,2] + h_avg/float(2)

            b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]	# at nbl replace values for umbrella cloud

        elif umbrella_flag == str("Fit"):

            lat_new,lon_new,height_new,emission_area_new,r_t,t_sec,T_sec=UmbrellaFitting(runname,data,b,t_sec,
            T_sec,vent_lat,vent_lon,i)
            
       
            if small_sources_flag:
            
                if r_t < r_limit:
                
                    b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]
                    
                else:    
                
                    r_sm = r_limit                   
                    utm_x, utm_y, utm_zone = convert_lat_lon_to_utm(lat_new, lon_new)
                    x_c = utm_x
                    y_c = utm_y           
                    emission_area_sm=np.pi*r_sm**2
                    x_coords, y_coords, source = generate_small_sources(r_sm,x_c,y_c,r_t,gaussian_source)
                    circle_data=[]
                    circle_data.append([lat_new, lon_new,r_t,"red"])
                    
                    bnew=[]
                    for index in range(len(x_coords)):
                        lat_pnt, lon_pnt = convert_utm_to_lat_lon(x_coords[index], y_coords[index], utm_zone)
                        circle_data.append([lat_pnt, lon_pnt, r_sm*(1/3),"blue"])
                        sourceTot=sum(source)
                        sourceNorm=source[index]/sourceTot
                        bnew.append([lat_pnt, lon_pnt,height_new,emission_area_sm,sourceNorm])
                    
                    if plot_fig:          
                        fig = create_circles_map(circle_data,lat_new, lon_new,r_t,runname + '_{0:03}'.format(i))
                        
                    b = np.delete(b, -1, 0)
                    b = np.vstack((b, bnew))

                    additional_data1_rows = np.tile(data1[-1], (len(x_coords)-1, 1))
                    data1 = np.vstack((data1, additional_data1_rows))                                                       
            else:
 
                b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]	


        elif umbrella_flag == str("False"):

            print("No umbrella")

  
        # b1 is an array containing lat,lon and height for time i repeated npart*n_sections times
        b1=[]
        num_release_pnts.append([len(b)])  

        for i0 in range(len(b)):    
            for i1 in range(npart):
                for i2 in range(n_sections):
                    b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3],b[i0,4]])

        b1=np.asarray(b1)
        b1=b1.reshape((-1,5))	

        # data3 is the array to be written in EMITTIMES for every time interval
        data3 = np.zeros((len(b)*npart*n_sections,5))

        for i0 in range(len(b)):

            i01 = i0*npart*n_sections

            for i1 in range(npart):

                for i2 in range(n_sections):

                    data3[i01+(i1*n_sections)+i2,0:4] = b1[i01+(i1*n_sections)+i2,0:4]
                  
                    data3[i01+(i1*n_sections)+i2,4] = data1[i0,(i1*n_sections)+i2] * b1[i01+(i1*n_sections)+i2,4]

	
        # mass released in one hour [kg]
        emission_rate = data3[:,4]*3600
    
        # released_mass_i: mass [kg] released during the simulation at i run time
        released_mass_i=np.sum(emission_rate*duration_h)
    
        released_mass=released_mass+released_mass_i
   


        with open('EMITTIMES.temp','a') as emittimes:


            if deltat_plumemom >= 0:

                emittimes.write(timei_str+' '+duration_hhhh+' XXX\n')
                

            else:

                if timei.hour != timei_old.hour:

                    for hs in dict_h: 
                        if hs.hour == timei_old.hour:
                            num_occurrence_done = int(dict_h[hs])
                            num_occurrence_to_append = num_occurrence - num_occurrence_done
                            data3_to_append = np.zeros((len(b)*npart*n_sections*num_occurrence_to_append,4))

                            timei_str_old = timei_old.strftime("%Y %m %d %H %M")
                        

                            for h1 in range(len(data3_to_append)):

                                emittimes.write(timei_str_old+' '+duration_hhmm+' '+
                                           str(vent_lat) + ' ' + str(vent_lon) +'  '+ 
                                           str(vent_height)+' 0.0 0.0 0.0\n')


                    emittimes.write(timei_str+' '+duration_hhhh+' XXX\n')	
                   

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
    t_sec = (d - datetime.datetime(2000,1,1)).total_seconds()

    duration_hhmm = str(d.strftime("%H%M"))

    timei_str = timei.strftime("%Y %m %d %H")
    timei_str_mm = timei.strftime("%Y %m %d %H %M")

    print ( 'Block',n_runs,duration_hhmm,t_sec)

    # name of the .hy file
    plume_hy = runname + '_{0:03}'.format(n_runs)+'_hy.csv'
    if os.path.isfile(str(plume_hy)):

        #data=np.loadtxt(plume_hy,skiprows=1)
        data=read_csv_file(plume_hy)
        data=data.reshape((-1,int(8+(npart*n_sections))))    

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

            b.append([lat_col, lon_col, height,emission_area,1])

        b = np.asarray(b)
        b = b.reshape((-1,5))

        if umbrella_flag == str("Model") :

            r_old_nbl = data[-1,3] # radius at NBL m
            x_old_nbl = data[-1,0]
            y_old_nbl = data[-1,1]

            x_old_vent = data[0,0]
            y_old_vent = data[0,1]

            umbrella_file = runname + '_{0:03}'.format(n_runs)+'.swu'
 
            a = np.loadtxt(umbrella_file,skiprows=1)
            a = np.asarray(a)

            lat_new = vent_lat + ((a[1]*10**-3)/float(100))
            lon_new = vent_lon + ((a[0]*10**-3)/float(100))
            emission_area_new = np.pi * a[2]**(2)
            h_avg = a[3]
            height_new = b[-1,2] + h_avg/float(2)

            b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]	# at nbl replace values for umbrella cloud 

        elif umbrella_flag == str("Fit"):

            lat_new,lon_new,height_new,emission_area_new,r_t,t_sec,T_sec=UmbrellaFitting(runname,data,b,t_sec,
        T_sec,vent_lat,vent_lon,n_runs)            
            
            
            if small_sources_flag:
            
                if r_t < r_limit:
                
                    b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]
                    
                else:  
                
                    r_sm = r_limit                     
                    utm_x, utm_y, utm_zone = convert_lat_lon_to_utm(lat_new, lon_new)
                    x_c = utm_x
                    y_c = utm_y           
                    emission_area_sm=np.pi*r_sm**2
                    x_coords, y_coords, source = generate_small_sources(r_sm,x_c,y_c,r_t,gaussian_source)
                    circle_data=[]
                    circle_data.append([lat_new, lon_new,r_t,"red"])
                    
                    bnew=[]
                    for index in range(len(x_coords)):
                        lat_pnt, lon_pnt = convert_utm_to_lat_lon(x_coords[index], y_coords[index], utm_zone)
                        circle_data.append([lat_pnt, lon_pnt, r_sm*(1/3),"blue"])
                        sourceTot=sum(source)
                        sourceNorm=source[index]/sourceTot
                        bnew.append([lat_pnt, lon_pnt,height_new,emission_area_sm,sourceNorm])
                    
                    if plot_fig:          
                        fig = create_circles_map(circle_data,lat_new, lon_new,r_t,runname + '_{0:03}'.format(n_runs))

                    b = np.delete(b, -1, 0)
                    b = np.vstack((b, bnew))

                    additional_data1_rows = np.tile(data1[-1], (len(x_coords)-1, 1))
                    data1 = np.vstack((data1, additional_data1_rows))                                                       
            else:
 
                b[-1,[0,1,2,3,4]] = [lat_new,lon_new,height_new,emission_area_new,1]	
            
            
        elif umbrella_flag == str("False"):

            print("No umbrella")	
	

        # b1 is an array containing lat,lon and height for time i repeated npart*n_sections times
        b1=[]
        num_release_pnts.append([len(b)]) 

        for i0 in range(len(b)):    
            for i1 in range(npart):
                for i2 in range(n_sections):
                    b1.append([b[i0,0],b[i0,1],b[i0,2],b[i0,3],b[i0,4]])

        b1=np.asarray(b1)
        b1=b1.reshape((-1,5))	

        # data3 is the array to be written in EMITTIMES for every time interval
        data3 = np.zeros((len(b)*npart*n_sections,5))

        for i0 in range(len(b)):

            i01 = i0*npart*n_sections

            for i1 in range(npart):

                for i2 in range(n_sections):

                    data3[i01+(i1*n_sections)+i2,0:4] = b1[i01+(i1*n_sections)+i2,0:4]
                    
                    data3[i01+(i1*n_sections)+i2,4] = data1[i0,(i1*n_sections)+i2] * b1[i01+(i1*n_sections)+i2,4]


        # mass released in one hour [kg]
        emission_rate = data3[:,4]*3600

        # released_mass_i: mass [kg] released during the simulation at i run time
        released_mass_i=np.sum(emission_rate*duration_h)

        released_mass=released_mass+released_mass_i
        


        with open('EMITTIMES.temp','a') as emittimes:


            if deltat_plumemom >= 0:

                emittimes.write(timei_str+' '+duration_hhhh+' XXX\n')
                

            else:


                if timei.hour != timei_old.hour:
                    for hs in dict_h: 
                        if hs.hour == timei_old.hour:
                            num_occurrence_done = int(dict_h[hs])
                            num_occurrence_to_append = num_occurrence - num_occurrence_done

                            data3_to_append = np.zeros((len(b)*npart*n_sections*num_occurrence_to_append,4))
                            timei_str_old = timei_old.strftime("%Y %m %d %H %M")

                            for h1 in range(len(data3_to_append)):
    
                                emittimes.write(timei_str_old+' '+duration_hhmm+' '+
                                           str(vent_lat) + ' ' + str(vent_lon) +'  '+ 
                                           str(vent_height)+' 0.0 0.0 0.0\n')

                    emittimes.write(timei_str+' '+duration_hhhh+' XXX\n')	
                    

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

        data3_to_append = np.zeros((len(b)*npart*n_sections*num_occurrence_to_append,4))

        with open('EMITTIMES.temp','a') as emittimes:

            for h1 in range(len(data3_to_append)):

                emittimes.write(timei_str_mm+' '+duration_hhmm+' '+
                          str(vent_lat) + ' ' + str(vent_lon) + ' ' +
                          str(vent_height)+' 0.0 0.0 0.0\n')
emittimes.close()

num_release_pnts=np.asarray(num_release_pnts)
num_release_pnts=num_release_pnts.reshape((-1,1))
max_num_release_pnts=np.max(num_release_pnts)
print("num release pnts ",num_release_pnts)

# Example usage
file_path = 'EMITTIMES.temp'  # Replace with your file path
output_file_path = 'EMITTIMES.part'
target_string = 'XXX'  # Replace with the specific string you are looking for
target_element_count = 6

line_numbers=find_line_numbers(file_path, target_string)
number_of_lines = compute_number_of_lines(line_numbers)
max_number_of_lines=max(number_of_lines)
process_file(file_path, target_string, number_of_lines)

with open(output_file_path,'w') as output_file:    
	output_file.write('YYYY MM DD HH    DURATION(hhhh) #RECORDS \nYYYY MM DD HH MM DURATION(hhmm) LAT LON HGT(m) RATE(/h) AREA(m2) HEAT(w) \n')
output_file.close()

multFact=(3600/deltat_plumemom)
print("mult Fact ",multFact)
write_emittimes(file_path,output_file_path, target_element_count,max_number_of_lines,vent_lat,vent_lon,vent_height,multFact)
output_file.close()




# write CONTROL file

if deltat_plumemom >= 3600:

    starttime_round_control = starttime_round.strftime("%Y %m %d %H %M")

else:

    starttime_round_control = starttime_round_c.strftime("%Y %m %d %H %M")

file_control=open('CONTROL.part','w')

file_control.writelines(starttime_round_control+'\n')
file_control.writelines('%d\n'%(max_num_release_pnts*num_occurrence))
for i in range(max_num_release_pnts*num_occurrence):
    file_control.writelines("%f %f %f\n"%(vent_lat,vent_lon,vent_height))
file_control.writelines(str(runtime_hh)+'\n')
file_control.writelines('0\n')
file_control.writelines(str(model_top)+'\n')

file_control.writelines(str(len(meteo_file))+'\n')
for meteo_file_i in meteo_file:
    file_control.writelines(str(meteo_file_dir)+'\n')
    file_control.writelines(meteo_file_i+'\n')

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
        if settling_formulation == "ganser":
            file_control.writelines('%f %f %f \n'%(1000.0*diam[j],density[i,j],-shapefactor[i,j]))#50.0 6.0 1.0
        #file_control.writelines('%f %f %f \n'%(0,0,0))#50.0 6.0 1.0
        elif settling_formulation== "stokes":
        #file_control.writelines('%f %f %f \n'%(0,0,0))#50.0 6.0 1.0
            file_control.writelines('%f %f %f \n'%(1000.0*diam[j],density[i,j],shapefactor[i,j]))#50.0 6.0 1.0        
        # Deposition velocity (m/s), Pollutant molecular weight (Gram/Mole), Surface Reactivity Ratio, Diffusivity  Ratio, Effective Henry's Constant
        file_control.writelines('0.0 0.0 0.0 0.0 0.0 \n')#0 0.0 0.0 0.0 0.0
        # file_control.writelines('1.0 0.0 0.0 0.0 0.0 \n')#0 0.0 0.0 0.0 0.0
        # file_control.writelines(str(particles_settling_velocity[i])+' 0.0 0.0 0.0 0.0 \n')#0 0.0 0.0 0.0 0.0
        file_control.writelines('0.0 0.0 0.0 \n')#0.0 1.0E+06 1.0E-06
        file_control.writelines('0\n')#0
        file_control.writelines('0.0\n')#0.0
file_control.close()



