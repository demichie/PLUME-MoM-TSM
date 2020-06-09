import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil

from extract_wind import write_atm
from input_file import *


def round_minutes(dt, direction, resolution):

    if ( dt.minute%resolution == 0 ):

        rounded_time = dt

    else: 

        new_minute = (dt.minute // resolution + (1 if direction == 'up' else 0)) * resolution

        rounded_time = dt + datetime.timedelta(minutes=new_minute - dt.minute)

    return rounded_time


# create a backup of the input file
src = 'input_file.py'
dst = runname+'.bak'

shutil.copyfile(src, dst)

phi1=np.ones(npart)*phi1
phi2=np.ones(npart)*phi2

cp_part=np.ones(npart)*cp_part
rho1 = np.ones(npart)*rho1
rho2 = np.ones(npart)*rho2

shapefactor = np.ones(npart)*shapefactor

mu = np.ones(npart)*mu
sigma = np.ones(npart)*sigma

solid_partial_mass_fraction=np.ones(npart)*solid_partial_mass_fraction

# create a second template with the parameters constant in time
f = open('plume_model.template','r')
filedata = f.read()
f.close()

filedata = filedata.replace("{npart}", str(npart) )

filedata = filedata.replace("{n_sections}", str(n_sections) )

filedata = filedata.replace("{deltaz_release}", str(deltaz_release) )

filedata = filedata.replace("{ncloud}", str(ncloud) )

filedata = filedata.replace("{vent_height}", str(vent_height) )

filedata = filedata.replace("{phi_min}", str(phi_min) )

filedata = filedata.replace("{delta_phi}", str(delta_phi) )

filedata = filedata.replace("{phi1}", ",".join(np.char.mod('%4f', phi1)) )
filedata = filedata.replace("{phi2}", ",".join(np.char.mod('%4f', phi2)) )

filedata = filedata.replace("{rho1}", ",".join(np.char.mod('%4f', rho1)) )
filedata = filedata.replace("{rho2}", ",".join(np.char.mod('%4f', rho2)) )

filedata = filedata.replace("{shapefactor}", ",".join(np.char.mod('%4f', shapefactor)) )

filedata = filedata.replace("{cp_part}", ",".join(np.char.mod('%4f', cp_part)) )

filedata = filedata.replace("{tmix0}", str(tmix0) )

filedata = filedata.replace("{ngas}", str(ngas) )

if distribution == "LOGNORMAL":

    filedata = filedata.replace("{distribution}", '"LOGNORMAL"')
    filedata = filedata.replace("{mu}", ",".join(np.char.mod('%4f', mu)) )
    filedata = filedata.replace("{sigma}", ",".join(np.char.mod('%4f', sigma)) )
    filedata = filedata.replace("{solid_partial_mass_fraction}", ",".join(np.char.mod('%f', solid_partial_mass_fraction)) )

    filedata = filedata.replace("{bin_partial_mass_fraction}", 'NaN')

elif distribution == "BIN":

    filedata = filedata.replace("{distribution}", '"BIN"')
    filedata = filedata.replace("{mu}", "NaN" )
    filedata = filedata.replace("{sigma}", "NaN" )
    filedata = filedata.replace("{solid_partial_mass_fraction}", ",".join(np.char.mod('%f', solid_partial_mass_fraction)) )

    filedata = filedata.replace("{bin_partial_mass_fraction}", ",".join(np.char.mod('%e', bin_partial_mass_fraction)) )


if ngas>0:

    filedata = filedata.replace("{rvolcgas}", ",".join(np.char.mod('%f', rvolcgas)) )
    filedata = filedata.replace("{cpvolcgas}", ",".join(np.char.mod('%f', cpvolcgas)) )
    filedata = filedata.replace("{volcgas_mol_wt}", ",".join(np.char.mod('%f', volcgas_mol_wt)) )
    filedata = filedata.replace("{volcgas_mass_fraction}", ",".join(np.char.mod('%f', volcgas_mass_fraction)) )

else:

    filedata = filedata.replace("{rvolcgas}", "NaN" )
    filedata = filedata.replace("{cpvolcgas}", "NaN" )
    filedata = filedata.replace("{volcgas_mol_wt}", "NaN" )
    filedata = filedata.replace("{volcgas_mass_fraction}", "NaN" )


filedata = filedata.replace("{water_mass_fraction0}", str(water_mass_fraction0) )

if water_flag == str("T"):

    filedata = filedata.replace("{water_flag}", 'T' )
    filedata = filedata.replace("{rho_lw}", str(rho_lw) )
    filedata = filedata.replace("{rho_ice}", str(rho_ice) )
    filedata = filedata.replace("{added_water_temp}", str(added_water_temp) )
    filedata = filedata.replace("{added_water_mass_fraction}", str(added_water_mass_fraction) )

else:

    filedata = filedata.replace("{water_flag}", 'F' )
    filedata = filedata.replace("{rho_lw}", str(rho_lw) )
    filedata = filedata.replace("{rho_ice}", str(rho_ice) )
    filedata = filedata.replace("{added_water_temp}", str(added_water_temp) )
    filedata = filedata.replace("{added_water_mass_fraction}", str(0) )
     

if umbrella_flag == str("T"):

    filedata = filedata.replace("{umbrella_flag}", 'T' )
    filedata = filedata.replace("{t_end}", str(t_end) )
    filedata = filedata.replace("{dt_output}", str(dt_output) )
    filedata = filedata.replace("{c_d}", str(c_d) )
    filedata = filedata.replace("{steady_flag}", str(steady_flag) )

else:

    filedata = filedata.replace("{umbrella_flag}", 'F' )
    filedata = filedata.replace("{t_end}", "NaN" )
    filedata = filedata.replace("{dt_output}", "NaN" )
    filedata = filedata.replace("{c_d}", "NaN" )
    filedata = filedata.replace("{steady_flag}", "NaN" )

time_format = "%y %m %d %H %M"

if deltat_plumemom >=3600:

    starttime_hhmm = datetime.datetime.strptime(starttime,time_format) 
    starttime_round = round_minutes(starttime_hhmm, 'down', 60) # arrotonda per difetto starttime

    endemittime_hhmm = datetime.datetime.strptime(endemittime,time_format)
    endemittime_round = round_minutes(endemittime_hhmm, 'up', 60) # arrotonda per eccesso endemittime

else:

    starttime_hhmm = datetime.datetime.strptime(starttime,time_format) 
    starttime_round = round_minutes(starttime_hhmm, 'down', 15)  # arrotonda per difettos endemittime

    endemittime_hhmm = datetime.datetime.strptime(endemittime,time_format)
    endemittime_round = round_minutes(endemittime_hhmm, 'down', 15)  # arrotonda per difettos endemittime

runtime=endemittime_round-starttime_round


print ( endemittime_round, starttime_round )
print ( runtime,deltat_plumemom )

n_runs = np.int(np.ceil( runtime.total_seconds() / deltat_plumemom ) ) # numero run di PlumeMoM

print ( 'runtime.total_seconds() ',runtime.total_seconds() )

print ( 'n_runs ', n_runs )

if 'plume_height' in locals():

    filedata = filedata.replace("{inversion_flag}", 'T' )
    filedata = filedata.replace("{log10_mfr}", 'NaN' )
    filedata = filedata.replace("{mfr}", 'NaN' )
    filedata = filedata.replace("{vent_radius}", 'NaN' )

    if isinstance(plume_height, list):
        if len(plume_height) == n_runs:
            plume_height = np.ones(n_runs)*plume_height
        else:
            print ( 'ERROR: check plume_height and n_runs: ',len(plume_height),n_runs )
            sys.exit()
    else:
        plume_height = np.ones(n_runs)*plume_height

    if 'mfr' in locals() or 'log10_mfr' in locals()  :
        print ( 'ERROR: not possible to fix log10_mfr, mfr and plume_height' )
        sys.exit()

if 'log10_mfr' in locals():

    filedata = filedata.replace("{inversion_flag}", 'F' )
    filedata = filedata.replace("{plume_height}", 'NaN' )
    filedata = filedata.replace("{mfr}", 'NaN' )
    filedata = filedata.replace("{vent_radius}", 'NaN' )

    if isinstance(log10_mfr, list):
        if len(log10_mfr) == n_runs:
            log10_mfr = np.ones(n_runs)*log10_mfr
        else:
            print ( 'ERROR: check log10_mfr and n_runs: ',len(log10_mfr),n_runs )
            sys.exit()
    else:
        log10_mfr = np.ones(n_runs)*log10_mfr

    if 'mfr' in locals() or 'plume_height' in locals()  :
        print ( 'ERROR: not possible to fix log10_mfr, mfr and plume_height' )
        sys.exit()

    if 'vent_velocity' in locals() and 'vent_radius' in locals()  :
        print ( 'ERROR: not possible to fix log10_mfr, vent_velocity and vent_radius' )
        sys.exit()


if 'mfr' in locals():

    filedata = filedata.replace("{inversion_flag}", 'F' )
    filedata = filedata.replace("{plume_height}", 'NaN' )
    filedata = filedata.replace("{log10_mfr}", 'NaN' )
    filedata = filedata.replace("{vent_radius}", 'NaN' )
    if isinstance(mfr, list):
        if len(mfr) == n_runs:
            mfr = np.ones(n_runs)*mfr
        else:
            print ( 'ERROR: check mrf and n_runs: ',len(mfr),n_runs )
            sys.exit()
    else:
        mfr = np.ones(n_runs)*mfr

    if 'log10_mfr' in locals() or 'plume_height' in locals()  :
        print ( 'ERROR: not possible to fix log10_mfr, mfr and plume_height' )
        sys.exit()

    if 'vent_velocity' in locals() and 'vent_radius' in locals()  :
        print ( 'ERROR: not possible to fix mfr, vent_velocity and vent_radius' )
        sys.exit()

if 'vent_velocity' and 'vent_radius' in locals():

    print("WARNING: Set both vent velocity and vent radius")

    if 'mfr' in locals():
        print("ERROR: not possible to fix mfr, log10_mfr, vent_velocity and vent_radius ")

    if 'log10_mfr' in locals():
        print("ERROR: not possible to fix mfr, log10_mfr, vent_velocity and vent_radius ")
 
    else:
        filedata = filedata.replace("{inversion_flag}", 'F' )
        filedata = filedata.replace("{plume_height}", 'NaN' )
        filedata = filedata.replace("{log10_mfr}", 'NaN' )
        filedata = filedata.replace("{mfr}", 'NaN' )

if 'vent_velocity' in locals():

    if isinstance(vent_velocity, list):
        if len(vent_velocity) == n_runs:
            vent_velocity = np.ones(n_runs)*vent_velocity
        else:
            print ( 'ERROR: check vent_velocity and n_runs: ',len(vent_velocity),n_runs )
            sys.exit()
    else:
        vent_velocity = np.ones(n_runs)*vent_velocity

if 'vent_radius' in locals():

    if isinstance(vent_radius, list):
        if len(vent_radius) == n_runs:
            vent_radius = np.ones(n_runs)*vent_radius
        else:
            print ( 'ERROR: check vent_radius and n_runs: ',len(vent_radius),n_runs )
            sys.exit()
    else:
        vent_radius = np.ones(n_runs)*vent_radius

f = open('plume_model.temp1','w')
f.write(filedata)

f.close()

col=[]

for i in range(n_runs):

    runnamenew = runname + '_{0:03}'.format(i+1)
    print ( '------------>runname',runnamenew )

    f = open('plume_model.temp1','r')
    filedata = f.read()
    f.close()

    run_flag = 0

    # create a third template with the parameters changing with time
    f = open('plume_model.temp2','w')

    filedata = filedata.replace("{runname}", '"'+str(runnamenew)+'"' )


    if 'plume_height' in locals():
        
        if plume_height[i] == 0:

            print ( "*** EMISSION STOP ***" )

            run_flag = 0

        else:

            filedata = filedata.replace("{plume_height}", str(plume_height[i]) )
            run_flag = 1

    if 'log10_mfr' in locals():

        if log10_mfr[i] == 0:
    
            print ( "*** EMISSION STOP ***" )

            run_flag = 0
    
        else:

            filedata = filedata.replace("{log10_mfr}", str(log10_mfr[i]) )
            run_flag = 1

    if 'mfr' in locals():

        if mfr[i] == 0:
    
            print ( "*** EMISSION STOP ***" )

            run_flag = 0
    
        else:

            filedata = filedata.replace("{mfr}", str(mfr[i]) )
            run_flag = 1

    if 'vent_velocity' in locals():

        filedata = filedata.replace("{vent_velocity}", str(vent_velocity[i]) )
        run_flag = 1

    if 'vent_radius' in locals():

        filedata = filedata.replace("{vent_radius}", str(vent_radius[i]) )
        run_flag = 1


    
    f.write(filedata)
    f.close()

    if run_flag == 0:

        pass

    else:

        timei =  datetime.datetime.strptime(starttime,time_format)+datetime.timedelta(seconds=i*deltat_plumemom)
	
        timei_str = timei.strftime("%y %m %d %H %M")

        print ( 'Time',timei_str )

        write_atm(timei_str)


        # append the atmospheric data to the input file
        filenames = ['plume_model.temp2', 'atm.txt']
        with open('plume_model.inp', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    
        subprocess.call(plumemom_dir+"/bin/PLUMEMoM", shell=True) 


#subprocess.call("rm plume_model.temp1", shell=True) 
#subprocess.call("rm plume_model.temp2", shell=True) 






