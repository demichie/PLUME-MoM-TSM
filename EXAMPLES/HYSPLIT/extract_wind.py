import numpy as np
import subprocess
import datetime
import os,sys,time
import shutil
from scipy import interpolate 

def calc_atm(profile,z_ground,TMPS,RH2M,prss,fields_list):
    """create atmosperic profile for plumemom """
    nlev=len(profile[:,0])

    print 'nlev',nlev

    grav = 9.81

    # gas constant for dry air (J kg-1 K-1)
    Rd = 287.04

    # gas constant for water vapor (J kg-1 K-1)
    Rw = 461.5

    if ( z_ground*TMPS*RH2M*prss == 0 ):

        # reference pressure (Pa)
        P0 = 1011*100.0
 
        T0 = 25.8+273.15
        RH0 = 25.8
        z0 = 24.0

    else:

        # surface pressure (Pa)
        P0 = prss*100.0

        # surface temperature (K)
        T0 = TMPS+273.15

        # relative humidity at 2m (percentage)
        RH0 = RH2M

        # surface height (m)
        z0 = z_ground

    #print 'P0',P0
    #print 'T0',T0
    #print 'RH0',RH0
    #print 'z0',z0

    Es0 = np.exp(21.4-(5351/T0))
    Q0 = RH0 / P0 * ( 0.622 * Es0 )
    Tv0 = T0 * ( 1 + 0.61*Q0)

    # pressure (Pa)
    P = profile[:,0]*100.0

    # temperature (deg C)
    T_idx = fields_list.index("TEMP") + 1

    # temperature (K)
    T = profile[:,T_idx] + 273.15

    # Saturation mixing ratio (hPa)
    Es = np.exp(21.4-(5351/T))

    # dry air density
    rho_da = P / ( Rd * T )

    if 'SPHU' in fields_list:

        # specific humidity in g/kg
        SPHU_idx = fields_list.index("SPHU") + 1

        # specific humidity in kg/kg
        Q = profile[:,SPHU_idx]/1000.0

	# relative humidity (in percentage)
        RH = (Q * P) / ( 0.622 * Es )

    else:

        RH_idx = fields_list.index("RELH") + 1

        RH = profile[:,RH_idx]

        Q = RH / P * ( 0.622 * Es )

    # Mixture density, Eq. (3) http://www.engineeringtoolbox.com/density-air-d_680.html
    rho = rho_da * ( 1 + Q ) / ( 1 + Q * Rw / Rd )
    # where Q is the specific humidity in kgkg-1

    if 'HGTSa' in fields_list:

        z_idx = fields_list.index("HGTS") + 1

        z = profile[:,z_idx]

    else:

        m = np.min(np.where(P<P0))

        # print 'm', m


        l = m-1
        
        z=np.zeros(nlev)

        # virtual temperature (Eq. 9 arl-224.pdf)
        Tv = T * ( 1 + 0.61*Q)

        if ( m > 0 ):

            # average virtual temperature
            Tvl0 = 0.5*(Tv[l]+Tv0)

            # height increment (Eq.10 arl-224.pdf)
            deltaz0 = np.log(P0/P[l]) *Rd *Tvl0 / grav

            z[l] = z0 + deltaz0

            # print 'z0,z[l],l',z0,z[l],l

            Tv12 = 0.5*(Tv[0:l]+Tv[1:l+1])
            # print 'Tv12',Tv12
            deltaz = np.log(P[1:l+1]/P[0:l]) * Rd * Tv12 / grav

            # print 'deltaz',deltaz

            z[0:l] = z[l] + np.cumsum(deltaz[::-1])[::-1] 

            # print 'z[0:l]',z[0:l+1]


        # average virtual temperature
        Tv0m = 0.5*(Tv0+Tv[m])

        # height increment (Eq.10 arl-224.pdf)
        deltaz0 = np.log(P0/P[m]) *Rd *Tv0m / grav


        z[m] = z0 + deltaz0

        Tv12 = 0.5*(Tv[m:-1]+Tv[m+1:])
        deltaz = np.log(P[m:-1]/P[m+1:]) *Rd *Tv12 / grav

        z[m+1:] = z[m] + np.cumsum(deltaz)

    print 'z',z

    # W->E component of horizontal velocity (m/s)
    U_idx = fields_list.index("UWND") + 1
    U = profile[:,U_idx]
   
    # S->N component of horizontal velocity (m/s)
    V_idx = fields_list.index("VWND") + 1
    V = profile[:,V_idx]

    # Return z (km), density (kg/m3) , P (hPa) , T (K) , specific humidity (kg/kg), WE vel (m/s) , SN vel (m/s) 
    return z/1000.0, rho, P/100, T, Q, U, V 

def write_atm(time_input):
    """create atmosperic profile for plumemom """

    from input_file import hysplit_dir, meteo_file, vent_lat, vent_lon

    run_path = './'

    profile = os.path.join(hysplit_dir,'exec','profile')

    subprocess.call(profile+" -d"+run_path+" -f"+meteo_file+" -y"
                    +str(vent_lat)+" -x"+str(vent_lon)+" -o0 -p00", shell=True) # wind profile at vent position at the beginnig of the wind file

    with open("profile_00.txt","r") as file1:
        line=file1.readlines()
        time_start_line=line[1]
        time_start=time_start_line[19:35]

        # add 0 to year if it is a single digit (ex: 6 for 2006 becomes 06)
        time_split = time_start.split()
        year = time_split[0]
        if len(year)==1:
            time_start = '0'+time_start[1:14]

    file1.close()

    for i in range(len(line)):
        if ( '2D Fields' in line[i] ):
            print 'search for 2D fields: success',i
            idx_2d_fields_line = i+1
            fields_list = line[idx_2d_fields_line].split()

    if 'SHGT' in fields_list:
  
        zground_idx = fields_list.index("SHGT") + 1
        print 'zground_idx',zground_idx
        values_list = line[idx_2d_fields_line+2].split()
        z_ground = float(values_list[zground_idx])
        print 'z_ground',z_ground

    else:

        z_ground = 0.0

    if 'TMPS' in fields_list:
  
        TMPS_idx = fields_list.index("TMPS") + 1
        print 'tmps_idx',TMPS_idx
        values_list = line[idx_2d_fields_line+2].split()
        TMPS = float(values_list[TMPS_idx])
        print 'tmps',TMPS

    else:

        TMPS = 0.0


    if 'RH2M' in fields_list:
  
        RH2M_idx = fields_list.index("RH2M") + 1
        print 'RH2M_idx',RH2M_idx
        values_list = line[idx_2d_fields_line+2].split()
        RH2M = float(values_list[RH2M_idx])
        print 'RH2M',RH2M

    else:

        RH2M = 0.0


    if 'PRSS' in fields_list:
  
        prss_idx = fields_list.index("PRSS") + 1
        print 'prss_idx',prss_idx
        values_list = line[idx_2d_fields_line+2].split()
        prss = float(values_list[prss_idx])
        print 'prss',prss

    else:

        prss = 0.0


    for i in range(len(line)):
        if ( '3D Field' in line[i] ):
            print 'search for 3D fields: success',i
            idx_fields_line = i+1
            fields_list = line[idx_fields_line].split()

    a= datetime.datetime.strptime(str(time_start[0:14]), "%y %m %d %H %M")
    b= datetime.datetime.strptime(str(time_input), "%y %m %d %H %M")

    delta = b-a
    time_difference_in_hours = int(delta.total_seconds() / 3600) # hours of difference between time and the start of the wind file

    subprocess.call(profile+ " -d"+run_path+" -f"+meteo_file
                    + " -y"+str(vent_lat)+" -x"+str(vent_lon)+" -o"
                    +str(time_difference_in_hours)+" -p01", shell=True)
    
    with open('atm_profile.txt','w') as file: # output file: wind profile at vent position at time
        with open("profile_01.txt") as fp:
            for i, line in enumerate(fp):
                if i > 13:
                    file.write(str(line[0:64])+'\n')
        fp.close()
    file.close()

             
    #subprocess.call("rm " +run_path+"profile_00.txt", shell=True) 
    #subprocess.call("rm " +run_path+"profile_01.txt", shell=True) 
    

    # read the file 'atm_profile.txt'
    arrays = [np.array(map(float, line.split())) for line in open('atm_profile.txt')]

    # compute the number of values in the first level
    n0 = len(arrays[0])
    b = np.zeros((len(arrays),n0))

    # search for the number of levels with the same number of values
    nfull = 0
    for i in range(len(arrays)):
        if (len(arrays[i]) == n0):
            b[i,:] = arrays[i]
            nfull = nfull+1

    # extract only the lines with the same number of values (n0)
    profile = np.asarray(b[0:nfull,:])

    # compute the missing fields
    a = calc_atm(profile,z_ground,TMPS,RH2M,prss,fields_list)
    a = np.asarray(a)

    if ( z_ground == 0.0 and prss > 0):

        f = interpolate.interp1d(a[2,:],1000*a[0,:])
        z_ground = f(prss)
        print 'z_ground',z_ground

    nrows = a.shape[1]

    f=open('atm.txt','w')
    f.write(str(nrows)+'\n')
    f.close()
    
    f=open('atm.txt','a')
    np.savetxt(f,a.transpose(),fmt='%.4e')    
    f.close()

    f1=open('meteo_ground_elev.txt','w')
    f1.write(str(z_ground)+'\n')
    f1.close()

    return 







