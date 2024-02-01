import numpy as np
import math
import pyproj
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def UmbrellaFitting(runname,data,b,t_sec,T_sec,vent_lat,vent_lon,i):
    
    print("** Umbrella Fitting**")
    r_old_nbl = data[-1,3] # radius at NBL m
    x_old_nbl = data[-1,0]
    y_old_nbl = data[-1,1]
    #print(r_old_nbl,x_old_nbl,y_old_nbl)

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

    T_mean = T_sec + int(t_sec/float(2)) # time at which the fitted values are calculated
    T_sec += t_sec # seconds from the beginning of the emission


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

        alpha=np.arctan2(v_atm_nbl,u_atm_nbl)
        x_new_nbl = d_t * math.cos(alpha)
        y_new_nbl = d_t * math.sin(alpha)


        lat_new = vent_lat + ((y_new_nbl*10**-3)/float(100))
        lon_new = vent_lon + ((x_new_nbl*10**-3)/float(100))
        height_new = b[-1,2]
        emission_area_new = np.pi * r_t**(2)
       
    """
    print("old values")
    print("d_old ",d_old)
    print("r_old ",r_old_nbl)
    print("u_atm_nbl ",u_atm_nbl)
    print("rho_mix_nbl ",rho_mix_nbl)
    print("mfr_nbl ",mfr_nbl)
    print("vfr_nbl ",vfr_nbl)
    print("alpha ",alpha)
    print("new values")
    print("Delta_d ",Delta_d)
    print("d_new_nbl ",d_new_nbl)
    print("x_new_nbl ",x_new_nbl)
    print("y_new_nbl ",y_new_nbl)
    print("r_new_nbl ",r_new_nbl)
    print("new coordinates")
    print("lat_new ",lat_new)
    print("lon_new ",lon_new)
    print("t_steady_fit ",t_steady_fit)
    print("f_t ",f_t)
    print("r_t ",r_t)
    print("d_t ",d_t)
    """

    umbr_fit = runname + '_{0:03}'.format(i)+'.umbrfit'
    with open(umbr_fit,"w") as file:

        file.writelines("Plume properties at NBL \n")
        file.writelines("u_atm_nbl           [m/s]    %f \n"%(u_atm_nbl))
        file.writelines("v_atm_nbl           [m/s]    %f \n"%(v_atm_nbl))
        file.writelines("rho_mix_nbl         [kg/m3]  %f \n"%(rho_mix_nbl))
        file.writelines("mfr_nbl             [kg/s]   %f \n"%(mfr_nbl))
        file.writelines("vfr_nbl             [m3/s]   %f \n"%(vfr_nbl))
        file.writelines("\n")

        file.writelines("Old Values \n")
        file.writelines("r_old_nbl           [m]      %f \n"%(r_old_nbl))
        file.writelines("d_old_nbl           [m]      %f \n"%(d_old))
        file.writelines("\n")

        file.writelines("Fitted Values at steady state\n")
        file.writelines("r_new_nbl           [m]      %f \n"%(r_new_nbl))
        file.writelines("d_new_nbl           [m]      %f \n"%(d_new_nbl))
        file.writelines("t_steady_fit        [s]      %f \n"%(t_steady_fit))
        file.writelines("\n")

        file.writelines("Fitted Values at time t \n")
        file.writelines("t                   [s]      %f \n"%(T_mean)) 
        file.writelines("r_t                 [m]      %f \n"%(r_t))
        file.writelines("d_t                 [m]      %f \n"%(d_t))
        file.writelines("f_t                 -        %f \n"%(f_t))

        file.writelines("lat_new             [deg]    %f \n"%(lat_new))
        file.writelines("lon_new             [deg]    %f \n"%(lon_new))
        file.writelines("height_new          [m]      %f \n"%(height_new))
        file.writelines("emission_area_new   [m2]     %f \n"%(emission_area_new))

    file.close()
    return lat_new,lon_new,height_new,emission_area_new,r_t,t_sec,T_sec
    
def generate_small_sources(r,x_c,y_c,r_c,gaussian_source):
    x_coords = []
    y_coords = []
    source = []

    sigma_c = r_c/3.0
    normal_c = 1/(2.0 * np.pi * sigma_c**2)

    r *= (2/3)
    
    nr = 3*int(np.ceil(r_c/r))
    
    points_x = range(-nr,nr)
    points_y = range(-nr,nr)

    
    for j in points_y:
        for i in points_x:

            if j % 2 == 0:  # Even row
                x_shift = r / 2
            else:
                x_shift = 0

            xp = i * r + x_shift
            yp = j * r * np.sqrt(3) / 2

            if ( (xp**2 + yp**2) <= r_c**2 ):
                
                x_coords.append(x_c+xp)
                y_coords.append(y_c+yp)

                if gaussian_source:

                    dst_c = np.sqrt((xp)**2+(yp)**2)
                    gauss_c = np.exp(-dst_c**2 / (2.0 * sigma_c**2)) * normal_c

                else:

                    gauss_c = 1.0

                source.append(gauss_c)
    print("num of small sources ",len(x_coords))        
    x_coords = np.array(x_coords)            
    y_coords = np.array(y_coords)
    source = np.array(source)

    return x_coords, y_coords, source 
    
# Function to create a map with multiple circles
def create_circles_map(circle_data,lat_new, lon_new,r_c,figname):
    print("create fig")
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    for circle in circle_data:
        center_lat, center_lon, radius_meters, color = circle
        num_points = 100
        circle_lat = [center_lat + radius_meters / 111320.0 * np.cos(2 * np.pi * i / num_points) for i in range(num_points)]
        circle_lon = [center_lon + radius_meters / 111320.0 * np.sin(2 * np.pi * i / num_points) for i in range(num_points)]
        ax.plot(circle_lon, circle_lat, marker='o', markersize=2, color=color, transform=ccrs.PlateCarree())

    # Set the map extent
    span=(r_c*2) * (1/100000)
    top = lat_new + span
    bottom = lat_new - span
    left = lon_new - span
    right = lon_new + span
    ax.set_extent([left, right, bottom, top], crs=ccrs.PlateCarree())

    ax.coastlines(resolution='50m')
    ax.set_title(figname)
    #plt.show()
    fig.savefig(figname+'.png', dpi=1200, bbox_inches='tight')
    plt.close()

    return fig
    
def convert_lat_lon_to_utm(latitude, longitude):
    utm_zone = int((longitude + 180) / 6) + 1  # Determine the UTM zone based on longitude
    utm_crs = f'+proj=utm +zone={utm_zone} +datum=WGS84'
    wgs84_crs = '+proj=longlat +datum=WGS84'

    transformer = pyproj.Transformer.from_crs(wgs84_crs, utm_crs, always_xy=True)
    utm_x, utm_y = transformer.transform(longitude, latitude)

    return utm_x, utm_y, utm_zone

def convert_utm_to_lat_lon(utm_x, utm_y, utm_zone):
    utm_crs = f'+proj=utm +zone={utm_zone} +datum=WGS84'
    wgs84_crs = '+proj=longlat +datum=WGS84'

    transformer = pyproj.Transformer.from_crs(utm_crs, wgs84_crs, always_xy=True)
    longitude, latitude = transformer.transform(utm_x, utm_y)

    return latitude, longitude


