hysplit_dir = "/home/Codes/hysplit.v5.0.1"         
plumemom_dir = "/home/Codes/PLUME-MoM-TSM"    #please do not modify the paths of PLUME-MoM and Hysplit
runname = 'Etna_test'
starttime   = "18 12 24 12 00" # Year,month,day,hour,minute
endemittime = "18 12 24 13 00"
endruntime  = "18 12 24 17 00"
deltat_plumemom = 3600  # seconds

lat = 37.73   # center latitude of the grid
lon = 15.00  # center longitude of the grid
model_top = 32000.0

meteo_file_dir="./"
meteo_file = ['extract_26568.bin'] # add meteo files separated by a comma, the first file is used to exract the atmospheric profile for plumemom
spacing_met_grid = 50 #km

spacing_lat = 0.05 # degrees between nodes of the sampling grid
spacing_lon = 0.05 # degrees between nodes of the sampling grid
span_lat = 5.00   # the total span of the grid in x direction. For instance, a span of 10 degrees would cover 5 degrees on each side of the center grid location
span_lon = 5.00   # the total span of the grid in y direction. For instance, a span of 10 degrees would cover 5 degrees on each side of the center grid location

vent_lat = 37.75       # vent latitude
vent_lon = 15.00       # vent longitude
vent_height = 3300   # vent height above sea level (it can be different from ground level of meteo data at vent lat,lon)
tmix0 = 1273 # mixture temperature in Kelvin
vent_velocity = 200.0

#options for mass flow rate or plume height (nbl - meters above the vent), select only one:
#log10_mfr = 7 
mfr = 1E5
#plume_height = 7000 #plume height at nbl. inversion: search for radius given vent_velocity
#compute mfr from vent_velocity and vent_radius: set vent_velocity and add vent_radius 
# N.B: write [X1,X2,X3,...,XN] to set a variable X (mass flow rate or plume height) with time (i.e. N plume-mom runs each with a different Xi).

#umbrella cloud parameters
umbrella_flag = "Fit" # Model: Umbrella expansion computed by the Shallow water model - Fit: Umbrella expansion computed by the fitting formula - False: No Umbrella
small_sources_flag = True # Valid only for umbrella Fit. 
gaussian_source = True
plot_fig = True
t_end = 3600 # duration of the umbrella simulation (seconds)
dt_output = 600 # saving output at dt_output time steps (seconds)
c_d = 0.1 # drag coefficient
steady_flag = "T" # stop the umbrella simulation when the steady state is reached

# volcanic gas parameters
ngas = 0   # in addition to H2O
rvolcgas = [189] # CO2 and SO2 R constant [J/kgK]
cpvolcgas = [844]
volcgas_mol_wt = [0.044]
volcgas_mass_fraction = [0.05]

#initial volcanic water mass fraction
water_mass_fraction0 = 0.03

#flag for water condensation - freezing - addition of external liquid water at the vent
water_flag = 'F'

#external water parametes
rho_lw =  1000.0
rho_ice =  920.0
added_water_temp =  273.0
added_water_mass_fraction =  0.1

# hysplit parameters
deltaz_release = 500.0
ncloud = 1

# HYSPLIT SETUP parameters
kmsl = 0  	# starting heights default to AGL=0 or MSL=1 *** NOTE: please do not change it (kmsl=0) *** 
ninit = 1  	# particle initialization(0-none; 1-once; 2-add; 3-replace)
ndump = 1  	# dump particles to/from file 0-none or nhrs-output intervall
ncycl = 1 	# pardump output cycle time
numpar = -100	# number of puffs or particles to released per cycle
maxpar = 1000000 # maximum number of particles carried in simulation
initd = 3 	# initial distribution, particle, puff, or combination.  0 = 3D particle (DEFAULT); 1 = Gh-THv; 2 = THh-THv; 3 = Gh-Pv; 4 = THh-Pv *** NOTE: please use initd=0 or initd=3 *** 
delt = 10	# hysplit integration step (minutes) 
pinpf = ''
kmixd = 0       # flag for boundary layer depth. Default value, see HYSPLIT user guide
kmix0 = 250     # minimum mixing depth. Default value, see HYSPLIT user guide
kzmix = 0       # Vertical Mixing Profile. Default value, see HYSPLIT user guide
kdef = 0        # Horizontal Turbulence. Default value, see HYSPLIT user guide
kbls = 1        # Boundary Layer Stability. Default value, see HYSPLIT user guide
kblt = 2        # Vertical Turbulence. Default value, see HYSPLIT user guide
cmass = 0       # Compute grid concentrations (cmass=0) or grid mass (cmass=1) *** NOTE: please do not change it (cmass=0) *** 

#SAMPLING INTERVAL
SI_TYPE = 0 # Avg:0 Now:1 Max:2 
SI_HOUR = 0 # hrs 
SI_MINUTE = 10 # min 

#HEIGHT OF EACH CONCENTRATION LEVEL (m-agl)
H_LEVELS = '0 30000'  

# particle parameters
npart = 1
n_sections = 10
phi_min = -4
delta_phi = 1
solid_partial_mass_fraction = 1
rcoef_flag = "T" # True or false

# particle distribution can be LOGNORMAL or BIN
distribution = "LOGNORMAL"

# parameters for LOGNORMAL distribution
mu = 0.0
sigma = 1.2

# parameters for BIN distribution
bin_partial_mass_fraction = [1.9E-03 , 1.1E-02, 4.3E-02 ,1.1E-01 ,2.1E-01, 2.5E-01, 2.1E-01 ,1.1E-01, 4.3E-02, 1.1E-02, 1.9E-03, 2.2E-04, 1.7E-05, 9.3E-07, 3.4E-08, 8.2E-10, 1.4E-11, 1.5E-13]

# parameters for particle densitysf_range = 30
phi1 = -1
rho1 = 2500
phi2 = 4
rho2 = 2500

cp_part = 1100 
sf_range = 30

shapefactor = 0.9
settling_formulation = "ganser" #Possible options: ganser and stokes

#SAMPLING POINTS: LAT - LON - SAMPLED_DEPOSIT (kgm-2)
P01=[36.00, 18.00, 1]
P02=[36.01, 18.01, 1]
P03=[35.99, 18.00, 1]
P04=[37.73, 15.00, 1]


