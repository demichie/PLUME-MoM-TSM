hysplit_dir = "/home/federica/hysplit/trunk"
plumemom_dir = "/home/federica/Scrivania/Codes/PLUME-MoM-TSM"
runname = 'Etna_default'
starttime="18 12 24 11 00" # Year,month,day,hour,minute
endemittime = "18 12 24 12 00"
endruntime = "18 12 24 17 00"
deltat_plumemom = 3600  # seconds

lat = 37.73   # center latitude of the grid
lon = 15.00  # center longitude of the grid
model_top = 32000.0
meteo_file = 'extract_26568.bin'

spacing_lat = 0.1 # degrees between nodes of the sampling grid
spacing_lon = 0.1 # degrees between nodes of the sampling grid
span_lat = 10.00   # the total span of the grid in x direction. For instance, a span of 10 degrees would cover 5 degrees on each side of the center grid location
span_lon = 10.00   # the total span of the grid in y direction. For instance, a span of 10 degrees would cover 5 degrees on each side of the center grid location


vent_lat = 37.73  	# vent latitude
vent_lon = 15.00       # vent longitude
vent_height = 3300    # vent height above sea level (it can be different from ground level of meteo data at vent lat,lon)
vent_velocity = 200.0
log10_mfr = 6.3 # options: log10_mfr or plume_height (plume height above the vent)


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
deltaz_release = 200.0
ncloud = 5


# SETUP parameters
kmsl = 0  	# starting heights default to AGL=0 or MSL=1 *** NOTE: please do not change it (kmsl=0) *** 
ninit = 1  	# particle initialization(0-none; 1-once; 2-add; 3-replace)
ndump = 1  	# dump particles to/from file 0-none or nhrs-output intervall
ncycl = 1 	# pardump output cycle time
numpar = 50000	# number of puffs or particles to released per cycle
maxpar = 1000000 # maximum number of particles carried in simulation
initd = 3 	# initial distribution, particle, puff, or combination.  0 = 3D particle (DEFAULT); 1 = Gh-THv; 2 = THh-THv; 3 = Gh-Pv; 4 = THh-Pv *** NOTE: please use initd=0 or initd=3 *** 
delt = 5	# hysplit integration step (minutes) 
pinpf = ''
kmixd = 0       # flag for boundary layer depth. Default value, see HYSPLIT user guide
kmix0 = 250     # minimum mixing depth. Default value, see HYSPLIT user guide
kzmix = 0       # Vertical Mixing Profile. Default value, see HYSPLIT user guide
kdef = 0        # Horizontal Turbulence. Default value, see HYSPLIT user guide
kbls = 1        # Boundary Layer Stability. Default value, see HYSPLIT user guide
kblt = 2        # Vertical Turbulence. Default value, see HYSPLIT user guide
cmass = 0       # Compute grid concentrations (cmass=0) or grid mass (cmass=1) *** NOTE: please do not change it (cmass=0) *** 

# CONTROL parameters

#SAMPLING INTERVAL
SI_TYPE = 0 # Avg:0 Now:1 Max:2 
SI_HOUR = 0 # hrs 
SI_MINUTE = 5 # min 

#HEIGHT OF EACH CONCENTRATION LEVEL (m-msl)
H_LEVELS = '0 30000'  


# particles parameters
npart = 1
n_sections = 11
phi_min = -4
delta_phi = 1
solid_partial_mass_fraction = 1
phi1 = -1
rho1 = 2500
phi2 = 4
rho2 = 2500
cp_part = 1100
shapefactor = 0.6
mu = 1
sigma = 1

#SAMPLING POINTS !WRONG FOR ETNA!
P01=[-1.431729583, -78.51587596, 10]
P02=[-1.453208011, -78.51823269, 10]
P03=[-1.362751501, -78.59156810, 10]

