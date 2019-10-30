This folder contains the input file for a test case with addition of external liquid water at the vent. The test case is taken from Koyaguchi and Woods [1996] (see Fig3).

&WATER_PARAMETERS
 RHO_LW=  1000.D0,
 RHO_ICE= 920.D0,
 ADDED_WATER_TEMP=  273.D0,
 ADDED_WATER_MASS_FRACTION= 0.03,

ADDED_WATER_MASS_FRACTION set the mass fraction of external liquid water with respect to the whole eruptive mixture. 

Mass flow rate of 5.7E6 kgs-1 (MFR0=5.7E6) and vent radius (R0 = 50 m) are fixed.

In order to reproduce the results of Koyaguchi and Woods, these values of external water can be testet.

ADDED_WATER_MASS_FRACTION
                         
0                         
0.03                       
0.10                      
0.20                      
0.30                     

To run the test case, please copy first (or link) the executable from the bin folder, and then run the model:

> ./PLUMEMoM

A Python script to plot the results is provided with the code in the folder UTILS. Copy the script in this folder and then run it:

> python plot_plume.py

The script also create an animation showing how the particle size distribution changes with the rise.

****
- Reference
Koyaguchi, T., & Woods, A. W. (1996). On the formation of eruption columns following explosive mixing of magma and surface‐water. Journal of Geophysical Research: Solid Earth, 101(B3), 5561-5574.
****
