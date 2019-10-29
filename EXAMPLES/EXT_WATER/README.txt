This folder contains the input file for a test case with addition of external liquid water at the vent. The test case is taken from Koyaguchi and Woods [1996] (see Fig3).

&WATER_PARAMETERS
 RHO_LW=  1000.D0,
 RHO_ICE= 920.D0,
 ADDED_WATER_TEMP=  273.D0,
 ADDED_WATER_MASS_FRACTION= 0.03,

ADDED_WATER_MASS_FRACTION set the mass fraction of external liquid water with respect to the whole eruptive mixture. 

Mass flow rate of 5.7E6 kgs-1 is computed from vent radius (R0 = 50 m) and initial mixture velocity (W0 = 178 ms-1).

Try to change ADDED_WATER_MASS_FRACTION adjusting W0 to mantain the reference mass flow rate of 5.7E6 kgs-1.

ADDED_WATER_MASS_FRACTION   |   W0
                            |
0                           |   100
0.03                        |   178
0.10                        |   295
0.20                        |   329
0.30                        |   260

To run the test case, please copy firt (or link) the executable from the bin folder, and then run the model:

> ./PLUMEMoM

A Python script to plot the results is provided with the code in the folder UTILS. Copy the script in this folder and then run it:

> python plot_plume.py

The script also create an animation showing how the particle size distribution of both original particles and aggregates (plotted with two different colors) changes with the rise.

****
- Reference
Koyaguchi, T., & Woods, A. W. (1996). On the formation of eruption columns following explosive mixing of magma and surface‐water. Journal of Geophysical Research: Solid Earth, 101(B3), 5561-5574.
****
