This folder contains the input file for a test cases with aggregation. A constant aggregation kernel is used, with value beta=10. 

&AGGREGATION_PARAMETERS
 AGGREGATION_MODEL = "constant",
 PARTICLES_BETA0 = 1.D1 ,

All the classes aggregate into a single class (the last one), which in this case is initialized in the namelist PARTICLES_PARAMETERS with a small mass fraction:

SOLID_PARTIAL_MASS_FRACTION= 0.999 0.001,    

The properties (shape factors, densities) of the aggregates can be different from those of the original particles.


To run the test case, please copy firt (or link) the executable from the bin folder, and then run the model:

> ./PLUMEMoM

A Python script to plot the results is provided with the code in the folder UTILS. Copy the script in this folder and then run it:

> python plot_plume.py

The script also create an animation showing how the particle size distribution of both original particles and aggregates (plotted with two different colors) changes with the rise.

