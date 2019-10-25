In this folder an inversion test case is included, in which the plume height, the initial velocity and an interval for the initial radius are prescribed in the INVERSION_PARAMETERS namelist.

&INVERSION_PARAMETERS
 HEIGHT_OBJ=               15000,
 W_MIN=                      40,
 W_MAX=                      200,
 /


To run the test case, please copy firt (or link) the executable from the bin folder, and then run the model:

> ./PLUMEMoM


