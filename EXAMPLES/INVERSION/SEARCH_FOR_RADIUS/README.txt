In this folder an inversion test case is included, in which the plume height, the initial radius and an interval for the velocity are prescribed in the INVERSION_PARAMETERS namelist.

&INVERSION_PARAMETERS
 HEIGHT_OBJ=               15000,
 R_MIN=                      10,
 R_MAX=                      200,
 /


To run the test case, please copy firt (or link) the executable from the bin folder, and then run the model:

> ./PLUMEMoM


