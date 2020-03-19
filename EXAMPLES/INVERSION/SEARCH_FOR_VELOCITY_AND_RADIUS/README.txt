In this folder an inversion test case is included, in which the plume height is prescribed. The user also
assigns intervals for both velocity and radius, and the code will search for N_VALUES couples (velocity,radius), with the radius varying uniformly in a log scale between the minimum and maximum values prescribed.



&INVERSION_PARAMETERS
 HEIGHT_OBJ=               15000,
 R_MIN=                      10,
 R_MAX=                      200,
 N_VALUES=          20,
 W_MIN=                      40,
 W_MAX=                      200,
 /


To run the test case, please copy first (or link) the executable from the bin folder, and then run the model:

> ./PLUMEMoM


