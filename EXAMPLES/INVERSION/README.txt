In this folder three inversion test case are included. If the flag INVERSION_FLAG is set to true ("T") in the namelist CONTROL_PARAMETERS, then the namelist INVERSION_PARAMETERS is read from the input file.
In the namelist, the desired plume height is given (HEIGHT_OBJ).
It is then possible to:
1 - assign the initial radius and an interval for the velocity for the search of the optimal value.
2 - assign the initial velocity and an interval for the radius for the search of the optimal value.
3 - assign intervals for both velocity and radius, and the code will search for N_VALUES couples, with radius varying uniformly in a log scale between the minimum and maximum values prescribed.

Example of inversion namelist:


&INVERSION_PARAMETERS
 HEIGHT_OBJ=               15000,
 R_MIN=                      10,
 R_MAX=                      200,
 N_VALUES=          20,
 W_MIN=                      40,
 W_MAX=                      200,
 /

In this case, intervals are prescribed for both radius and velocity, and there is no need to give initial values for the two parameters in the namelist INITIAL_VALUES

&INITIAL_VALUES
 T_MIX0=  1053.0000000000000     ,
 INITIAL_NEUTRAL_DENSITY=F,
 WATER_MASS_FRACTION0=  0.05,
 VENT_HEIGHT=  1500.0000000000000     ,
 DS0=  0.10000000000000000     ,
 N_GAS=          0,
 /

Three different test cases are presented, for the different inversions presented above.

