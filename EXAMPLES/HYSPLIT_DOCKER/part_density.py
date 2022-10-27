import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from extract_wind import write_atm


def calc_density(diam_phi):

    from input_file import phi1, phi2, rho1, rho2, npart, n_sections

    # diam is in millimiters while diam1 and diam2 are in meters

    phi1=np.ones(npart)*phi1
    phi2=np.ones(npart)*phi2

    rho1 = np.ones(npart)*rho1
    rho2 = np.ones(npart)*rho2

    density = np.zeros((npart,n_sections))

    for i in range(npart):

        for j in range(n_sections):

            if ( diam_phi[j] <= phi1[i] ):

                density[i,j] = rho1[i]

            elif ( phi1[i] < diam_phi[j] < phi2[i] ):
 
                density[i,j] = rho1[i] + ( diam_phi[j] - phi1[i] ) / ( phi2[i] - phi1[i] ) * ( rho2[i] - rho1[i] )
       
            elif ( diam_phi[j] >= phi2[i] ):

                density[i,j] = rho2[i]

	#print 'diam',diam[i],diam1[i],diam2[i],diam_phi[i],density[i]

       
    return density
