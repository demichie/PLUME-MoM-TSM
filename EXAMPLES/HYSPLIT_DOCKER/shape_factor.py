import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from extract_wind import write_atm


def calc_shapefactor(diam_phi):

    try:

        from input_file import shapefactor, npart, n_sections
     
        shapefactor=np.asarray(shapefactor)
        shapefactor=shapefactor.reshape((npart,1))
        shapefactor = np.ones((npart,n_sections))*shapefactor

    except:

        from input_file import phi1, phi2, shape1, shape2, npart, n_sections

        # diam is in millimiters while diam1 and diam2 are in meters

        phi1=np.ones(npart)*phi1
        phi2=np.ones(npart)*phi2

        shape1 = np.ones(npart)*shape1
        shape2 = np.ones(npart)*shape2

        shapefactor = np.zeros((npart,n_sections))

        for i in range(npart):

            for j in range(n_sections):
  
                if ( diam_phi[j] <= phi1[i] ):

                    shapefactor[i,j] = shape1[i]

                elif ( phi1[i] < diam_phi[j] < phi2[i] ):
 
                    shapefactor[i,j] = shape1[i] + ( diam_phi[j] - phi1[i] ) / ( phi2[i] - phi1[i] ) * ( shape2[i] - shape1[i] )
       
                elif ( diam_phi[j] >= phi2[i] ):

                    shapefactor[i,j] = shape2[i]

	    #print 'diam',diam[i],diam1[i],diam2[i],diam_phi[i],density[i]
   
    return shapefactor
