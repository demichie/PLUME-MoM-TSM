import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from part_density import calc_density
from input_file import *

with open('SETUP.part','w') as setup:
    setup.write('&SETUP  \n')
    setup.write('initd='+str(initd)+'\n')
    setup.write('kmsl='+str(kmsl)+'\n')
    setup.write('ninit='+str(ninit)+'\n')
    setup.write('ndump='+str(ndump)+'\n')
    setup.write('ncycl='+str(ncycl)+'\n')
    setup.write('delt='+str(delt)+'\n')
    setup.write("efile = 'EMITTIMES.part', \n")
    setup.write('numpar='+str(numpar)+'\n')
    setup.write('maxpar='+str(maxpar)+'\n')
    setup.write("pinpf ='"+pinpf+"'\n")
    setup.write("poutf ='pdump_part_"+runname+"'\n")
    setup.write("kmixd="+str(kmixd)+"\n")
    setup.write("kmix0="+str(kmix0)+"\n")
    setup.write("kzmix="+str(kzmix)+"\n")
    setup.write("kdef="+str(kdef)+"\n")
    setup.write("kbls="+str(kbls)+"\n")
    setup.write("kblt="+str(kblt)+"\n")
    setup.write("cmass="+str(cmass)+"\n")
    setup.write('/ \n')
setup.close()


with open('SETUP.gas','w') as setup:
    setup.write('&SETUP  \n')
    setup.write('initd='+str(initd)+'\n')
    setup.write('kmsl='+str(kmsl)+'\n')
    setup.write('ninit='+str(ninit)+'\n')
    setup.write('ndump='+str(ndump)+'\n')
    setup.write('ncycl='+str(ncycl)+'\n')
    setup.write('delt='+str(delt)+'\n')
    setup.write("efile = 'EMITTIMES.gas', \n")
    setup.write('numpar='+str(numpar)+'\n')
    setup.write('maxpar='+str(maxpar)+'\n')
    setup.write("pinpf ='"+pinpf+"'\n")
    setup.write("poutf ='pdump_gas_"+runname+"'\n")
    setup.write("kmixd="+str(kmixd)+"\n")
    setup.write("kmix0="+str(kmix0)+"\n")
    setup.write("kzmix="+str(kzmix)+"\n")
    setup.write("kdef="+str(kdef)+"\n")
    setup.write("kbls="+str(kbls)+"\n")
    setup.write("kblt="+str(kblt)+"\n")
    setup.write('/ \n')
setup.close()


with open('ASCDATA.CFG','w') as ascdata:
    ascdata.write('-90.0  -180.0 \n')	
    ascdata.write('1.0     1.0 \n')
    ascdata.write('180     360 \n')
    ascdata.write('2 \n')
    ascdata.write('0.2 \n')
    ascdata.write("'"+hysplit_dir+"/bdyfiles/'\n")
ascdata.close()


