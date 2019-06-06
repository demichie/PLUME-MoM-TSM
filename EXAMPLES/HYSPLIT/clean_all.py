#!/usr/bin/python


import glob, os
from input_file import *

filelist = glob.glob(runname+'*')
for f in filelist:
    try:
        os.remove(f)
    except (OSError):
        pass

filelist = glob.glob('atm*')
for f in filelist:
    os.remove(f)


filelist = glob.glob('cdump*')
for f in filelist:
    os.remove(f)

filelist = glob.glob('pdump*')
for f in filelist:
    os.remove(f)

filelist = glob.glob('profile*')
for f in filelist:
    os.remove(f)

filelist = glob.glob('PARDUMP*')
for f in filelist:
    os.remove(f)

filelist = glob.glob('MESSAGE*')
for f in filelist:
    os.remove(f)

filelist = glob.glob('*.CFG')
for f in filelist:
        os.remove(f)

filelist = glob.glob('*.txt')
for f in filelist:
        os.remove(f)

filelist = glob.glob('*part')
for f in filelist:
        os.remove(f)


filelist = glob.glob('*gas')
for f in filelist:
        os.remove(f)


filelist = glob.glob('*.pdf')
for f in filelist:
    os.remove(f)


filelist = glob.glob('*.ps')
for f in filelist:
    os.remove(f)

filelist = glob.glob('*.inp')
for f in filelist:
    os.remove(f)


filelist = glob.glob('*.pyc')
for f in filelist:
    os.remove(f)

filelist = glob.glob('*~')
for f in filelist:
    os.remove(f)

filelist = glob.glob('sample_dep*')
for f in filelist:
    os.remove(f)

filelist = [ 'PARDUMP' , 'WARNING' , 'STARTUP' , 'VMSDIST', 'con2stn.txt' , 'plume_model.temp1' , 'plume_model.temp2' , 'check_part.out', 'check_gas.out']

for f in filelist:
    try:
        os.remove(f)
    except OSError:
        pass
                       
     

