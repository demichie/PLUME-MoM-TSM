from time import gmtime, strftime
import datetime
from input_file import *

jobstart = strftime("%Y-%m-%d %H:%M:%S", gmtime())

# create a second template with the parameters constant in time
f = open('MAPTEXT.template','r')
filedata = f.read()
f.close()

filedata = filedata.replace("{jobstart}", str(jobstart) )

filedata = filedata.replace("{runname}", str(runname) )

filedata = filedata.replace("{lat}", str(lat) )
filedata = filedata.replace("{lon}", str(lon) )
filedata = filedata.replace("{vent_height}", str(vent_height) )

time_format = "%y %m %d %H %M"




runtime = datetime.datetime.strptime(endruntime,time_format) - datetime.datetime.strptime(starttime,time_format)
emittime = datetime.datetime.strptime(endemittime,time_format) - datetime.datetime.strptime(starttime,time_format)

time_format = "%y %m %d %H %M"


released_mass = emittime.total_seconds() * 10**log10_mfr * ( 1.0 - gas_mass_fraction )


filedata = filedata.replace("{released_mass}", str(released_mass) )

filedata = filedata.replace("{starttime}", str(datetime.datetime.strptime(starttime,time_format)) )

filedata = filedata.replace("{endemittime}", str(datetime.datetime.strptime(endemittime,time_format)) )

filedata = filedata.replace("{meteo_file}", str(meteo_file) )


f = open('MAPTEXT.CFG','w')
f.write(filedata)

f.close()

