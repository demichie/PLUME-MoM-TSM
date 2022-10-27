#!/bin/sh

echo "### $0 ###"


#-------------------------------------------------------------
result=$(grep -i 'hysplit_dir' input_file.py | cut -c 15-42)

temp="${result%\"}"
result="${temp#\"}"
temp="${result%\'}"
result="${temp#\'}"

MDL="$result"

result=$(grep -i 'runname' input_file.py | cut -c 11-)

temp="${result%\"}"
result="${temp#\"}"

temp="${result%\'}"
result="${temp#\'}"

DUMP_GAS="cdump_gas_$result"

#----------------------------------------------------------

${MDL}/exec/con2asc -i$DUMP_GAS -t -x -z

mv CON2ASC.OUT CON2ASC.GAS

python3 plot_gas_air.py 

rm CON2ASC.GAS
rm *.gas
