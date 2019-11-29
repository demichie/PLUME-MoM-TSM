#!/bin/sh

echo "### $0 ###"


#-------------------------------------------------------------
result=$(grep -i 'hysplit_dir' input_file.py | cut -c 15-)

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

DUMP_PART="cdump_part_$result"

#----------------------------------------------------------

${MDL}/exec/con2asc -i$DUMP_PART -t -x -z

mv CON2ASC.OUT CON2ASC.AIR

python plot_concentration_FL.py 

rm CON2ASC.AIR
rm *.air
