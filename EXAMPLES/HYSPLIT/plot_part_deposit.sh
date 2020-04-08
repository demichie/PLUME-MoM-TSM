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

DUMP_ACC_PART="cdumpcum_part_$result"

#----------------------------------------------------------

${MDL}/exec/con2asc -i$DUMP_ACC_PART -t -x -z

mv CON2ASC.OUT CON2ASC.GROUND

python3 plot_part_deposit.py 

rm CON2ASC.GROUND
rm *.gnd
