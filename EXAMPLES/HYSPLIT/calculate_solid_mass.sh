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

echo "-------------- check mass ---------------"

${MDL}/exec/con2asc -i$DUMP_PART -t -x -z

mv CON2ASC.OUT CON2ASC.AIR

${MDL}/exec/con2asc -i$DUMP_ACC_PART -t -x -z

mv CON2ASC.OUT CON2ASC.GROUND

python calculate_solid_mass.py

rm CON2ASC.GROUND
rm CON2ASC.AIR


echo
echo "Files created: 1. partial_mass.part (loading of the single particle classes) 2.total_mass.part (total solid particle loading in the domain)" 




