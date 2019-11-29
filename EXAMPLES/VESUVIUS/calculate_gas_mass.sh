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

ngas=$(grep -i 'ngas' input_file.py | cut -c 8-11)

DUMP_GAS="cdump_gas_$result"


#----------------------------------------------------------

if [ $ngas -gt 0 ]
  then
  echo "-------------- check mass ---------------"

  ${MDL}/exec/con2asc -i$DUMP_GAS -t -x -z

  mv CON2ASC.OUT CON2ASC.AIR

  python calculate_gas_mass.py

  rm CON2ASC.AIR

  echo
  echo "Files created 1. partial_mass.gas (gas mass for the single gas species) 2. total_mass.gas (total gas mass in the domain)" 

else

   echo "-------------- ngas: $ngas ---------------"

fi
