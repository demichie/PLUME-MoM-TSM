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

ngas=$(grep -i 'ngas' input_file.py | cut -c 8-11)

DUMP_PART="cdump_part_$result"

DUMP_ACC_PART="cdumpcum_part_$result"

DUMP_SUM_PART="cdumpsum_part_$result"

PDUMP_PART="pdump_part_$result"

DUMP_GAS="cdump_gas_$result"

DUMP_ACC_GAS="cdumpcum_gas_$result"

DUMP_SUM_GAS="cdumpsum_gas_$result"

PDUMP_GAS="pdump_gas_$result"

#----------------------------------------------------------

python3 run_plumemom.py || exit

python3 create_hysplit_emittimes_control_part.py || exit

python3 create_hysplit_setup_ascdata.py || exit

echo "-------------- particles dispersion simulation ---------------"

echo $MDL
${MDL}/exec/hycs_std part  


if [ $ngas -gt 0 ] 
  then
  echo "-------------- gas dispersion simulation ---------------"
  python3 create_hysplit_emittimes_control_gas.py || exit
  ${MDL}/exec/hycs_std gas  

  else
  echo "-------------- skipped gas simulation - ngas: $ngas ---------------"

fi

rm *.txt

echo "-------------- start postprocessing ---------------"

${MDL}/exec/concacc -i$DUMP_PART -o$DUMP_ACC_PART

${MDL}/exec/concsum -i$DUMP_ACC_PART -o$DUMP_SUM_PART

if [ $ngas -gt 0 ]
  then
  ${MDL}/exec/concacc -i$DUMP_GAS -o$DUMP_ACC_GAS
  ${MDL}/exec/concsum -i$DUMP_ACC_GAS -o$DUMP_SUM_GAS
fi

echo "-------------- end postprocessing ---------------"


