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

#DUMP_ACC_GAS="cdumpcum_gas_$result"

#DUMP_SUM_GAS="cdumpsum_gas_$result"

PDUMP_GAS="pdump_gas_$result"


echo "-------------- exporting plots ---------------"


#python create_maptext.py !check!

echo "'PARTICLES &','### $0 ### &'" >LABELS.CFG

${MDL}/exec/parxplot -i$PDUMP_PART -k1 -z20 -j${MDL}/graphics/arlmap -oparxplot_part.ps

${MDL}/exec/par2asc -i$PDUMP_PART -oPARDUMP_PART.txt 
    
${MDL}/exec/concplot -i$DUMP_PART -j${MDL}/graphics/arlmap -s0 -z20 -d1 -ukg -oconcplot_part.ps

${MDL}/exec/concplot -i$DUMP_ACC_PART -j${MDL}/graphics/arlmap -s0 -t0 -z20 -d1 -ukg -oconcplot_cum_part.ps

if [ $ngas -gt 0 ] 
   then

   echo "'GAS &','### $0 ### &'" >LABELS.CFG

   ${MDL}/exec/parxplot -i$PDUMP_GAS -k1 -z20 -j${MDL}/graphics/arlmap -oparxplot_gas.ps

   ${MDL}/exec/par2asc -i$PDUMP_GAS -oPARDUMP_GAS.txt 
    
   ${MDL}/exec/concplot -i$DUMP_GAS -j${MDL}/graphics/arlmap -s0 -z20 -d1 -ukg -oconcplot_gas.ps

   else

   echo "ngas: "$ngas

fi


rm -f LABELS.CFG

echo "-------------- convert ps to pdf ---------------"

filelist=`(find . -name \*.ps)`

for i in $filelist; do

        ps2pdf $i

        rm $i
done

#python check_deposit.py





