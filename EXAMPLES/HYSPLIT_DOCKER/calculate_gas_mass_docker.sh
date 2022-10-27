#!/bin/sh

docker exec --workdir /home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT_DOCKER plumemom_hysplit ./calculate_gas_mass.sh

#docker cp plumemom_hysplit_v3:/home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT/total_mass.gas $PWD 
#docker cp plumemom_hysplit_v3:/home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT/partial_mass.gas $PWD 

