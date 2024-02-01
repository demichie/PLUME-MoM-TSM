#!/bin/sh

docker exec -i --workdir /home/Codes/PLUME-MoM-TSM/EXAMPLES/HYSPLIT_DOCKER plumemom_hysplit2024 ./plot_part_air.sh

#docker cp plumemom_hysplit_v3:/home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT/total_mass.part $PWD 
#docker cp plumemom_hysplit_v3:/home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT/partial_mass.part $PWD 

