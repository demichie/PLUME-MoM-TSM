#!/bin/sh

#docker cp /home/ash/Scrivania/Applicazioni/PLUME-MoM-TSM-DOCKER/EXAMPLES/DOCKER/input_file.py plumemom_hysplit:/home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT
docker exec --workdir /home/Codes/PLUME-MoM-TSM/EXAMPLES/PLUMEMoM_DOCKER plumemom /home/Codes/PLUME-MoM-TSM/bin/PLUMEMoM 
