#!/bin/sh

#docker cp /home/ash/Scrivania/Applicazioni/PLUME-MoM-TSM-DOCKER/EXAMPLES/DOCKER/input_file.py plumemom_hysplit:/home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT
docker exec --workdir /home/Codes/PLUME-MoM-TSM/EXAMPLES/HYSPLIT_DOCKER plumemom_hysplit2024 ./run_plumemom_hysplit.sh
