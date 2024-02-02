#!/bin/sh    
InpVar=$1
RUNNAME="${InpVar}"
echo $RUNNAME
docker exec --workdir /home/Codes/PLUME-MoM-TSM/EXAMPLES/PLUMEMoM_DOCKER plumemom python /home/Codes/PLUME-MoM-TSM/UTILS/plot_plume.py -run "$RUNNAME"
