#!/bin/sh

docker run --name plumemom_hysplit -d -t -p 8501:8501 -v "$(pwd)":/home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT_DOCKER federicapardini/plumemom_hysplit:v6 

