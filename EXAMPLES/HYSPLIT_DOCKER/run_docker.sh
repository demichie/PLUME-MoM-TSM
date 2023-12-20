#!/bin/sh

docker run --name plumemom_hysplitDec2023 -d -t -p 8501:8501 -v "$(pwd)":/home/Codes/PLUME-MoM-TSM/EXAMPLES/HYSPLIT_DOCKER federicapardini/plumemom_hysplit:4 

