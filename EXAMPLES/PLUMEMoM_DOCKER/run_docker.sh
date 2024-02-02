#!/bin/sh
docker run --name plumemom -d -t -p 8501:8501 -v "$(pwd)":/home/Codes/PLUME-MoM-TSM/EXAMPLES/PLUMEMoM_DOCKER federicapardini/plumemom:2024v3
