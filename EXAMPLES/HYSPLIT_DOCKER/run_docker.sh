#!/bin/sh

xhost +local:root

docker run --name plumemom_hysplit -d -t --net=host -e DISPLAY -v /tmp/.X11-unix  -v "$(pwd)":/home/Codes/PLUME-MoM-TSM-master/EXAMPLES/HYSPLIT_DOCKER federicapardini/plumemom_hysplit:v1 

