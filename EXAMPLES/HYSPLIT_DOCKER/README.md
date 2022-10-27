## PLUME-MoM-TSM and HYSPLIT with Docker

- To run a PLUME-MoM-TSM/HYSPLIT simulation with Docker check that Docker is installed on your machine. To do this run `docker info` on a terminal screen. If the command is not found, please install Docker.

- Execute `./run_docker.sh` to run the container *plumemom_hysplit* from the image *federicapardini/plumemom_hysplit:v1* and check that the container is running with `docker ps`.

- Edit *input_file.py* and run the simulation with `./run_plumemom_hysplit_docker.sh`.

- Execute the post-processing routines:
  - `./create_plots_docker.sh` (HYSPLIT plotting routines for solid particles and gas concentrations)

  - `./plot_part_deposit_docker.sh` to plot the deposit 

  - `./plot_part_air_docker.sh` to plot atmospheric particle concentration

  - `./plot_gas_air_docker.sh` to plot atmospheric concentration of volcanic gasses different than water vapour

  - `./extract_samples_docker.sh` to extract particles mass loading at specific lat-lon locations

  - `./calculate_part_mass_docker.sh` to compute the mass of solid particles in the domain 

  - `./calculate_gas_mass_docker.sh` to compute the mass of volcanic gases in the domain

- To clean the working directory execute `python clean_all.py`

- To stop and remove the container *plumemom_hysplit* execute `./end_docker.sh` 

The scripts provided in this folder have been written by F. Pardini<sup>(1)</sup> and M. de' Michieli Vitturi<sup>(1)</sup>

(1) Istituto Nazionale di Geofisica e Vulcanologia - Sezione di Pisa</br>    federica.pardini@ingv.it mattia.demichielivitturi@ingv.it













