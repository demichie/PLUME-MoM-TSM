The scripts provided in this folder have been written by F.Pardini<sup>(1)</sup> and M.de' Michieli Vitturi<sup>(1)</sup>

To run a PLUME-MoM-TSM/HYSPLIT simulation with Docker:

1) be sure that Docker is installed on your machine. To do this run `$ docker info` on a terminal screen. If the command is not found, please install Docker.

2) Execute `./run_docker.sh` to run the container *plumemom_hysplit* from the image *federicapardini/plumemom_hysplit:v1* and check that the container is running with `docker ps`.

3) Edit *input_file.py* 

4) Execute `./run_plumemom_hysplit_docker.sh`   - To run the PLUMEMoM/HYSPLIT simulation with Docker

5) Execute `./create_plots_docker.sh`           - HYSPLIT plotting routines for solid particles and gas concentrations

6) Execute `./plot_part_deposit_docker.sh`       - To plot the deposit 

7) Execute `./plot_part_air_docker.sh`           - To plot atmospheric particle concentration

8) Execute `./plot_gas_air_docker.sh`            - To plot atmospheric concentration of volcanic gasses different than water vapour

7) Execute `./extract_samples_docker.sh`        - To extract particles mass loading at specific lat-lon locations (see input_file.py)

10) Execute `./calculate_part_mass_docker.sh`    - To compute the mass of solid particles in the domain 

11) Execute `./calculate_gas_mass_docker.sh`     - To compute the mass of volcanic gases in the domain

12) Execute `./python clean_all.py`                      - To clean the working directory

From point 5 to point 12 the execution order can be varied by the user

To stop and remove the container *plumemom_hysplit* execute `./end_docker.sh` 

(1) Istituto Nazionale di Geofisica e Vulcanologia
    Sezione di Pisa
    federica.pardini@ingv.it
    mattia.demichielivitturi@ingv.it

 


