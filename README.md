# Welcome to PLUME-MoM-TSM

&nbsp;

<img src="https://github.com/demichie/PLUME-MoM-TSM/raw/master/PLUME-MoM-TSM.png?invert_in_darkmode" align=middle width=617.85445pt/>

PLUME-MoM-TSM is a FORTRAN90 code designed to solve the equations for a steady-state integral volcanic plume model, describing the rise in the atmosphere of a mixture of gas and volcanic ash during an eruption. 

The model describes the steady-state dynamics of a plume in a 3-D coordinate system, and the two-size moment (TSM) method is adopted to describe changes in grain-size distribution along the plume associated with particle loss from plume margins and with particle aggregation. For this reason, the new version is named PLUME-MoM-TSM. 

For the first time in a plume model, the full Smoluchowski coagulation equation is solved, allowing to quantify the formation of aggregates during the rise of the plume. In addition, PLUME-MOM-TSM allows to model the phase change of water, which can be either magmatic, added at the vent as liquid from external sources, or incorporated through ingestion of moist atmospheric air. 

Finally, the code includes the possibility to simulate the initial spreading of the umbrella cloud intruding from the volcanic column into the atmosphere. A transient shallow water system of equations models the intrusive gravity current, allowing to compute the upwind spreading.

### Authors and Contributors

Mattia de' Michieli Vitturi (@demichie)

Federica Pardini (@federicapardini)

### Installation and execution

Download the PLUME-MoM-TSM package and create the executable with the following commands from a terminal:

>./configure
>
>make
>
>make install

This will create the executable and copy it in the bin folder. You can test the executable copying it in the EXAMPLES folder and running it.

### Documentation

A wiki page describing the model is available at:

[https://github.com/demichie/PLUME-MoM-TSM/wiki](https://github.com/demichie/PLUME-MoM-TSM/wiki) 

Doxygen generated documentation of the code can be found at:

[http://demichie.github.io/PLUME-MoM-TSM/html/](http://demichie.github.io/PLUME-MoM-TSM/html/) 

### Acknowledgments

The development of PLUME-MoM-TSM has been partially funded by the Italian MIUR project Premiale Ash-RESILIENCE, the MIUR FISR project "Sale Operative Integrate e Reti di Monitoraggio del Futuro", and the European project EUROVOLC (grant number 731070). The authors also thank Alessandro Tadini, Kyle Mohr and Silvia Giansante for testing of the code and for the feedbacks provided.

