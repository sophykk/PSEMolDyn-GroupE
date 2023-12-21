MolSim
===
# Group E #
Members:
* Layla
* Markus
* Sophy

# Code #
* Link:     https://github.com/sophykk/PSEMolDyn-GroupE.git
* Branch:   master
* Compiler: g++ 11.4.0
* Commit ID: last commit on the master branch
* Build and execute instructions:
 ```
 mkdir build
 cd build
 ccmake ..
 make
 ./MolSim ../input/<file_name>

Files: -w4t2test.xml (Task 2, the test simulation)
       -w4t2.xml (Task 2, the big simulation)
       -w4t3.xml (Task 3)
 
 Doxygen documentation: 
 make doc_doxygen
 To open documentation in the browser:
 doxys_documentation/html/index.html

UnitTests:
ctest
 
 then for the visualisation go into ParaView and open the beforehand generated MD_vtk folder/group, hit apply and use the glyph filter with the according configurations on it
```
# Report #
## Task 1 ##

- implemented Thermostat for controlling the temperature in the system
- current temperature of the system is calculated with the Kinetic Energy of the particles
- velocity scaling factor is calculated for the particles to achieve the target temperature
- systems that have no initial velocities are initialized with Brownian Motion
- delta T is applied, so the temperature change at each application of Thermostat does not exceed the chosen threshold
- MolSim is adapted accordinly, Thermostat is applied periodically at the given interval. Thermostat is applied before Velocity calculation, since it influences the velocity scaling
- XML input is adapted to read the Thermostat parameters

## Task 2 ##

- implemented a new particle container, which realizes the linked cell algo
- LinkedCellContainer contains a grid with cells made out of a special vector structure
- particleList is checked and based on the position the particles is assigned a cell
- generated different sizes of cuboids and calculated the runtime per iteration using the LinkedCellContainer and the BasicParticleContainer
- measured the time for 1000, 2000, 4000 and 8000 particles comparing performance for Basic and Linked-Cell Containers:

<img src="Basic.png" alt="Time measurements for Basic Container" title="Basic" width="200"/> <img src="Linked Cell.png" alt="Time measurements for Linked Cell Container" title="Linked Cell" width="200"/>

### Visualization Plot for the Time Measurements: ###
<img src="./measurements_pic.png" alt="Corresponding Graph" title="Comparison graph" style="width: 70%;">


  
## Task 3 ##

- implemented a new output writer TXTWriter, which writes on the first line of a txt file the number of the particles and then starting with the new line, the values for every particle on its own line; the values of the particles are separated by a space and written in this order: position, velocity, force, old force, mass, sigma, epsilon, gravitational acceleration and type
- implemented a new file reader ParticlesFileReader, which reads out of a file, given as a parameter, the number of the particles and then reads every line one by one and initializes a particle with the corresponding values; the particle is then added to the vector given as a parameter
- updated the main function to execute the simulation until timestamp 15; after that, the drop is generated, the values are saved in a file, then the simulation continues with all particles read from the saved file; for the part after the drop is added to the simulation, the thermostat is switched off
- conducted the experiment of a falling drop by using periodic boundaries on the left and right sides and refleting boundaries for the up and down sides 
