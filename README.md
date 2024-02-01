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

Files: -w5t1.xml (Task 1)
       -w5t3.xml (Task 3)
       -w5t4.xml (Task 4)
 
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

- task2 can be seen on the parallelization branch
- added the "-pg" flag to the compiler flags
- after the programm finishes executing, the most time consuming parts of the code can be seen in the analysis.txt file, by executing the following command in the terminal: gprof MolSim gmon.out > analysis.txt
- extended the xsd schema, so the parallelization strategy and number of threads can be read from the xml file
- if no strategy is given in the xml file, the code will be executed without using one of the parallelization strategies
- implemented a runtime logic in the main function, which outputs the total time in microseconds of executing the programm, used for creating the plot of the speed up
- strategy 1:
   - allowed a shared acces of the grid to all threads in the initGrid() function of the LinkedCellContainer 
   - each thread receives its own copies of the loop indexes
- strategy 2:
   - collapsed 3 loops into a single one for parallelization in the calculateF() function of the LinkedCellContainer
- unfortunately, after trying out different parallelization strategies, we couldnt find a suitable parallelization strategy to speedup our simulation
- for both strategies, the runtime was slower than executing it without a parallelization strategy
    
  
## Task 3 ##

- extended the Rayleigh-Taylor Instability to 3D
- Cuboid Generator generates the cuboids correctly in the beginning of the simualtion.
- we had some problems with disappearing particles. Particle container became smaller at every iteration. Found a problematic place in calcF function in LinkedCell Container. Some particles were having negative X position after initGrid() function.
- Solved problem with disappearing particles, but simulation was not running correctly. Particles do not move after some iteration.
- Rendering takes too much time, could not run the whole simulation.

## Task 4 ##

- extended the thermostat class with a new method, that ...
- added a new variable to the cuboid generator, which allocates to each particle if it is a wall or not
- if a particle is a wall, no force is applied to it and its position shouldn't change
- created an csv file exporter, which creates 50 bins of equal size and adds the particles, based on the x-axis, to the corresponding files
- then for every 10.000 iterations, a csv file is created, which outputs on each line the bin index, the density and the average velocity of the bin 


