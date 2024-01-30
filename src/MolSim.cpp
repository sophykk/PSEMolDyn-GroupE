#include "inputOutput/inputReader/CuboidFileReader.h"
#include "Containers/BasicParticleContainer.h"
#include "Containers/LinkedCellContainer.h"
#include "Forces/GravitationalForce.h"
#include "Forces/LennardJonesForce.h"
#include "inputOutput/inputReader/XMLFileReader.h"
#include "inputOutput/inputReader/ParticlesFileReader.h"
#include "inputOutput/outputWriter/TXTWriter.h"

#include <spdlog/spdlog.h>
#include <iostream>
#include <chrono>

#include "Thermostat.h"

#include "utils/MaxwellBoltzmannDistribution.h"

bool isWithinCuboid(const std::array<double, 3>& position, const std::array<double, 3>& corner, const std::array<double, 3>& dimensions) {
    return position[0] >= corner[0] && position[0] < corner[0] + dimensions[0] &&
           position[1] >= corner[1] && position[1] < corner[1] + dimensions[1] &&
           position[2] >= corner[2] && position[2] < corner[2] + dimensions[2];
}

int main(int argc, char *argsv[]) {


    double end_time;
    double delta_t;
    std::string modelType;
    std::string containerType;
    std::string objectType;
    int plotInterval;
    spdlog::level::level_enum log_level = spdlog::level::debug;
    bool calcRunTime = false;
    bool checkpointing;
    bool isMembrane = false;
    outputWriter::TXTWriter checkpointingWriter;
    ParticlesFileReader checkpointingReader;
    std::string checkpointingFile = "checkpointing";
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;

    XMLFileReader xmlReader;

    // Read Simulation parameters from the file
    xmlReader.readSimulationParams(argsv[1], end_time, delta_t, modelType, containerType, objectType, plotInterval, checkpointing);

    spdlog::set_level(log_level);

    //std::cout << "XML Parameters " << end_time << delta_t << modelType << containerType << objectType << plotInterval << checkpointing << std::endl;

    spdlog::info("Application started");
    spdlog::info("Hello from MolSim for PSE!");

    // Select to chosen force model
    std::unique_ptr<ForceBase> forceModel;
    if (modelType == "gravitation") {
        forceModel = std::make_unique<GravitationalForce>();
    }
    else if (modelType == "lennardJones") {
        forceModel = std::make_unique<LennardJonesForce>();
    } else {
        spdlog::error("Unknown force model selected: {}", modelType);
        exit(-1);
    }

    // Select to chosen container type
    std::unique_ptr<ParticleContainerBase> particleContainer;
    if (containerType == "basic") {
        particleContainer = std::make_unique<BasicParticleContainer>(*forceModel);
    } else if (containerType == "linkedCells") {
        // Read the parameters for the LinkedCell Container from the file

        //what are those parameters?

        std::vector<double> d;
        double c;
        std::array<char, 6> b;
        double g;
        int k;
        double r0;
        double pullUpF;

        std::cout << "d:";
        for (const auto& value : d) {
            std::cout << value << " ";
        }
        std::cout << ", c:" << c << ", b:";
        for (const char& elem : b) {
            std::cout << elem;
        }
        std::cout << ", g:" << g << ", k:" << k
                  << ", r0:" << r0 << ", pullUpF:" << pullUpF << std::endl;



        xmlReader.readLinkedCellParams(argsv[1], d, c, b, g, isMembrane);
        if(isMembrane){


            xmlReader.readMembraneParams(argsv[1], k, r0, pullUpF);
            particleContainer = std::make_unique<LinkedCellContainer>(*forceModel, d, c, b, g, isMembrane, k, r0, pullUpF);

        } else {

           /** here */
            particleContainer = std::make_unique<LinkedCellContainer>(*forceModel, d, c, b, g, isMembrane);

            std::cout << "Particles size in the beginning: " << particleContainer->getParticles().size() << ", "
                      << std::endl;
            //particleContainer->getParticles().size()
        }
    } else {
        spdlog::error("Unknown container type selected: {}", containerType);
        exit(-1);
    }

    if(objectType == "sphere"){
        // Read Spheres from the file

        auto generators = xmlReader.readSpheres(argsv[1]);

        // Loop over all generators and let them create particles in the container
        for (auto &gen: generators) {
            gen.generateParticles(*particleContainer);

        }
    }
    else{

        // Read Cuboids from the file
        /** here */

        auto generators = xmlReader.readCuboids(argsv[1]);

        std::cout << "Number of generators:" << generators.size() << std::endl;

        // Loop over all generators and let them create particles in the container
        int i = 1;
        for (auto &gen: generators) {
            gen.generateParticles(*particleContainer);
            std::cout << "particleContainer->size() for generator:"<< i << std::endl;
            std::cout << particleContainer->size() << std::endl;
            i++;
        }

        int j = 0;



        for (const auto& particle : particleContainer->getParticles()) {

            /** Printing all particles position, force, velocity  -> terminal_dumb file
             * Output written down in the text file, you can check out the data for 100k particles - good luck! :)
             * DEBUGGING IS FUN
             *
            auto position = particle.getX();
            auto force = particle.getF();
            auto velocity = particle.getV();
            std::cout << "Particle " << j << std::endl;
            // Printing position, force, and velocity
            std::cout << "Particle position: (" << position[0] << ", " << position[1] << ", " << position[2] << ")\n";
            std::cout << "Particle force: (" << force[0] << ", " << force[1] << ", " << force[2] << ")\n";
            std::cout << "Particle velocity: (" << velocity[0] << ", " << velocity[1] << ", " << velocity[2] << ")\n";
            std::cout << "----------------------------------------\n";
            j++; */
        }

        std::array<double, 3> firstCuboidCorner = {0.6, 0.6, 0.6};
        std::array<double, 3> secondCuboidCorner = {0.6, 24.6, 0.6};
        std::array<double, 3> cuboidDimensions = {(50 * 1.2), (20 * 1.2), (50 * 1.2)}; // Width, Height, Depth

        int firstCuboidCount = 0, secondCuboidCount = 0, outsideCount = 0;

        for (const auto& particle : particleContainer->getParticles()) {
            auto position = particle.getX();

            if (isWithinCuboid(position, firstCuboidCorner, cuboidDimensions)) {
                firstCuboidCount++;
            } else if (isWithinCuboid(position, secondCuboidCorner, cuboidDimensions)) {
                secondCuboidCount++;
            } else {
                outsideCount++;
            }
        }

        std::cout << "First Cuboid Particle Count: " << firstCuboidCount << "\n";
        std::cout << "Second Cuboid Particle Count: " << secondCuboidCount << "\n";
        std::cout << "Particles Outside Cuboids: " << outsideCount << "\n";

    }


    // Thermostat setup
    double initialTemperature;
    int thermostatInterval;

    //if not the membrane simulation is ran, the thermostat should be initialized
    if(!isMembrane){
        // Read thermostat parameters out of the file
        xmlReader.readThermostatParams(argsv[1], initialTemperature, thermostatInterval);
        std::cout << "Thermostat Parameters: " << initialTemperature << ", " << thermostatInterval << std::endl;
    }

    double targetTemperature = initialTemperature;
    double maxTempChange = std::numeric_limits<double>::infinity();
    bool useBrownianMotion = true;

    Thermostat thermostat(*particleContainer, initialTemperature, thermostatInterval, useBrownianMotion);

    std::cout << "Particles size after allocating Thermostat: " << particleContainer->getParticles().size() << std::endl;

    // Velocity is 0 initially for all simulations, so we useBrownianMotion
    if(useBrownianMotion && !isMembrane){
        std::cout << "useBrownianMotion && !isMembrane" << std::endl;
        thermostat.initializeWithBrownianMotion(*particleContainer);
    }

    // Calculate initial forces
    particleContainer->plotParticles(0);  // works


    std::cout << "Particles size before calculateF(): " << particleContainer->getParticles().size() << std::endl;
    particleContainer->calculateF();
    std::cout << "Particles size after calculateF111(): " << particleContainer->getParticles().size()  << std::endl;

    int iteration = 0;
    double current_time = 0;

    //Main simulation loop
    while (current_time < end_time) {

        if(checkpointing && current_time == 15){

            // Read Sphere from the file
            auto generators = xmlReader.readSpheres(argsv[1]);

            // Loop over all generators and let them create particles in the container
            for (auto &gen: generators) {
                gen.generateParticles(*particleContainer);
            }

            //write particles to a file
            checkpointingWriter.plotParticles(particleContainer->getParticles(), checkpointingFile);

            // Read the liquid and the drop out from the file
            checkpointingReader.readFile(particleContainer->getParticles(), checkpointingFile);
        }

        if(calcRunTime){
            start = std::chrono::high_resolution_clock::now();
        }

        // calculate new x
        //spdlog::info("this is calcX in MolSim");
        //particleContainer->calculateX(delta_t);

        // reset forces
        //spdlog::info("this is resetF in MolSim");
        particleContainer->resetF(); // resets forces on all particles but saves them as old_f

        // calculate new f
        //spdlog::info("this is calcF in MolSim");
        particleContainer->calculateF();

        if((!checkpointing || current_time < 15) && !isMembrane){
            if (iteration % thermostatInterval == 0) {
                thermostat.applyThermostat(*particleContainer);
            }
        }

        // calculate new v
        //spdlog::info("this is calcV in MolSim");
        particleContainer->calculateV(delta_t);
        particleContainer->calculateX(delta_t);
        // apply the thermostat


        if(calcRunTime){
            // Record end time
            end = std::chrono::high_resolution_clock::now();

            // Calculate duration
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

            // Print the duration in microseconds
            std::cout << "Runtime: " << duration.count() << " microseconds\n";
        }

        iteration++;
        particleContainer->plotParticles(1);

        if (iteration % plotInterval == 0) {
           // particleContainer->plotParticles(iteration);
        }
        spdlog::info("Iteration {} finished.", iteration);

        current_time += delta_t;
    }

    spdlog::info("output written. Terminating...");

    return 0;
}