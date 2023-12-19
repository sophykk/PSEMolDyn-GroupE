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

//
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
    outputWriter::TXTWriter checkpointingWriter;
    ParticlesFileReader checkpointingReader;
    std::string checkpointingFile = "checkpointing";
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;

    XMLFileReader xmlReader;

    // Read Simulation parameters from the file
    xmlReader.readSimulationParams(argsv[1], end_time, delta_t, modelType, containerType, objectType, plotInterval, checkpointing);

    spdlog::set_level(log_level);

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
        std::vector<double> d;
        double c;
        std::array<char, 4> b;
        double g;
        xmlReader.readLinkedCellParams(argsv[1], d, c, b, g);

        particleContainer = std::make_unique<LinkedCellContainer>(*forceModel, d, c, b, g);
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
        auto generators = xmlReader.readCuboids(argsv[1]);

        // Loop over all generators and let them create particles in the container
        for (auto &gen: generators) {
            gen.generateParticles(*particleContainer);
        }
    }

    // Thermostat setup
    double initialTemperature;
    int thermostatInterval;

    // Read thermostat parameters out of the file
    xmlReader.readThermostatParams(argsv[1], initialTemperature, thermostatInterval);

    double targetTemperature = initialTemperature;
    double maxTempChange = std::numeric_limits<double>::infinity();
    bool useBrownianMotion = true;

    Thermostat thermostat(initialTemperature, targetTemperature, maxTempChange, thermostatInterval, useBrownianMotion);

    // Velocity is 0 initially for all simulations, so we useBrownianMotion
    if(useBrownianMotion){
        thermostat.initializeWithBrownianMotion(*particleContainer);
    }

    // Calculate initial forces
    particleContainer->calculateF();

    int iteration = 0;
    double current_time = 15;

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
        particleContainer->calculateX(delta_t);

        // reset forces
        particleContainer->resetF();

        // calculate new f
        particleContainer->calculateF();

        // apply the thermostat
        if(!checkpointing || current_time < 15){
            if (iteration % thermostatInterval == 0) {
                thermostat.applyThermostat(*particleContainer, iteration);
            }
        }

        // calculate new v
        particleContainer->calculateV(delta_t);

        if(calcRunTime){
            // Record end time
            end = std::chrono::high_resolution_clock::now();

            // Calculate duration
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

            // Print the duration in microseconds
            std::cout << "Runtime: " << duration.count() << " microseconds\n";
        }

        iteration++;
        if (iteration % plotInterval == 0) {
            particleContainer->plotParticles(iteration);
        }
        spdlog::info("Iteration {} finished.", iteration);

        current_time += delta_t;
    }

    spdlog::info("output written. Terminating...");

    return 0;
}