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
    bool useParallelization;
    bool isMembrane = false;
    outputWriter::TXTWriter checkpointingWriter;
    ParticlesFileReader checkpointingReader;
    std::string checkpointingFile = "checkpointing";
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;

    XMLFileReader xmlReader;

    // Read Simulation parameters from the file
    xmlReader.readSimulationParams(argsv[1], end_time, delta_t, modelType, containerType, objectType, plotInterval, checkpointing, useParallelization);

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
        std::array<char, 6> b;
        double g;
        int k;
        double r0;
        double pullUpF;
        xmlReader.readLinkedCellParams(argsv[1], d, c, b, g, isMembrane);
        if(isMembrane){
            xmlReader.readMembraneParams(argsv[1], k, r0, pullUpF);
            particleContainer = std::make_unique<LinkedCellContainer>(*forceModel, d, c, b, g, useParallelization, isMembrane, k, r0, pullUpF);
        } else {
            particleContainer = std::make_unique<LinkedCellContainer>(*forceModel, d, c, b, g, useParallelization, isMembrane);
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
        auto generators = xmlReader.readCuboids(argsv[1], isMembrane);

        // Loop over all generators and let them create particles in the container
        for (auto &gen: generators) {
            gen.generateParticles(*particleContainer);
        }
    }
    // Thermostat setup
    double initialTemperature;
    int thermostatInterval;

    //if not the membrane simulation is run, the thermostat should be initialized
    if(!isMembrane){
        // Read thermostat parameters out of the file
        xmlReader.readThermostatParams(argsv[1], initialTemperature, thermostatInterval);
    }

    //double targetTemperature = initialTemperature;
    //double maxTempChange = std::numeric_limits<double>::infinity();
    bool useBrownianMotion = true;

    Thermostat thermostat(*particleContainer, initialTemperature, thermostatInterval, useBrownianMotion);

    // Velocity is 0 initially for all simulations, so we useBrownianMotion
    if(useBrownianMotion && !isMembrane){
        thermostat.initializeWithBrownianMotion(*particleContainer);
    }

    // Calculate initial forces
    //particleContainer->calculateF();

    int iteration = 0;
    double current_time = 0;

    //Main simulation loop
    //while (current_time < end_time) {
    start = std::chrono::high_resolution_clock::now();
    while(iteration <= 200){

        // calculate new x
        //spdlog::info("this is calcX in MolSim");
        particleContainer->calculateX(delta_t);

        // reset forces
        //spdlog::info("this is resetF in MolSim");
        particleContainer->resetF();

        // calculate new f
        //spdlog::info("this is calcF in MolSim");
        particleContainer->calculateF();

        // apply the thermostat
        if(!isMembrane) {
            if (iteration % thermostatInterval == 0) {
                thermostat.applyThermostatExtension(*particleContainer);
            }
        }

        // calculate new v
        //spdlog::info("this is calcV in MolSim");
        particleContainer->calculateV(delta_t);

        iteration++;
        if (iteration % plotInterval == 1) {
            particleContainer->plotParticles(iteration);
        }
        spdlog::info("Iteration {} finished.", iteration);

        current_time += delta_t;
    }

    spdlog::info("output written. Terminating...");
    // Record end time
    end = std::chrono::high_resolution_clock::now();

    // Calculate duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Print the duration in microseconds
    std::cout << "Runtime: " << duration.count() << " microseconds\n";

    return 0;
}