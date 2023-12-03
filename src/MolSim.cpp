#include "inputOutput/inputReader/CuboidFileReader.h"
#include "Containers/BasicParticleContainer.h"
#include "Containers/LinkedCellContainer2.h"
#include "Forces/GravitationalForce.h"
#include "Forces/LennardJonesForce.h"
#include "inputOutput/inputReader/XMLFileReader.h"

#include <spdlog/spdlog.h>
#include <iostream>

//
int main(int argc, char *argsv[]) {


    double end_time;
    double delta_t;
    std::string modelType;
    std::string containerType;
    int plotInterval;
    double epsilon;
    double sigma;
    spdlog::level::level_enum log_level = spdlog::level::debug;

    XMLFileReader xmlReader;

    // Read Simulation parameters from the file
    xmlReader.readSimulationParams(argsv[1], end_time, delta_t, modelType, containerType, plotInterval);

    spdlog::set_level(log_level);

    spdlog::info("Application started");
    spdlog::info("Hello from MolSim for PSE!");

    // Read cuboids from the file
    auto generators = xmlReader.readCuboids(argsv[1]);

    // Select to chosen force model
    std::unique_ptr<ForceBase> forceModel;
    if (modelType == "gravitation") {
        forceModel = std::make_unique<GravitationalForce>();
    }
    else if (modelType == "lennardJones") {
        // Read sigma and epsilon from teh file
        xmlReader.readLennardJonesForceParams(argsv[1], sigma,epsilon);
        forceModel = std::make_unique<LennardJonesForce>(sigma, epsilon);
    } else {
        spdlog::error("Unknown force model selected: {}", modelType);
        exit(-1);
    }

    // Select to chosen container type
    std::unique_ptr<ParticleContainerBase> particleContainer;
    if (containerType == "basic") {
        particleContainer = std::make_unique<BasicParticleContainer>(*forceModel);
    } else if (containerType == "linkedCells") {
        std::vector<double> d{180.0, 90.0, 1.0};
        double c = 3.0;
        //boundary 'r' or 'o'
        particleContainer = std::make_unique<LinkedCellContainer2>(*forceModel, d, c, 'r');

        // spdlog::error("Linked Cells container is not yet implemented!");
        // exit(-1);
    } else {
        spdlog::error("Unknown container type selected: {}", containerType);
        exit(-1);
    }

    // Loop over all generators and let them create particles in the container
    for (auto &gen: generators) {
        gen.generateParticles(*particleContainer);
    }

    // Calculate initial forces

    //spdlog::debug("This is before calculate F");

    particleContainer->calculateF();

    // Main simulation loop
    int iteration = 0;
    double current_time = 0;
    while (current_time < end_time) {

        // calculate new x
        particleContainer->calculateX(delta_t);

        // reset forces
        particleContainer->resetF();

        // calculate new f
        particleContainer->calculateF();

        // calculate new v
        particleContainer->calculateV(delta_t);

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
