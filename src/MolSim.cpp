#include "inputOutput/inputReader/CuboidFileReader.h"
#include "Containers/BasicParticleContainer.h"
#include "Forces/GravitationalForce.h"
#include "Forces/LennardJonesForce.h"

#include <spdlog/spdlog.h>

int main(int argc, char *argsv[]) {

    double end_time = 5;
    double delta_t = 0.0002;
    std::string modelType = "lennardJones";
    std::string containerType = "basic";
    spdlog::level::level_enum log_level = spdlog::level::debug;
    int plotInterval = 10;

    ///////////////////////////////////////////////////
    // Read xml file here and set parameters/options !!
    spdlog::set_level(log_level);

    spdlog::info("Application started");
    spdlog::info("Hello from MolSim for PSE!");

    CuboidFileReader cuboidReader;
    auto generators = cuboidReader.readFile(argsv[1]);

    // Select to chosen force model
    std::unique_ptr<ForceBase> forceModel;
    if (modelType == "gravitation") {
        forceModel = std::make_unique<GravitationalForce>();
    }
    else if (modelType == "lennardJones") {
        forceModel = std::make_unique<LennardJonesForce>(1.0, 5.0);
    } else {
        spdlog::error("Unknown force model selected: {}", modelType);
        exit(-1);
    }

    // Select to chosen container type
    std::unique_ptr<ParticleContainerBase> particleContainer;
    if (containerType == "basic") {
        particleContainer = std::make_unique<BasicParticleContainer>(*forceModel);
    }
    else if (containerType == "linkedCells") {
        spdlog::error("Linked Cells container is not yet implemented!");
        exit(-1);
    } else {
        spdlog::error("Unknown container type selected: {}", containerType);
        exit(-1);
    }

    // Loop over all generators and let them create particles in the container
    for (auto& gen : generators) {
        gen.generateParticles(*particleContainer);
    }

    // Calculate initial forces
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
