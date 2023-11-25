#include "inputOutput/inputReader/FileReader.h"
#include "utils/ArrayUtils.h"
#include "Containers/ParticleContainer.h"
#include "ParticleGenerator.h"
#include "Formulas.h"
#include "Forces/GravitationalForce.h"
#include "Forces/LennardJonesForce.h"

#include <vector>
#include "inputOutput/outputWriter/VTKWriter.h"
#include <spdlog/spdlog.h>

int main(int argc, char *argsv[]) {

    double end_time = 5;
    double delta_t = 0.0002;
    double current_time = 0;
    std::string modelType = "lennardJones";

    std::unique_ptr<ForceBase> forceModel;

    if (modelType == "gravitation") {
        forceModel = std::make_unique<GravitationalForce>();
    }
    else if (modelType == "lennardJones") {
        forceModel = std::make_unique<LennardJonesForce>(1.0, 5.0);
    } else {
        spdlog::error("Unknown force model selected: {}", modelType);
    }

    FileReader cuboidReader;
    auto generators = cuboidReader.readFile(argsv[1]);

    spdlog::info("Application started");
    spdlog::info("Hello from MolSim for PSE!");

    ParticleContainer p1(*forceModel);

    for (auto& gen : generators) {
        gen.generateParticles(p1);
    }

    int iteration = 0;

    p1.calculateF();

    while (current_time < end_time) {

        // calculate new x
        p1.calculateX(delta_t);

        // reset forces
        p1.resetF();

        // calculate new f
        p1.calculateF();

        // calculate new v
        p1.calculateV(delta_t);

        iteration++;
        if (iteration % 10 == 0) {
            p1.plotParticles(iteration);
        }
        spdlog::info("Iteration {} finished.", iteration);

        current_time += delta_t;
    }

    spdlog::info("output written. Terminating...");
    return 0;
}
