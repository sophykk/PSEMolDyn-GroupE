#include "inputOutput/inputReader/FileReader.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "Formulas.h"

#include <iostream>
#include <vector>
#include "inputOutput/outputWriter/VTKWriter.h"
#include <spdlog/spdlog.h>


int main(int argc, char *argsv[]) {

    double end_time = 5;
    double delta_t = 0.0002;
    double current_time = 0;


    FileReader cuboidReader;
    auto generators = cuboidReader.readFile(argsv[1]);

    spdlog::info("Application started");
    spdlog::info("Hello from MolSim for PSE!");

    ParticleContainer p1;

    for (auto& gen : generators) {
        gen.generateParticles(p1);
    }

    int iteration = 0;

    std::vector<Particle> particles = p1.getParticles();

    Formulas::calculateBM(p1);
    p1.calculateF(1,5);

    while (current_time < end_time) {

        // reset forces
        p1.resetF();

        // calculate new x
        p1.calculateX(delta_t);

        // calculate new f
        p1.calculateF(1.0 ,5.0);

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


void calculateF(ParticleContainer cont) {

    for (auto &p: cont.getParticles()) {
        std::array<double, 3> F_i{0., 0., 0.};
        for (auto &p2: cont.getParticles()) {
            // formula: Fij = ((mi * mj) / ||xi −xj||^3) * (xj − xi)
            std::array<double, 3> F_ij{};
            if (&p != &p2) {
                auto mul = p.getM() * p2.getM() * (p2.getX() - p.getX());

                for (int i = 0; i < 3; ++i) {
                    F_ij[i] = mul[i] / pow(Formulas::secondNorm((p.getX() - p2.getX())), 3.0);
                    F_i[i] += F_ij[i];
                }
            }
        }
        p.addF(F_i);
    }
}
