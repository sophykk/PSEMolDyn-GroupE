
#include "FileReader.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "Formulas.h"

#include <iostream>
#include <vector>
#include "outputWriter/VTKWriter.h"
#include <spdlog/spdlog.h>


int main(int argc, char *argsv[]) {

    spdlog::info("Application started");
    spdlog::info("Hello from MolSim for PSE!");

    int N1 = 3, N2 = 3, N3 = 1;
    double h = 1.0, m = 1.0;
    std::array<double, 3> v = {0.0, 0.0, 0.0};
    std::array<double, 3> v2 = {0.0, -10.0, 0.0};
    double x = 0.0, y = 0.0, z = 0.0;
    int type = 1;
    double end_time = 5;
    double delta_t = 0.0002;
    double current_time = 0;

    ParticleGenerator particleGenerator1(40, 8, 1, 1.1225, 1.0, v, 0.0, 0.0, 0.0, type);
    ParticleGenerator particleGenerator2(8, 8, 1, 1.1225, 1.0, v2, 15.0, 15.0, 0.0, type);

    int iteration = 0;

    // Get the ParticleContainer from ParticleGenerator
    ParticleContainer p1 = particleGenerator1.getParticleContainer();
    ParticleContainer p2 = particleGenerator2.getParticleContainer();

    for(auto &p : p2.getParticles()){
        p1.addParticle(p);
    }

    std::vector<Particle> particles = p1.getParticles();

    /*spdlog::info("before createParticlePairs");
    p1.createParticlePairs();
    spdlog::info("after creating particlePairs");*/

    /*for (auto &p:  p1.getParticlePairs()) {
        p.first.setF(Formulas::calculateLJForce(const_cast<std::array<double, 3> &>(p.first.getX()),
                                                 const_cast<std::array<double, 3> &>(p.second.getX()), 1, 5));
        p.second.setF(Formulas::calculateLJForce(const_cast<std::array<double, 3> &>(p.second.getX()),
                                                  const_cast<std::array<double, 3> &>(p.second.getX()), 1, 5));
    }*/

    p1.calculateF(1,5);
    //calculateF(p1);

    Formulas::calculateBM(p1);

    while (current_time < end_time) {
        // calculate new x
        p1.calculateF(1.0 ,5.0);
        //spdlog::info("im while rein");
        p1.calculateV(delta_t);

        p1.calculateX(delta_t);
        // calculate new f

        // calculate new v


        iteration++;
        if (iteration % 10 == 0) {
            p1.plotParticles(iteration);
        }
        spdlog::info("Iteration {} finished.", iteration);

        current_time += delta_t;

    }

    spdlog::info("output written. Terminating...");

     //Logger::getLogger()->info("output written. Terminating...");

    // Get the particles from ParticleContainer
    //std::vector<Particle> particles1 = p1.getParticles();
    //std::vector<Particle> particles2 = p2.getParticles();

    // Output information about each particle
    /*for (auto &particle : particles1) {
        std::cout << particle << std::endl;
    }
    for (auto &particle : particles2) {
        std::cout << particle << std::endl;
    }*/

    //plotParticles(1, p1);
    //plotParticles(1,p2);

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
        p.setF(F_i);
    }
}
