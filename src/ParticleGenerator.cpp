//
// Created by Layla Zadina on 09.11.2023.
//
#include "ParticleGenerator.h"
#include "Particle.h"

ParticleGenerator::ParticleGenerator(std::array<int, 3> N1, double h1, double m1, std::array<double, 3> v1,
                                     std::array<double, 3> x1, int type1)
        : N(N1), h(h1), m(m1), v(v1), x(x1), type(type1) {
/** Number of particles per dimension: N1 × N2 × N3
// Distance h of the particles (mesh width of the grid)
// Mass m of one particle
// Initial velocity v of the particles (3 Components)
 The coordinate of the lower left front-side corner (x,y,z)
 */
}

void ParticleGenerator::generateParticles(ParticleContainer& particleContainer) {
    /**
    * create a N1xN2xN3 grid of particles with mass m, initial velocity v, left corner coordinate (x,y,z)
   */
    for (int i = 0; i < N[0]; ++i) {
        for (int j = 0; j < N[1]; ++j) {
            for (int k = 0; k < N[2]; ++k) {
                Particle particle({x[0] + i * h, x[1] + j * h, x[2] + k * h}, v, m, type);
                Formulas::calculateBM(particle);
                particleContainer.addParticle(particle);
            }
        }
    }
}




