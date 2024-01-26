//
// Created by Layla Zadina on 09.11.2023.
//
#include "CuboidParticleGenerator.h"
#include "Particle.h"
#include "Formulas.h"
#include <iostream>

/**
 * @param N1 Number of particles per dimension: N1 × N2 × N3
 * @param h distance between two particle(mesh width of the grid)
 * @param m mass of one particle
 * @param v initial velocity v of the particles (3 Components)
 * @param x coordinate of the lower left front-side corner (x,y,z)
 */

CuboidParticleGenerator::CuboidParticleGenerator(std::array<int, 3> N1, double h1, double m1, std::array<double, 3> v1,
                                                 std::array<double, 3> x1, double sigma1, double epsilon1,
                                                 double gGrav1, int type1)
        : N(N1), h(h1), m(m1), v(v1), x(x1), sigma(sigma1), epsilon(epsilon1), gGrav(gGrav1), type(type1) {
}

/**
 * @brief create a N1xN2xN3 grid of particles with mass m, initial velocity v, left corner coordinate (x,y,z)
 * */
void CuboidParticleGenerator::generateParticles(ParticleContainerBase &particleContainer) {

    for (int xIndex = 0; xIndex < N[0]; ++xIndex) {
        for (int yIndex = 0; yIndex < N[1]; ++yIndex) {
            for (int zIndex = 0; zIndex < N[2]; ++zIndex) {
                Particle particle({x[0] + xIndex * h, x[1] + yIndex * h, x[2] + zIndex * h}, v, m, gGrav, sigma,
                                  epsilon, type, xIndex, yIndex);
                Formulas::addMB(particle);
                particleContainer.addParticle(particle);
            }
        }
    }
}




