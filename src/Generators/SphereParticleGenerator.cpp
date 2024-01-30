//
// Created by markus on 11/30/23.
//
#include "SphereParticleGenerator.h"
#include "Particle.h"
#include "Formulas.h"
#include <cmath>

/**
 * @param h distance between two particle(mesh width of the grid)
 * @param m mass of one particle
 * @param v initial velocity v of the particles (3 Components)
 * @param x coordinate of the lower left front-side corner (x,y,z)
 * @param r radius of the sphere
 */

SphereParticleGenerator::SphereParticleGenerator(std::array<double, 3> x1, double h1, double m1, std::array<double, 3> v1,
                                                 double r1, double sigma1, double epsilon1, double gGrav1, int type1, bool isWall1)
        : x(x1), h(h1), m(m1), v(v1), r(r1), sigma(sigma1), epsilon(epsilon1), gGrav(gGrav1), type(type1), isWall(isWall1) {
}

/**
 * @brief create a grid of particles with mass m, initial velocity v within a sphere of radius r and position of the center x
*/
void SphereParticleGenerator::generateParticles(ParticleContainerBase& particleContainer) {

    // Define the number of Particles in the grid
    int numParticlesX = static_cast<int>(std::ceil(2.0 * r / h));
    int numParticlesY = static_cast<int>(std::ceil(2.0 * r / h));

    // Iterate through each particle
    for (int i = 0; i < numParticlesX; ++i) {
        for (int j = 0; j < numParticlesY; ++j) {
            // Calculate the position of the current particle within the grid
            double x1 = x[0] - r + i * h;
            double y1 = x[1] - r + j * h;

            // Check if the particle is inside the circular boundary
            if (std::hypot(x1 - x[0], y1 - x[1]) <= r) {
                Particle particle({x1, y1, x[2]}, v, m, gGrav, sigma, epsilon, type, isWall);
                Formulas::addMB(particle);
                particleContainer.addParticle(particle);
            }
        }
    }
}
