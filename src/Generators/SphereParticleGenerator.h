//
// Created by markus on 11/30/23.
//


#ifndef PSEMOLDYN_GROUPE_SPHEREPARTICLEGENERATOR_H
#define PSEMOLDYN_GROUPE_SPHEREPARTICLEGENERATOR_H

#include <array>
#include "Containers/BasicParticleContainer.h"

class SphereParticleGenerator {

private:
    /** Coordinates of the center: x1 × x2 × x3
    // Distance h of the particles (mesh width of the grid)
    // Mass m of one particle
    // Initial velocity v of the particles (3 Components)
    // Radius r of the Sphere
     */

    std::array<double, 3> x;
    double h;
    double m;
    std::array<double, 3> v;
    double r;
    double sigma;
    double epsilon;
    double gGrav;
    int type;

public:
    /**
     * create a N1xN2xN3 grid of particles with mass m, initial velocity v, left corner coordinate (x,y,z)
    */

    SphereParticleGenerator(std::array<double, 3> x1, double h1, double m1, std::array<double, 3> v1, double r1,
                            double sigma, double epsilon, double gGrav, int type1);

    void generateParticles(ParticleContainerBase &particleContainer);
};

#endif //PSEMOLDYN_GROUPE_SPHEREPARTICLEGENERATOR_H

