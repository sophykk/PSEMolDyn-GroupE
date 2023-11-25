//
// Created by Layla Zadina on 09.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_PARTICLEGENERATOR_H
#define PSEMOLDYN_GROUPE_PARTICLEGENERATOR_H

#include <array>
#include "Containers/ParticleContainer.h"

class ParticleGenerator {

private:
    /** Number of particles per dimension: N1 × N2 × N3
    // Distance h of the particles (mesh width of the grid)
    // Mass m of one particle
    // Initial velocity v of the particles (3 Components)
     The coordinate of the lower left front-side corner (x,y,z)
     */

    std::array<int, 3> N;
    double h;
    double m;
    std::array<double, 3> v;
    std::array<double, 3> x;
    int type;

public:
    /**
     * create a N1xN2xN3 grid of particles with mass m, initial velocity v, left corner coordinate (x,y,z)
    */

    ParticleGenerator(std::array<int, 3> N1, double h1, double m1, std::array<double, 3> v1,
                      std::array<double, 3> x1, int type1);

    void generateParticles(ParticleContainer& particleContainer);
};


#endif //PSEMOLDYN_GROUPE_PARTICLEGENERATOR_H
