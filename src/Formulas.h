//
// Created by Layla Zadina on 09.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_FORMULAS_H
#define PSEMOLDYN_GROUPE_FORMULAS_H

#include <array>
#include "ParticleContainer.h"

class Formulas {

public:

    /**
   * L2 norm of (xi, xj)
   */
    static double secondNorm(const std::array<double, 3> &arr1);

    /**
   * Lennard-Jones potential
   */
    [[maybe_unused]] double calculatePotential(std::array<double, 3> &xi, std::array<double, 3> &xj, double sigma, double eps);

    /**
  * Lennard-Jones force
  */
    static void calculateLJForce(std::vector<std::pair<Particle*, Particle*>>& pairVector, double sigma, double eps);

    static void calcF(ParticleContainer& cont, double sigma, double eps);

    static void calculateBM(ParticleContainer& pc);
};


#endif //PSEMOLDYN_GROUPE_FORMULAS_H
