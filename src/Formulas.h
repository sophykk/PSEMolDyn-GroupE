//
// Created by Layla Zadina on 09.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_FORMULAS_H
#define PSEMOLDYN_GROUPE_FORMULAS_H

#include <array>
#include "particle/ParticleContainer.h"

class Formulas {

public:

    /**
   * L2 norm of (xi, xj)
   */
    static double secondNorm(const std::array<double, 3> &arr1);

    /**
   * Lennard-Jones potential
   */
    double calculatePotential(std::array<double, 3> &xi, std::array<double, 3> &xj, double sigma, double eps);

    /**
  * Lennard-Jones force
  */
    static std::array<double, 3> calculateLJForce(std::array<double, 3> &xi, std::array<double, 3> &xj, double sigma, double eps);

    static void calcF(ParticleContainer& cont, double sigma, double eps);

    static void calculateBM(ParticleContainer& pc);
};


#endif //PSEMOLDYN_GROUPE_FORMULAS_H