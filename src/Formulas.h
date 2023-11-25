//
// Created by Layla Zadina on 09.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_FORMULAS_H
#define PSEMOLDYN_GROUPE_FORMULAS_H

#include <array>
#include "Containers/BasicParticleContainer.h"

namespace Formulas {

    /**
   * L2 norm of (xi, xj)
   */
    double secondNorm(const std::array<double, 3> &arr1);

    /**
   * Lennard-Jones potential
   */
    double calculatePotential(std::array<double, 3> &xi, std::array<double, 3> &xj, double mi, double mj);

    /**
  * Lennard-Jones force
  */
    std::array<double, 3> calculateLJForce(const std::array<double, 3> &xi, const std::array<double, 3> &xj,
                                                  double sigma, double eps);

    std::array<double, 3> calculateGravitationForce(const std::array<double, 3> &xi,
                                                           const std::array<double, 3> &xj, double mi, double mj);

    void calculateBM(Particle& p);

    std::array<double, 3> verletXStep(const std::array<double, 3>& x_old, const std::array<double, 3>& v,
                                              const std::array<double, 3>& f, const double m, const double delta_t);

    std::array<double, 3> verletVStep(const std::array<double, 3>& v, const std::array<double, 3>& f_old,
                                             const std::array<double, 3>& f, const double m, const double delta_t);
} // end namespace Formulas


#endif //PSEMOLDYN_GROUPE_FORMULAS_H
