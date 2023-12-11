//
// Created by Layla Zadina on 09.11.2023.
//
#include "Formulas.h"
#include "utils/ArrayUtils.h"
#include <cmath>
#include <array>
#include "utils/MaxwellBoltzmannDistribution.h"

namespace Formulas {
/**
   * L2 norm (xi,xj)
   */
    double secondNorm(const std::array<double, 3> &arr1) {
        double norm = 0.0;
        for (size_t i = 0; i < 3; i++) {
            norm += arr1[i] * arr1[i];
        }
        return std::sqrt(norm);
    }


    void calculateBM(Particle &p) {
        auto bm = maxwellBoltzmannDistributedVelocity(0.1, 2);
        p.setV(p.getV() + bm);
    }

    std::array<double, 3> verletXStep(const std::array<double, 3> &x_old, const std::array<double, 3> &v,
                                      const std::array<double, 3> &f, const double m, const double delta_t) {
        return x_old + delta_t * v + (delta_t * delta_t) / (2.0 * m) * f;
    }

    std::array<double, 3> verletVStep(const std::array<double, 3> &v, const std::array<double, 3> &f_old,
                                      const std::array<double, 3> &f, const double m, const double delta_t) {
        return v + delta_t / (2 * m) * (f_old + f);
    }

}
