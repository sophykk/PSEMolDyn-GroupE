//
// Created by markus on 11/25/23.
//

#include "LennardJonesForce.h"
#include "Formulas.h"
#include "utils/ArrayUtils.h"

LennardJonesForce::LennardJonesForce(double sigma, double eps) : sigma(sigma), eps(eps) {}

std::array<double, 3> LennardJonesForce::calculateForce(Particle &p1, Particle &p2) {
    auto force = (-24 * eps / pow(Formulas::secondNorm(p1.getX() - p2.getX()), 2.0)) *
                 (pow((sigma / Formulas::secondNorm(p1.getX() - p2.getX())), 6.0) -
                  2 * pow((sigma / Formulas::secondNorm(p1.getX() - p2.getX())), 12.0)) * (p1.getX() - p2.getX());
    return force;
}
