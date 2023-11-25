//
// Created by markus on 11/25/23.
//

#include "GravitationalForce.h"
#include "Formulas.h"

std::array<double, 3> GravitationalForce::calculateForce(Particle& p1, Particle& p2) {
    auto mul = p1.getM() * p2.getM() / pow(Formulas::secondNorm((p1.getX() - p2.getX())), 3.0);
    auto force = mul * (p1.getX() - p2.getX());

    return force;
}