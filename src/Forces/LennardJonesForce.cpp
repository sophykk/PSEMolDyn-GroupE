//
// Created by markus on 11/25/23.
//

#include "LennardJonesForce.h"
#include "Formulas.h"
#include "utils/ArrayUtils.h"

LennardJonesForce::LennardJonesForce() = default;

void LennardJonesForce::setSigma(const double nSigma) {
    sigma = nSigma;
}

void LennardJonesForce::setEps(const double nEps) {
    eps = nEps;
}

std::array<double, 3> LennardJonesForce::calculateForce(Particle &p1, Particle &p2) {
    if(p1.getType() != p2.getType()){
        setSigma((p1.getSigma() + p2.getSigma()) / 2);
        setEps(std::sqrt(p1.getEps() + p2.getEps()));
    }
    else{
        setSigma(p1.getSigma());
        setEps(p1.getEps());
    }
    auto force = (-24 * eps / pow(Formulas::secondNorm(p1.getX() - p2.getX()), 2.0)) *
                 (pow((sigma / Formulas::secondNorm(p1.getX() - p2.getX())), 6.0) -
                  2 * pow((sigma / Formulas::secondNorm(p1.getX() - p2.getX())), 12.0)) * (p1.getX() - p2.getX());
    force[2] = 0.0;
    return force;
}
