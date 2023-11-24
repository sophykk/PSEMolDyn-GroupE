//
// Created by Layla Zadina on 09.11.2023.
//

#include "Formulas.h"
#include "utils/ArrayUtils.h"
#include <cmath>
#include <array>
#include "ParticleContainer.h"


//TODO: refactoring of the MolSim class to extract the formulas?

/**
   * L2 norm (xi,xj)
   */

//template <double, std::size_t N>
double Formulas::secondNorm(const std::array<double, 3> &arr1) {
    double norm = 0.0;
    for (size_t i = 0; i < 3; i++) {
        norm += arr1[i] * arr1[i] ;
    }
    return std::sqrt(norm);
}

/**
   * Lennard-Jones potential
   */

[[maybe_unused]] double Formulas::calculatePotential(std::array<double, 3> &xi, std::array<double, 3> &xj, double sigma, double eps) {
    //U(xi, xj) = 4*eps((sigma/L2_norm(xi,xj))^12 - (sigma/L2_norm(xi,xj))^6)
    auto potential = 4 * eps * (pow((sigma / Formulas::secondNorm(xi - xj)), 12.0) -
                                pow((sigma / Formulas::secondNorm(xi -xj)), 6.0));
    return potential;
}

/**
   * Lennard-Jones force
   */

[[maybe_unused]] std::array<double, 3> Formulas::calculateLJForce(std::array<double, 3> &xi, std::array<double, 3> &xj, double sigma, double eps) {
    //

    auto force = (-24 * eps / pow(Formulas::secondNorm(xi - xj), 2.0)) *
                 (pow((sigma / Formulas::secondNorm(xi - xj)), 6.0) -
                  2 * pow((sigma / Formulas::secondNorm(xi - xj)), 12.0)) * (xi - xj);
    return force;

}

void Formulas::calcF(ParticleContainer cont,double sigma, double eps){
    for (auto &p: cont.getParticles()) {
        //std::array<double, 3> F_i{0., 0., 0.};
        for (auto &p2: cont.getParticles()) {
            // formula: Fij = ((mi * mj) / ||xi −xj||^3) * (xj − xi)
            //std::array<double, 3> F_ij{};
            if (&p != &p2) {
                auto force = (-24 * eps / pow(Formulas::secondNorm(p.getX() - p2.getX()), 2.0)) *
                             (pow((sigma / Formulas::secondNorm(p.getX() - p2.getX())), 6.0) -
                              2 * pow((sigma / Formulas::secondNorm(p.getX() - p2.getX())), 12.0)) * (p.getX() - p2.getX());
                p.setF(force);
                }
            }

        }

    }


