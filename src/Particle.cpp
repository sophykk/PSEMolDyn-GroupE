/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <iostream>
#include "utils/ArrayUtils.h"
#include <spdlog/spdlog.h>

Particle::Particle(int type_arg) {
    type = type_arg;
    spdlog::debug("Particle generated!");
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    gGrav = 0.0;
    sigma = 0.0;
    epsilon = 0.0;
}

Particle::Particle(const Particle &other) {
    x = other.x;
    v = other.v;
    f = other.f;
    old_f = other.old_f;
    m = other.m;
    type = other.type;
    gGrav = other.gGrav;
    sigma = other.sigma;
    epsilon = other.epsilon;
    //spdlog::debug("Particle generated by copy!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, double gGrav_arg,
                   double sigma_arg, double epsilon_arg, int type_arg) : x(x_arg), v(v_arg), m(m_arg), gGrav(gGrav_arg),
                                                                         sigma(sigma_arg), epsilon(epsilon_arg),
                                                                         type(type_arg) {
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    // spdlog::debug("Particle generated!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, std::array<double, 3> f_arg,
                   std::array<double, 3> oldf_arg, double m_arg, double gGrav_arg, double sigma_arg, double epsilon_arg,
                   int type_arg) : x(x_arg), v(v_arg), f(f_arg), old_f(oldf_arg), m(m_arg), gGrav(gGrav_arg),
                                   sigma(sigma_arg), epsilon(epsilon_arg), type(type_arg) {
}

Particle::~Particle() { //spdlog::debug("Particle destructed!"); }
}

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<double, 3> &Particle::getF() const { return f; }

const std::array<double, 3> &Particle::getOldF() const { return old_f; }

double Particle::getM() const { return m; }

double Particle::getGGrav() const { return gGrav; }

double Particle::getSigma() const { return sigma; }

double Particle::getEps() const { return epsilon; }


void Particle::setX(const std::array<double, 3> &x_new) {
    x = x_new;
}

void Particle::setV(const std::array<double, 3> &v_new) {
    v = v_new;
}


void Particle::addF(const std::array<double, 3> &f_add) {
    f = f + f_add;
    f[1] = f[1] + (getM() * gGrav);
}

void Particle::resetF() {
    old_f = f;
    f = {0.0, 0.0, 0.0};
}

int Particle::getType() const { return type; }

std::string Particle::toString() const {
    std::stringstream stream;
    stream << "Particle: X:" << x << " v: " << v << " f: " << f
           << " old_f: " << old_f << " sigma: " << sigma << " epsilon: " << epsilon << " gGrav: " << gGrav << " type: " << type;
    return stream.str();
}

bool Particle::operator==(const Particle &other) const {
    return (x == other.x) and (v == other.v) and (f == other.f) and
           (type == other.type) and (m == other.m) and (old_f == other.old_f) and
            (sigma == other.sigma) and (epsilon == other.epsilon) and (gGrav == other.gGrav);
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
    stream << p.toString();
    return stream;
}
