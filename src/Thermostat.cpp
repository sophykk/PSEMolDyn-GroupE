//
// Created by Layla Zadina on 08.12.2023.
//
#include "Thermostat.h"
#include <cmath>
#include "utils/MaxwellBoltzmannDistribution.h"

Thermostat::Thermostat(ParticleContainerBase &pc, double init, int n, bool useBM) : initT(init), nThermostat(n) {
    targetT = initT;
    type = 2;
    deltaT = INFINITY;
    numDimensions = pc.getDimension();
}

Thermostat::Thermostat(ParticleContainerBase &pc, double init, double target, double delta, int n, bool useBM) :
        initT(init), targetT(target), deltaT(delta), nThermostat(n) {
    if (initT < targetT) {
        type = 0;
    } else if (initT > targetT) {
        type = 1;
    } else {
        type = 2;
    }
    numDimensions = pc.getDimension();
}

void Thermostat::initializeWithBrownianMotion(ParticleContainerBase &particleContainer) {
    for (auto &p1: particleContainer.getParticles()) {
        p1.setV(maxwellBoltzmannDistributedVelocity(sqrt(initT / p1.getM()), numDimensions));
    }
}


double Thermostat::calculateNewTemperature(ParticleContainerBase &particleContainer) {
    auto eKin = 0.0;
    for (auto &p1: particleContainer.getParticles()) {
        if (!p1.getIsWall()) {
            eKin += (p1.getM() *
                     (p1.getV()[0] * p1.getV()[0] + p1.getV()[1] * p1.getV()[1] + p1.getV()[2] * p1.getV()[2])) / 2;
        }
    }
    auto divisor = ((numDimensions * particleContainer.getParticles().size()) / 2) * kBoltzmann;
    return eKin / divisor;
}

void Thermostat::applyThermostat(ParticleContainerBase &particleContainer) {
    //if either heating or cooling
    if (type != 2 && ((type == 0 && initT < targetT) || (type == 1 && initT > targetT))) {
        auto tNew = calculateNewTemperature(particleContainer);

        auto scaling = sqrt(tNew / initT);
        if (deltaT != INFINITY) {
            if (std::abs(tNew - initT) > deltaT) {
                if (initT < targetT) {
                    scaling = sqrt((initT + deltaT) / initT);
                    tNew = initT + deltaT;
                } else {
                    scaling = sqrt((initT - deltaT) / initT);
                    tNew = initT - deltaT;
                }
            }
        }
        for (auto &p1: particleContainer.getParticles()) {
            p1.setV({p1.getV()[0] * scaling, p1.getV()[1] * scaling, p1.getV()[2] * scaling});
        }

        initT = tNew;
    }
}

void Thermostat::applyThermostatExtension(ParticleContainerBase &particleContainer) {
    if (type != 2 && ((type == 0 && initT < targetT) || (type == 1 && initT > targetT))) {
        // 1. determine avg velocity
        std::array<double, 3> avgV = {0.0, 0.0, 0.0};
        int n = 0;
        for (auto &p1: particleContainer.getParticles()) {
            if (!p1.getIsWall()) {
                avgV = {avgV[0] + p1.getV()[0], avgV[1] + p1.getV()[1], avgV[2] + p1.getV()[2]};
                n++;
            }
        }
        auto x = 1 / n;
        avgV = {avgV[0] * x, avgV[1] * x, avgV[2] * x};

        // 2. difference btw p's velocity and avg velocity
        for (auto &p1: particleContainer.getParticles()) {
            if (!p1.getIsWall()) {
                p1.setV({p1.getV()[0] - avgV[0], p1.getV()[1] - avgV[1], p1.getV()[2] * avgV[2]});
            }
        }

        // 3. only use difference for calculating temperature
        auto tNew = calculateNewTemperature(particleContainer);

        auto scaling = sqrt(tNew / initT);
        if (deltaT != INFINITY) {
            if (std::abs(tNew - initT) > deltaT) {
                if (initT < targetT) {
                    scaling = sqrt((initT + deltaT) / initT);
                    tNew = initT + deltaT;
                } else {
                    scaling = sqrt((initT - deltaT) / initT);
                    tNew = initT - deltaT;
                }
            }
        }

        // 4. apply scaling to difference and add to avg velocity
        for (auto &p1: particleContainer.getParticles()) {
            if (!p1.getIsWall()) {
                p1.setV({(p1.getV()[0] * scaling) + avgV[0], (p1.getV()[1] * scaling) + avgV[1],
                         (p1.getV()[2] * scaling) + avgV[2]});
            }
        }
        initT = tNew;
    }
}




