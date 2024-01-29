//
// Created by Layla Zadina on 08.12.2023.
//
#include "Thermostat.h"
#include <cmath>
#include "utils/MaxwellBoltzmannDistribution.h"
#include <iostream>

Thermostat::Thermostat(ParticleContainerBase &pc, double init, int n, bool useBM) : initT(init), nThermostat(n) {
    targetT = initT;
    deltaT = INFINITY;
    numDimensions = pc.getDimension();

    std::cout << "Hello from Thermostat!" << std::endl;
    std::cout << "numDimensions " << numDimensions << std::endl;

    if (useBM) {
        initializeWithBrownianMotion(pc);
    }
}

Thermostat::Thermostat(ParticleContainerBase &pc, double init, double target, double delta, int n, bool useBM) :
        initT(init), targetT(target), deltaT(delta), nThermostat(n) {
    numDimensions = pc.getDimension();
    std::cout << "Hello from Thermostat!" << std::endl;
    std::cout << "numDimensions " << numDimensions << std::endl;
    if (useBM) {
        initializeWithBrownianMotion(pc);
    }
}

void Thermostat::initializeWithBrownianMotion(ParticleContainerBase &particleContainer) {
    for (auto &p1: particleContainer.getParticles()) {
        p1.setV(maxwellBoltzmannDistributedVelocity(sqrt(initT / p1.getM()), numDimensions));
    }
}


double Thermostat::calculateNewTemperature(ParticleContainerBase &particleContainer) {
    auto eKin = 0.0;
    for (auto &p1: particleContainer.getParticles()) {
        eKin += (p1.getM() *
                 (p1.getV()[0] * p1.getV()[0] + p1.getV()[1] * p1.getV()[1] + p1.getV()[2] * p1.getV()[2])) / 2;
    }
    auto divisor = ((numDimensions * particleContainer.getParticles().size()) / 2) * kBoltzmann;
    return eKin / divisor;
}

void Thermostat::applyThermostat(ParticleContainerBase &particleContainer) {
    if (initT != targetT) {
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

/*
/** Apply the thermostat to a LinkedCellContainer
void Thermostat::applyThermostat(ParticleContainerBase &particleContainer, int currentStep) {

    /** Systems that have no initial velocities need to be initialized
     * with Brownian Motion to have a non-zero temperature.


    double currentT = calculateCurrentTemperature(particleContainer);
    double scalingFactor = std::sqrt(targetT / currentT);
    // Ensuring the temperature change does not exceed deltaT
    scalingFactor = std::min(scalingFactor, std::sqrt((currentT + deltaT) / currentT));
    scaleVelocities(particleContainer, scalingFactor);
    //std::cout << "Tcurrent: " << Tcurrent << std::endl;
}
*/




