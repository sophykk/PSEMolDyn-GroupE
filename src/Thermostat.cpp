//
// Created by Layla Zadina on 08.12.2023.
//
#include "Containers/LinkedCellContainer.h"
#include "Thermostat.h"
#include <cmath>
#include <array>
#include "utils/MaxwellBoltzmannDistribution.h"
#include <iostream>

// Constructor
Thermostat::Thermostat(double Tinit, double Ttarget, double deltaT, int nthermostat, bool useBrownianMotion)
        : Tinit(Tinit), Ttarget(Ttarget == 0 ? Tinit : Ttarget), deltaT(deltaT), nthermostat(nthermostat), useBrownianMotion(useBrownianMotion) {}

/** Apply the thermostat to a LinkedCellContainer2 */
void Thermostat::applyThermostat(ParticleContainerBase& particleContainer, int currentStep) {
    std::cout << "Particles size in apply Thermostat: " << particleContainer.size() << ", "
     << std::endl;

        /** Systems that have no initial velocities need to be initialized
         * with Brownian Motion to have a non-zero temperature. */


        double Tcurrent = calculateCurrentTemperature(particleContainer);
        if (std::abs(Tcurrent - Ttarget) < std::numeric_limits<double>::epsilon()) return; // Skip if temperature is already at target

        double newTemp;
        if (Tcurrent > Ttarget) {
            newTemp = std::max(Ttarget, Tcurrent - deltaT);
        } else {
            newTemp = std::min(Ttarget, Tcurrent + deltaT);
        }

        double scalingFactor = std::sqrt(newTemp / Tcurrent);
        scaleVelocities(particleContainer, scalingFactor);
        std::cout << "Tcurrent: " << Tcurrent << std::endl;
}


/** Helper function
 * Calculate the current temperature of the system */
double Thermostat::calculateCurrentTemperature(ParticleContainerBase& particleContainer) const {
    double Ekin = 0.0;
    // apply formula (2) from the worksheet to calculate the Kinetic Energy of the particles
    auto& particles = particleContainer.getParticles();

    for (auto& particle : particles) {
        double vSquared = 0.0;
        for (auto velocity : particle.getV()) {
            vSquared += velocity * velocity;
        }

        Ekin += particle.getM() * vSquared / 2.0;
    }
    //apply formula (1) from the worksheet to calculate T
    std::cout << "Ekin: " << Ekin << ", "
              << "numDimensions: " << numDimensions << ", "
              << "kBoltzmann: " << kBoltzmann << ", "
              << "Particles size: " << particles.size() << ", "
              << std::endl;
    return (2 * Ekin) / (numDimensions * kBoltzmann * particles.size());
}


/** Helper function
 * Initialize particles' velocities with Brownian motion */
void Thermostat::initializeWithBrownianMotion(ParticleContainerBase& particleContainer) {
        auto& particles = particleContainer.getParticles();
        for (auto& particle : particles) {
            auto bmVelocity = maxwellBoltzmannDistributedVelocity(std::sqrt(Tinit / particle.getM()), numDimensions);
            particle.setV(bmVelocity);
        }

}



/** Helper function
 * Scale the velocities of the particles */
void Thermostat::scaleVelocities(ParticleContainerBase& particleContainer, double scalingFactor) {
    auto& particles = particleContainer.getParticles();
    for (auto& particle : particles) {
        const auto& currentVelocity = particle.getV();
        std::array<double, 3> scaledVelocity;

        for (size_t i = 0; i < currentVelocity.size(); ++i) {
            scaledVelocity[i] = currentVelocity[i] * scalingFactor;
        }

        particle.setV(scaledVelocity);
    }
}


