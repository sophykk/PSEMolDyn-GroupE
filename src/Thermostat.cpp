//
// Created by Layla Zadina on 08.12.2023.
//
#include "Containers/LinkedCellContainer2.h"
#include "Thermostat.h"
#include <cmath>
#include <array>
#include "utils/MaxwellBoltzmannDistribution.h"

// Constructor
Thermostat::Thermostat(double Tinit, double Ttarget, double deltaT, int nthermostat, bool useBrownianMotion)
        : Tinit(Tinit), Ttarget(Ttarget == 0 ? Tinit : Ttarget), deltaT(deltaT), nthermostat(nthermostat), useBrownianMotion(useBrownianMotion) {}

/** Apply the thermostat to a LinkedCellContainer2 */
void Thermostat::applyThermostat(LinkedCellContainer2& particleContainer, int currentStep) {
    if (currentStep % nthermostat == 0) {
        /** Systems that have no initial velocities need to be initialized
         * with Brownian Motion to have a non-zero temperature. */
        if (useBrownianMotion) {
            initializeWithBrownianMotion(particleContainer);
        }

        double Tcurrent = calculateCurrentTemperature(particleContainer);
        double scalingFactor = std::sqrt(Ttarget / Tcurrent);
        // Ensuring the temperature change does not exceed deltaT
        scalingFactor = std::min(scalingFactor, std::sqrt((Tcurrent + deltaT) / Tcurrent));
        scaleVelocities(particleContainer, scalingFactor);
    }
}

/** Helper function
 * Calculate the current temperature of the system */
double Thermostat::calculateCurrentTemperature(LinkedCellContainer2& particleContainer) const {
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
    return (2 * Ekin) / (numDimensions * kBoltzmann * particles.size());
}


/** Helper function
 * Initialize particles' velocities with Brownian motion */
void Thermostat::initializeWithBrownianMotion(LinkedCellContainer2& particleContainer) {
        auto& particles = particleContainer.getParticles();
        for (auto& particle : particles) {
            auto bmVelocity = maxwellBoltzmannDistributedVelocity(std::sqrt(Tinit / particle.getM()), 3);
            particle.setV(bmVelocity);
        }

}



/** Helper function
 * Scale the velocities of the particles */
void Thermostat::scaleVelocities(LinkedCellContainer2& particleContainer, double scalingFactor) {
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


