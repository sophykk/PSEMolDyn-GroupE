//
// Created by Layla Zadina on 08.12.2023.
//
#include "Containers/LinkedCellContainer2.h"
#include "Thermostat.h"
#include <cmath>
#include <algorithm>

// Constructor
Thermostat::Thermostat(double Tinit, double Ttarget, double deltaT, int nthermostat, bool useBrownianMotion)
        : Tinit(Tinit), Ttarget(Ttarget == 0 ? Tinit : Ttarget), deltaT(deltaT), nthermostat(nthermostat), useBrownianMotion(useBrownianMotion) {}

/** Apply the thermostat to a LinkedCellContainer2 */
void Thermostat::applyThermostat(LinkedCellContainer2& particleContainer, int currentStep) {
    if (currentStep % nthermostat == 0) {
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

/** Calculate the current temperature of the system */
double Thermostat::calculateCurrentTemperature(const LinkedCellContainer2& particleContainer) const {
    double Ekin = 0.0;
    // apply formula (2) from the worksheet to calculate the Kinetic Energy of the particles
    const auto& particles = particleContainer.getParticles();
    for (const auto& particle : particles) {
        double vSquared = std::inner_product(particle.getV().begin(), particle.getV().end(), particle.getV().begin(), 0.0);
        Ekin += particle.getM() * vSquared / 2.0;
    }
    //apply formula (1) from the worksheet to calculate T 
    return (2 * Ekin) / (3 * particles.size()); // Assuming 3 dimensions
}

// Initialize particles' velocities with Brownian motion
void Thermostat::initializeWithBrownianMotion(LinkedCellContainer2& particleContainer) {
    auto& particles = particleContainer.getParticles();
    for (auto& particle : particles) {
        double fi = std::sqrt(Tinit / particle.getM());
        // Apply Maxwell-Boltzmann distribution scaling here
    }
}

// Scale the velocities of the particles
void Thermostat::scaleVelocities(LinkedCellContainer2& particleContainer, double scalingFactor) {
    auto& particles = particleContainer.getParticles();
    for (auto& particle : particles) {
        auto& v = particle.getV();
        std::transform(v.begin(), v.end(), v.begin(), [scalingFactor](double velocity) { return velocity * scalingFactor; });
    }
}

