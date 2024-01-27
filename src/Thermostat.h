//
// Created by Layla Zadina on 08.12.2023.
//

#ifndef PSEMOLDYN_GROUPE_THERMOSTAT_H
#define PSEMOLDYN_GROUPE_THERMOSTAT_H

#include <vector>
#include <cmath>
#include "Particle.h"
#include "Containers/LinkedCellContainer.h"
#include <limits>


class Thermostat {
private:
    static constexpr int numDimensions = 2;
    static constexpr double kBoltzmann = 1.0;
    double Tinit;          // Initial temperature
    double Ttarget;        // Target temperature
    double Tcurrent;
    double deltaT;         // Max temperature change per application
    int nthermostat;       // Number of time steps for periodic application

    /** useBrownianMotion = True if initial Velocities == 0, False otherwise
     * initialize Thermostat with this flag*/
    bool useBrownianMotion;// Flag for initializing with Brownian Motion

public:
    /** Constructor with initial temperature, target temperature, deltaT, and nthermostat */
    Thermostat(double Tinit, double Ttarget = 0, double deltaT = std::numeric_limits<double>::infinity(), int nthermostat = 1, bool useBrownianMotion = false);

    /** Function to apply the thermostat to a LinkedCell particle container */
    void applyThermostat(ParticleContainerBase& particleContainer, int currentStep);

    /** Helper method to calculate the current temperature of the system */
    double calculateCurrentTemperature(ParticleContainerBase& particleContainer) const;

    /** Function to initialize particles' velocities with Brownian motion */
    void initializeWithBrownianMotion(ParticleContainerBase& particleContainer);

    /** Function to scale the velocities of the particles */
    void scaleVelocities(ParticleContainerBase& particleContainer, double scalingFactor);

    double getCurrentT() const;

    double getTargetT() const;
};


#endif //PSEMOLDYN_GROUPE_THERMOSTAT_H
