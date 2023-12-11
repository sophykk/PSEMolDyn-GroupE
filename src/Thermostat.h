//
// Created by Layla Zadina on 08.12.2023.
//

#ifndef PSEMOLDYN_GROUPE_THERMOSTAT_H
#define PSEMOLDYN_GROUPE_THERMOSTAT_H

#include <vector>
#include <cmath>
#include "Particle.h"
#include "Containers/LinkedCellContainer2.h"


class Thermostat {
private:
    double Tinit;          // Initial temperature
    double Ttarget;        // Target temperature
    double deltaT;         // Max temperature change per application
    int nthermostat;       // Number of time steps for periodic application
    bool useBrownianMotion;// Flag for initializing with Brownian Motion

public:
    /** Constructor with initial temperature, target temperature, deltaT, and nthermostat */
    Thermostat(double Tinit, double Ttarget = 0, double deltaT = std::numeric_limits<double>::infinity(), int nthermostat = 1, bool useBrownianMotion = false);

    /** Function to apply the thermostat to a LinkedCell particle container */
    void applyThermostat(LinkedCellContainer2& particleContainer, int currentStep);

    /** Helper method to calculate the current temperature of the system */
    double calculateCurrentTemperature(const LinkedCellContainer2& particleContainer) const;

    /** Function to initialize particles' velocities with Brownian motion */
    void initializeWithBrownianMotion(LinkedCellContainer2& particleContainer);

    /** Function to scale the velocities of the particles */
    void scaleVelocities(LinkedCellContainer2& particleContainer, double scalingFactor);
};


#endif //PSEMOLDYN_GROUPE_THERMOSTAT_H
