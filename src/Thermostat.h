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
    int numDimensions;
    static constexpr double kBoltzmann = 1.0;
    double initT;          // Initial temperature
    double targetT;        // Target temperature
    double deltaT;         // Max temperature change per application
    int nThermostat;       // Number of time steps for periodic application

    /** useBrownianMotion = True if initial Velocities == 0, False otherwise
     * initialize Thermostat with this flag*/
    //bool useBrownianMotion;// Flag for initializing with Brownian Motion

public:

    Thermostat(ParticleContainerBase &pc, double init, double target, double delta, int n, bool useBM);

    Thermostat(ParticleContainerBase &pc, double init, int n, bool useBM);

    /*
    /** Constructor with initial temperature, target temperature, deltaT, and nthermostat
    Thermostat(double initT, double targetT = 0, double deltaT = std::numeric_limits<double>::infinity(),
               int nthermostat = 1, bool useBrownianMotion = false);

    /** Function to apply the thermostat to a LinkedCell particle container*/
    void applyThermostat(ParticleContainerBase &particleContainer);

    /** Helper method to calculate the current temperature of the system */
    double calculateNewTemperature(ParticleContainerBase &particleContainer);

    /** Function to initialize particles' velocities with Brownian motion */
    void initializeWithBrownianMotion(ParticleContainerBase &particleContainer);

};


#endif //PSEMOLDYN_GROUPE_THERMOSTAT_H
