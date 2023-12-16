#include "Containers/BasicParticleContainer.h"
#include "Containers/LinkedCellContainer2.h"
#include "Forces/GravitationalForce.h"
#include "Forces/LennardJonesForce.h"
#include "Thermostat.h"
#include <spdlog/spdlog.h>
#include <iostream>
#include <chrono>
#include <memory>
#include <vector>

int main(int argc, char *argv[]) {
    // Hard-coded simulation parameters
    double end_time = 100.0;   // Duration of the simulation
    double delta_t = 0.1;      // Time step
    int plotInterval = 10;     // Interval for plotting

    // Hard-coded force model parameters (Lennard-Jones potential)
    double epsilon = 1.0;  // Depth of the potential well
    double sigma = 1.0;    // Finite distance at which inter-particle potential is zero
    auto forceModel = std::make_unique<LennardJonesForce>(sigma, epsilon);

    // Hard-coded parameters for LinkedCellContainer2
    std::vector<double> domainSize = {10.0, 10.0, 10.0}; // Domain size
    double cutoffRadius = 1.0;                           // Cutoff radius
    char boundaryCondition = 'o';

    // Setup LinkedCellContainer2
    auto particleContainer = std::make_unique<LinkedCellContainer2>(*forceModel, domainSize, cutoffRadius, boundaryCondition);

    // Thermostat setup
    double initialTemperature = 300; // Initial temperature
    double targetTemperature = 300;  // Target temperature
    double maxTempChange = 10;       // Maximum change in temperature per step
    int thermostatInterval = 10;     // Apply thermostat every 10 iterations
    bool initializeWithBrownianMotion = true; // Initialize with Brownian motion

    Thermostat thermostat(initialTemperature, targetTemperature, maxTempChange, thermostatInterval, initializeWithBrownianMotion);

    // Main simulation loop
    spdlog::info("Starting simulation loop");
    int iteration = 0;
    double current_time = 0;
    while (current_time < end_time) {
        // Calculate new positions
        particleContainer->calculateX(delta_t);

        // Reset forces
        particleContainer->resetF();

        // Calculate new forces
        particleContainer->calculateF();

        // Apply thermostat
        thermostat.applyThermostat(*particleContainer, iteration);

        // Calculate new velocities
        particleContainer->calculateV(delta_t);

        // Log current temperature
        double currentTemperature = thermostat.calculateCurrentTemperature(*particleContainer);
        spdlog::info("Current Temperature at iteration {}: {}", iteration, currentTemperature);

        // Plotting (if required)
        if (iteration % plotInterval == 0) {
            particleContainer->plotParticles(iteration);
        }

        // Update iteration and time
        iteration++;
        current_time += delta_t;
    }

    spdlog::info("Simulation completed");
    return 0;
}
