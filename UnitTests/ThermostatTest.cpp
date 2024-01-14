//
// Created by Layla Zadina on 19.12.2023.
//
#include <gtest/gtest.h>
#include <cmath>

#include "Thermostat.h" // Replace with your actual thermostat header file
#include "Containers/ParticleContainerBase.h"
#include "Containers/LinkedCellContainer.h"



class ThermostatTest : public testing::Test {
protected:
    static double calculateCurrentTemperature(const LinkedCellContainer& container) {
        // Assuming the temperature calculation is based on kinetic energy
        double totalKineticEnergy = 0.0;
        for (const auto& particle : container.getParticles()) {
            auto velocity = particle.getV();
            double speedSquared = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
            totalKineticEnergy += 0.5 * particle.getM() * speedSquared;
        }

        size_t particleCount = container.size();
        return (2 * totalKineticEnergy) / (3 * particleCount); // 3 for 3 degrees of freedom per particle
    }

TEST(ThermostatTest, IncreaseTemperatureTest) {
double initTemperature = 5.0;
double targetTemperature = 15.0;
double deltaT = 2.0;
double eps = 0.001;

Thermostat thermostat(initTemperature, targetTemperature, deltaT);
ParticleContainerBase particles;
initializeParticles(particles);

double startTemp = calculateCurrentTemperature(particles);
thermostat.applyThermostat(particles, 0);
double endTemp = calculateCurrentTemperature(particles);

EXPECT_NEAR(endTemp, std::min(startTemp + deltaT, targetTemperature), eps);
}

TEST(ThermostatTest, DecreaseTemperatureTest) {
double initTemperature = 20.0;
double targetTemperature = 10.0;
double deltaT = 3.0;
double eps = 0.001;

Thermostat thermostat(initTemperature, targetTemperature, deltaT);
ParticleContainerBase particles;
initializeParticles(particles);

double startTemp = calculateCurrentTemperature(particles);
thermostat.applyThermostat(particles, 0); // Assuming the thermostat is applied immediately
double endTemp = calculateCurrentTemperature(particles);

EXPECT_NEAR(endTemp, std::max(startTemp - deltaT, targetTemperature), eps);
}

TEST(ThermostatTest, HoldingTemperatureTest) {
double initTemperature = 10.0;
double targetTemperature = 10.0;
double deltaT = 2.0;
double eps = 0.001;

Thermostat thermostat(initTemperature, targetTemperature, deltaT);
ParticleContainerBase particles;
initializeParticles(particles);

double startTemp = calculateCurrentTemperature(particles);
thermostat.applyThermostat(particles, 0);
double endTemp = calculateCurrentTemperature(particles);

EXPECT_NEAR(endTemp, startTemp, eps);

}

