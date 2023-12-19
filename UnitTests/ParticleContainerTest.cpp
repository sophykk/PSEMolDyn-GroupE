//
// Created by markus on 11/7/23.
//
#include <gtest/gtest.h>
#include "../src/Containers/BasicParticleContainer.h"
#include "../src/Containers/LinkedCellContainer.h"
#include "../src/Particle.h"
#include "../src/Forces/LennardJonesForce.h"
#include "../src/Formulas.h"
#include <array>
#include <cmath>

/** BasicParticleContainer Tests: */

// testing the default constructor of the BasicParticleContainer
TEST(ParticleContainerTest, DefaultConstructor) {
    LennardJonesForce force;
    BasicParticleContainer container(force);
    EXPECT_EQ(container.size(), 0);
}

// testing the constructor with a list as a parameter
TEST(ParticleContainerTest, ConstructorWithList) {
    LennardJonesForce force;
    std::vector<Particle> particlesList = {Particle(), Particle()};
    BasicParticleContainer container(force, particlesList);
    EXPECT_EQ(container.size(), particlesList.size());
}

// testing the adding function of the BasicParticleContainer
TEST(ParticleContainerTest, AddParticle) {
    LennardJonesForce force;
    BasicParticleContainer container(force);
    Particle particle;
    container.addParticle(particle);
    EXPECT_EQ(container.size(), 1);
}

/** LinkedCellContainer Tests: */

TEST(ParticleContainerTest, LinkedCellContainerDefaultConstructor) {
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 3.0;
    std::array<char, 4> boundaryCon{'r', 'r', 'r', 'r'}; // Reflecting
    double gGrav = 12;
    LennardJonesForce force;

    LinkedCellContainer container(force, domainSize, cutoffRadius, boundaryCon, gGrav);

    EXPECT_EQ(container.size(), 0);
    EXPECT_EQ(container.getBoundaryCon(1), boundaryCon[1]);
}

/**This test checks if the constructor correctly
 * initializes the container with a given list of particles*/
TEST(ParticleContainerTest, LinkedCellContainer2ConstructorWithParticles) {
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 3.0;
    std::array<char, 4> boundaryCon{'r', 'r', 'r', 'r'}; // Reflecting
    double gGrav = 12;
    LennardJonesForce force;
    std::vector<Particle> particlesList = {Particle(), Particle()};

    LinkedCellContainer container(force, particlesList, domainSize, cutoffRadius, boundaryCon, gGrav);

    EXPECT_EQ(container.size(), particlesList.size());
}

/** test addParticle method correctly adds a particle to the container*/
TEST(ParticleContainerTest, LinkedCellContainer2AddParticle) {
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 3.0;
    std::array<char, 4> boundaryCon{'r', 'r', 'r', 'r'}; // Reflecting
    double gGrav = 12;

    LennardJonesForce force;

    LinkedCellContainer container(force, domainSize, cutoffRadius, boundaryCon, gGrav);

    Particle particle;
    container.addParticle(particle);

    EXPECT_EQ(container.size(), 1);
}

/** test if the initGrid method correctly initializes the grid. */
TEST(ParticleContainerTest, LinkedCellContainer2InitGrid) {
    LennardJonesForce force;
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 10.0;
    std::array<char, 4> boundaryCon{'r', 'r', 'r', 'r'}; // Reflecting
    std::vector<Particle> particlesList = {
            Particle({5.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0 ,12.0, 1.0, 5.0, 1),
            Particle({15.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 12.0, 1.0,5.0, 1)
    };
    double gGrav = 12.0;

    LinkedCellContainer container(force, particlesList, domainSize, cutoffRadius, boundaryCon, gGrav);
    container.initGrid();

    auto grid = container.getGrid();
    EXPECT_EQ(grid.size(), 18); // Based on domainSize and cutoffRadius
    EXPECT_EQ(grid[0].size(), 9);
    EXPECT_EQ(grid[0][0].size(), 1);
    EXPECT_EQ(grid[1][0].size(), 1);
}

/** test if the calculateF correctly calculates and applies forces between 2 particles*/
TEST(ParticleContainerTest, LinkedCellContainer2CalculateForce) {
    LennardJonesForce force; // sigma = 1 and eps = 5
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 10.0;
    std::array<char, 4> boundaryCon{'r', 'r', 'r', 'r'}; // Reflecting
    std::vector<Particle> particlesList = {
            Particle({5.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 12.0, 1.0, 5.0, 1),
            Particle({15.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 12.0, 1.0, 5.0, 1)
    };
    double gGrav = 12.0;

    LinkedCellContainer container(force, particlesList, domainSize, cutoffRadius, boundaryCon, gGrav);
    container.calculateF();

    auto particles = container.getParticles();
    // Expected force calculated with external tool WolframAlpha using the Lennard-Jones Formula
    std::array<double, 3> expectedForce = {1.1999976e-05, 0.0, 0.0};

    // Check forces on the first particle are as expected
    EXPECT_NEAR(particles[0].getF()[0], expectedForce[0], 1e-7);
    EXPECT_NEAR(particles[0].getF()[1], expectedForce[1], 1e-7);
    EXPECT_NEAR(particles[0].getF()[2], expectedForce[2], 1e-7);

    // Check forces on the second particle suing N3L
    EXPECT_NEAR(particles[1].getF()[0], -expectedForce[0], 1e-7);
    EXPECT_NEAR(particles[1].getF()[1], -expectedForce[1], 1e-7);
    EXPECT_NEAR(particles[1].getF()[2], -expectedForce[2], 1e-7);
}

/** test for resetF */
TEST(ParticleContainerTest, LinkedCellContainer2ResetForces) {
    LennardJonesForce force;
    std::vector<double> domainSize = {100.0, 100.0, 1.0};
    double cutoffRadius = 5.0;
    std::array<char, 4> boundaryCon{'r', 'r', 'r', 'r'}; // Reflecting
    std::vector<Particle> particlesList = {/* Initialize some particles */};
    double gGrav = 12.0;

    LinkedCellContainer container(force, particlesList, domainSize, cutoffRadius, boundaryCon, gGrav);

    // Simulate a step to ensure some forces are calculated
    container.calculateF();
    container.resetF();

    for (const auto& particle : container.getParticles()) {
        // Check if the force on each particle is reset to zero
        auto force = particle.getF();
        EXPECT_EQ(force[0], 0.0);
        EXPECT_EQ(force[1], 0.0);
        EXPECT_EQ(force[2], 0.0);
    }
}

/** Test for getting Reflecting/Outflow Boundary Condition */
TEST(ParticleContainerTest, LinkedCellContainer2BoundaryConditions) {
    LennardJonesForce force;
    std::vector<double> domainSize = {100.0, 100.0, 1.0};
    double cutoffRadius = 5.0;
    std::array<char, 4> reflectingBoundaryCon{'r', 'r', 'r', 'r'}; // Reflecting
    std::array<char, 4> openBoundaryCon{'o', 'o', 'o', 'o'}; // Outflow
    double gGrav;

    LinkedCellContainer reflectingContainer(force, domainSize, cutoffRadius, reflectingBoundaryCon, gGrav);
    EXPECT_EQ(reflectingContainer.getBoundaryCon(0), reflectingBoundaryCon[0]);

    LinkedCellContainer openContainer(force, domainSize, cutoffRadius, openBoundaryCon, gGrav);
    EXPECT_EQ(openContainer.getBoundaryCon(2), openBoundaryCon[2]);

}













