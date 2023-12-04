//
// Created by markus on 11/7/23.
//
#include <gtest/gtest.h>
#include "../src/Containers/BasicParticleContainer.h"
#include "../src/Containers/LinkedCellContainer2.h"
#include "../src/Particle.h"
#include "../src/Forces/LennardJonesForce.h"
#include "../src/Formulas.h"

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

/** LinkedCellContainer2 Tests: */

TEST(ParticleContainerTest, LinkedCellContainer2DefaultConstructor) {
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 3.0;
    char boundaryCon = 'r'; // Reflecting boundary condition
    LennardJonesForce force;

    LinkedCellContainer2 container(force, domainSize, cutoffRadius, boundaryCon);

    EXPECT_EQ(container.size(), 0);
    EXPECT_EQ(container.getBoundaryCon(), boundaryCon);
}

/**This test checks if the constructor correctly
 * initializes the container with a given list of particles*/
TEST(ParticleContainerTest, LinkedCellContainer2ConstructorWithParticles) {
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 3.0;
    char boundaryCon = 'r';
    LennardJonesForce force;
    std::vector<Particle> particlesList = {Particle(), Particle()};

    LinkedCellContainer2 container(force, particlesList, domainSize, cutoffRadius, boundaryCon);

    EXPECT_EQ(container.size(), particlesList.size());
}

/** test addParticle method correctly adds a particle to the container*/
TEST(ParticleContainerTest, LinkedCellContainer2AddParticle) {
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 3.0;
    char boundaryCon = 'r';
    LennardJonesForce force;

    LinkedCellContainer2 container(force, domainSize, cutoffRadius, boundaryCon);

    Particle particle;
    container.addParticle(particle);

    EXPECT_EQ(container.size(), 1);
}

/** test if the initGrid method correctly initializes the grid. */
TEST(ParticleContainerTest, LinkedCellContainer2InitGrid) {
    LennardJonesForce force;
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 10.0;
    char boundaryCon = 'r';
    std::vector<Particle> particlesList = {
            Particle({5.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 1),
            Particle({15.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 1)
    };

    LinkedCellContainer2 container(force, particlesList, domainSize, cutoffRadius, boundaryCon);
    container.initGrid();

    auto grid = container.getGrid();
    EXPECT_EQ(grid.size(), 18); // Based on domainSize and cutoffRadius
    EXPECT_EQ(grid[0].size(), 9);
    EXPECT_EQ(grid[0][0].size(), 1);
    EXPECT_EQ(grid[1][0].size(), 1);
}

/** test if the calculateF correctly calculates and applies forces between 2 particles*/
TEST(ParticleContainerTest, LinkedCellContainer2CalculateForce) {
    LennardJonesForce force(1.0, 5.0); // sigma = 1 and eps = 5
    std::vector<double> domainSize = {180.0, 90.0, 1.0};
    double cutoffRadius = 10.0;
    char boundaryCon = 'r';
    std::vector<Particle> particlesList = {
            Particle({5.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 1),
            Particle({15.0, 5.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 1)
    };

    LinkedCellContainer2 container(force, particlesList, domainSize, cutoffRadius, boundaryCon);
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
    char boundaryCon = 'r';
    std::vector<Particle> particlesList = {/* Initialize some particles */};

    LinkedCellContainer2 container(force, particlesList, domainSize, cutoffRadius, boundaryCon);

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

/** Test for getting Refelcting/Outflow Boundary Condition */
TEST(ParticleContainerTest, LinkedCellContainer2BoundaryConditions) {
    LennardJonesForce force;
    std::vector<double> domainSize = {100.0, 100.0, 1.0};
    double cutoffRadius = 5.0;
    char reflectingBoundaryCon = 'r'; // Reflecting
    char openBoundaryCon = 'o'; // Outflow

    LinkedCellContainer2 reflectingContainer(force, domainSize, cutoffRadius, reflectingBoundaryCon);
    EXPECT_EQ(reflectingContainer.getBoundaryCon(), reflectingBoundaryCon);

    LinkedCellContainer2 openContainer(force, domainSize, cutoffRadius, openBoundaryCon);
    EXPECT_EQ(openContainer.getBoundaryCon(), openBoundaryCon);

}












