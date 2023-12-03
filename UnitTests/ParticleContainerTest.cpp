//
// Created by markus on 11/7/23.
//
#include <gtest/gtest.h>
#include "../src/Containers/BasicParticleContainer.h"
#include "../src/Containers/LinkedCellContainer2.h"
#include "../src/Particle.h"
#include "../src/Forces/LennardJonesForce.h"

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


// testing the default constructor of the LinkedCellContainer
TEST(ParticleContainerTestLC, DefaultConstructor) {
    LennardJonesForce force;
    std::vector<double> dSize;
    double cRad;
    char b;
    LinkedCellContainer2 container(force, dSize, cRad, b);
    EXPECT_EQ(container.size(), 0);
}

// testing the constructor with a list as a parameter
TEST(ParticleContainerTestLC, ConstructorWithList) {
    LennardJonesForce force;
    std::vector<Particle> particlesList = {Particle(), Particle()};
    std::vector<double> dSize;
    double cRad;
    char b;
    LinkedCellContainer2 container(force, particlesList, dSize, cRad, b);
    EXPECT_EQ(container.size(), particlesList.size());
}

// testing the adding function of the LinkedCellContainer
TEST(ParticleContainerTestLC, AddParticle) {
    LennardJonesForce force;
    std::vector<double> dSize;
    double cRad;
    char b;
    LinkedCellContainer2 container(force, dSize, cRad, b);
    Particle particle;
    container.addParticle(particle);
    EXPECT_EQ(container.size(), 1);
}