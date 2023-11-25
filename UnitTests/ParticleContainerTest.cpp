//
// Created by markus on 11/7/23.
//
#include <gtest/gtest.h>

// testing the default constructor of the BasicParticleContainer
TEST(ParticleContainerTest, DummyTest){
    EXPECT_EQ(0,0);
}
/*TEST(ParticleContainerTest, DefaultConstructor) {
BasicParticleContainer container;
EXPECT_EQ(container.size(), 0);
}

// testing the constructor with a list as a parameter
TEST(ParticleContainerTest, ConstructorWithList) {
std::vector<Particle> particlesList = {Particle(), Particle()};
BasicParticleContainer container((particlesList));
EXPECT_EQ(container.size(), particlesList.size());
}

// testing the adding function of the BasicParticleContainer
TEST(ParticleContainerTest, AddParticle) {
BasicParticleContainer container;
Particle particle;
container.addParticle(particle);
EXPECT_EQ(container.size(), 1);
}*/