//
// Created by sophy on 02.11.2023.
//

#include "ParticleContainer.h"
#include <list>
#include <vector>
#include <utility> //for pair
#include <algorithm>
#include <cstddef>
#include <iostream>

ParticleContainer::ParticleContainer() = default;

std::vector<Particle>& ParticleContainer::getParticles() {
    return particleList;
}


std::vector<std::pair<Particle*, Particle*>>& ParticleContainer::getParticlePairs() {
    return particlePairs;
}

void ParticleContainer::addParticle(Particle &particle) {
    particleList.push_back(particle);
}


/*void ParticleContainer::createParticlePairs() {
    //std::cout<<"Inside the function";
    for (auto& p1 : particleList) {
        for (auto& p2 : particleList) {
            if (&p1 != &p2) {
                std::pair<Particle*, Particle*> currentPair = std::make_pair(&p1, &p2);
                std::pair<Particle*, Particle*> reversePair = std::make_pair(&p2, &p1);

                if (std::find(particlePairs.begin(), particlePairs.end(), currentPair) == particlePairs.end() &&
                    std::find(particlePairs.begin(), particlePairs.end(), reversePair) == particlePairs.end()) {
                    particlePairs.push_back(currentPair);
                }
            }
        }
        //std::cout<<"doing something";
    }
}*/

void ParticleContainer::createParticlePairs() {
    for (int i = 0; i < particleList.size(); ++i) {
        for (int j = i + 1; j < particleList.size(); ++j) {
            // Check if the pair is not repeated
            if (i != j) {
                // Create a pair and add it to the vector
                particlePairs.emplace_back(&particleList[i], &particleList[j]);
            }
        }
    }
}

std::size_t ParticleContainer::size() const {
    return particleList.size();
}



