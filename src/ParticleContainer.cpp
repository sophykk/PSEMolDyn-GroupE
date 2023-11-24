//
// Created by sophy on 02.11.2023.
//

#include "ParticleContainer.h"
#include <list>
#include <vector>
#include <utility> //for pair
#include <algorithm>
#include <cstddef>

ParticleContainer::ParticleContainer() {
    createParticlePairs();
}

ParticleContainer::ParticleContainer(std::vector <Particle> pVector) : particleList(pVector){
    createParticlePairs();
}

std::vector<Particle>& ParticleContainer::getParticles() {
    return particleList;
}


std::vector<std::pair<Particle&, Particle&>>& ParticleContainer::getParticlePairs() {
    return particlePairs;
}

void ParticleContainer::addParticle(Particle &particle) {
    particleList.push_back(particle);
}


void ParticleContainer::createParticlePairs() {
    for (auto& p1 : particleList) {
        for (auto& p2 : particleList) {
            if (&p1 != &p2) {
                std::pair<Particle&, Particle&> currentPair = std::make_pair(std::ref(p1), std::ref(p2));
                std::pair<Particle&, Particle&> reversePair = std::make_pair(std::ref(p2), std::ref(p1));

                if (std::find(particlePairs.begin(), particlePairs.end(), currentPair) == particlePairs.end() &&
                    std::find(particlePairs.begin(), particlePairs.end(), reversePair) == particlePairs.end()) {
                    particlePairs.push_back(currentPair);
                }
            }
        }
    }
}

std::size_t ParticleContainer::size() const {
    return particleList.size();
}



