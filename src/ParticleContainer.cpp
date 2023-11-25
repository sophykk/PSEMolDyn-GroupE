//
// Created by sophy on 02.11.2023.
//

#include "ParticleContainer.h"
#include <list>
#include <vector>
#include <utility> //for pair
#include <algorithm>
#include <cstddef>
#include "Formulas.h"

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

void ParticleContainer::resetF() {
    for (auto &p: getParticles()) {
        p.resetF();
    }
}

void ParticleContainer::calculateF(double sigma, double eps){
    for (auto &p: getParticles()) {
        //std::array<double, 3> F_i{0., 0., 0.};
        for (auto &p2: getParticles()) {
            // formula: Fij = ((mi * mj) / ||xi −xj||^3) * (xj − xi)
            //std::array<double, 3> F_ij{};
            if (&p != &p2) {
                auto force = (-24 * eps / pow(Formulas::secondNorm(p.getX() - p2.getX()), 2.0)) *
                             (pow((sigma / Formulas::secondNorm(p.getX() - p2.getX())), 6.0) -
                              2 * pow((sigma / Formulas::secondNorm(p.getX() - p2.getX())), 12.0)) * (p.getX() - p2.getX());
                //std::cout<<"force:"<<force;
                p.addF(force);
                //std::cout<<p;
            }
        }

    }

}

void ParticleContainer::calculateX(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: xi(tn+1) = xi(tn) + ∆t · vi(tn) + (∆t)^2 * (Fi(tn)/2mi)
        // Calculate xi(tn+1)

        auto xi_tn1 = p.getX() + delta_t * p.getV() + (delta_t * delta_t) / (2.0 * p.getM()) * p.getF();

        p.setX(xi_tn1);  // Update the position
    }
}

void ParticleContainer::calculateV(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: vi(tn+1) = vi(tn) + ∆t * ((Fi(tn) + Fi(tn+1))/ 2mi)
        // Calculate the velocity at time tn+1

        auto vi_tn1 = p.getV() + delta_t / (2 * p.getM()) * (p.getOldF() + p.getF());
        p.setV(vi_tn1);
    }
}

void ParticleContainer::plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(size());
    for (auto &p: getParticles()) {
        writer.plotParticle(p);
    }
    writer.writeFile(out_name, iteration);
}