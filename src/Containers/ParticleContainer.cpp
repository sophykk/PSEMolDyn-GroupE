//
// Created by sophy on 02.11.2023.
//
#include "ParticleContainer.h"
#include <list>
#include <vector>
#include "Formulas.h"

std::vector<Particle>& ParticleContainer::getParticles() {
    return particleList;
}

void ParticleContainer::addParticle(Particle &particle) {
    particleList.push_back(particle);
}

std::size_t ParticleContainer::size() const {
    return particleList.size();
}

void ParticleContainer::resetF() {
    for (auto &p: getParticles()) {
        p.resetF();
    }
}

void ParticleContainer::calculateF(){
    for (auto p1 = getParticles().begin(); p1 < getParticles().end(); p1++) {
        for (auto p2 = p1 + 1; p2 < getParticles().end(); p2++) {
            auto force = forceModel.calculateForce(*p1, *p2);
            p1->addF(force);
            p2->addF(-1.0 * force);
        }
    }
}

void ParticleContainer::calculateX(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: xi(tn+1) = xi(tn) + ∆t · vi(tn) + (∆t)^2 * (Fi(tn)/2mi)
        // Calculate xi(tn+1)
        auto xi_tn1 = Formulas::verletXStep(p.getX(), p.getV(), p.getF(), p.getM(), delta_t);
        p.setX(xi_tn1);  // Update the position
    }
}

void ParticleContainer::calculateV(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: vi(tn+1) = vi(tn) + ∆t * ((Fi(tn) + Fi(tn+1))/ 2mi)
        // Calculate the velocity at time tn+1
        auto vi_tn1 = Formulas::verletVStep(p.getV(), p.getOldF(), p.getF(), p.getM(), delta_t);
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