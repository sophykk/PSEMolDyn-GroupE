//
// Created by sophy on 02.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_BASICPARTICLECONTAINER_H
#define PSEMOLDYN_GROUPE_BASICPARTICLECONTAINER_H

#pragma once

#include "Particle.h"
#include <vector>
#include <cstddef> // for std::size_t
#include "ParticleContainerBase.h"


class BasicParticleContainer : public ParticleContainerBase {
private:

    /**
   * All the particles
   */
    std::vector<Particle> particleList;

public:

    BasicParticleContainer(ForceBase &model) : ParticleContainerBase(model) {};

    BasicParticleContainer(ForceBase &model, std::vector<Particle> &particles) : particleList(particles),
                                                                                 ParticleContainerBase(model) {};

    void addParticle(Particle &particle);

    std::vector<Particle> &getParticles();

    int getDimension() const;

    std::vector<double> getDomainSize();

    std::size_t size() const;

    void resetF();

    void calculateF();

    void calculateX(double delta_t);

    void calculateV(double delta_t);

    void plotParticles(int iteration);
};

#endif //PSEMOLDYN_GROUPE_BASICPARTICLECONTAINER_H
