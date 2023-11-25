//
// Created by sophy on 02.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_PARTICLECONTAINER_H
#define PSEMOLDYN_GROUPE_PARTICLECONTAINER_H

#pragma once

#include <list>
#include "Particle.h"
#include <vector>
#include <utility>
#include <cstddef> // for std::size_t
#include "utils/ArrayUtils.h"
#include "inputOutput/outputWriter/VTKWriter.h"
#include "ParticleContainerBase.h"


class ParticleContainer : public ParticleContainerBase {
private:

    /**
   * All the particles
   */
    std::vector <Particle> particleList;

public:

    ParticleContainer(ForceBase& model) : ParticleContainerBase(model){};

    void addParticle(Particle &particle);

    std::vector<Particle>& getParticles();

    std::size_t size() const;

    void resetF();

    void calculateF();

    void calculateX(double delta_t);

    void calculateV(double delta_t);

    void plotParticles(int iteration);
};

#endif //PSEMOLDYN_GROUPE_PARTICLECONTAINER_H
