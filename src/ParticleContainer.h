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
#include "outputWriter/VTKWriter.h"


class ParticleContainer {

private:

    /**
   * All the particles
   */
    std::vector <Particle> particleList;

    /**
  * Pairs of all particles
  * a.e. for 4 particles p1,p2,p3,p4, there are 6 particle pairs: (p1,p2),(p1,p3),(p1,p4),(p2,p3),(p2,p4),(p3,p4).
  */
    std::vector <std::pair<Particle&, Particle&>> particlePairs;

public:

    ParticleContainer();

    ParticleContainer(std::vector <Particle> pVector);

    void addParticle(Particle &particle);

    std::vector<Particle>& getParticles();

    std::vector<std::pair<Particle&, Particle&>>& getParticlePairs();

    std::size_t size() const;

    void createParticlePairs();

    void calculateF(double sigma, double eps);

    void calculateX(double delta_t);

    void calculateV(double delta_t);

    void plotParticles(int iteration);
};

#endif //PSEMOLDYN_GROUPE_PARTICLECONTAINER_H
