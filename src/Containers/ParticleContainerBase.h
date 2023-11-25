//
// Created by sophy on 02.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_PARTICLECONTAINERBASE_H
#define PSEMOLDYN_GROUPE_PARTICLECONTAINERBASE_H

#pragma once

#include <list>
#include "Particle.h"
#include <vector>
#include <utility>
#include <cstddef> // for std::size_t
#include "utils/ArrayUtils.h"
#include "inputOutput/outputWriter/VTKWriter.h"
#include "Formulas.h"
#include "Forces/ForceBase.h"

class ParticleContainerBase {

protected:
    ForceBase& forceModel;

public:
    ParticleContainerBase(ForceBase& model) : forceModel(model) {}

    virtual void addParticle(Particle &particle) = 0;

    std::vector<Particle>& getParticles();

    virtual std::size_t size() const = 0;

    virtual void resetF() = 0;

    virtual void calculateF() = 0;

    virtual void calculateX(double delta_t) = 0;

    virtual void calculateV(double delta_t) = 0;

    virtual void plotParticles(int iteration) = 0;
};

#endif //PSEMOLDYN_GROUPE_PARTICLECONTAINERBASE_H
