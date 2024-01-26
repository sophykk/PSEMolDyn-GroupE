//
// Created by sophy on 02.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_PARTICLECONTAINERBASE_H
#define PSEMOLDYN_GROUPE_PARTICLECONTAINERBASE_H

#include <vector>
#include <cstddef> // for std::size_t
#include "Forces/ForceBase.h"
#include "Particle.h"

/**
 * Base class for particle containers
 * All particle containers should inherit from this class
 */
class ParticleContainerBase {

protected:
    /**
     * Reference to the used force model, e.g. gravitation or lennard-jones
     */
    ForceBase &forceModel;

public:
    /**
     * Abstract base class for all particle containers.
     * @param model Force model to be used for calculating forces on all particles.
     */
    ParticleContainerBase(ForceBase &model) : forceModel(model) {}

    /**
     * Default Destructor
     */
    virtual ~ParticleContainerBase() = default;

    /**
     * Add a particle to the container.
     * @param particle reference to the particle to add
     */
    virtual void addParticle(Particle &particle) = 0;

    /**
     * Get a flat vector of all particles
     * @return Vector of all particles
     */
    virtual std::vector<Particle> &getParticles() = 0;

    /**
     * Get dimension
     */
    virtual int getDimension() const = 0;

    /**
     * @return Number of particles in the container.
     */
    virtual std::size_t size() const = 0;

    /**
     * Swap current particle forces to old particle forces.
     * Set current particle forces to zero.
     */
    virtual void resetF() = 0;

    /**
     * Calculate particle forces.
     */
    virtual void calculateF() = 0;

    /**
     * Calculate new particle positions.
     * @param delta_t
     */
    virtual void calculateX(double delta_t) = 0;

    /**
     * Calculate new particle velocities.
     * @param delta_t
     */
    virtual void calculateV(double delta_t) = 0;

    /**
     * Write a VTK file with current particle information.
     * @param iteration Current simulation iteration
     */
    virtual void plotParticles(int iteration) = 0;
};

#endif //PSEMOLDYN_GROUPE_PARTICLECONTAINERBASE_H
