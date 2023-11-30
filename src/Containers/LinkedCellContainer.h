//
// Created by sophy on 26.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER_H
#define PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER_H

#pragma once

#include "Particle.h"
#include <vector>
#include <array>
#include <unordered_map>
#include <cstddef> // for std::size_t
#include "ParticleContainerBase.h"


class LinkedCellContainer : public ParticleContainerBase {

    /** ArrayHash is efficient hashing for each cell
     * takes the cell std::array<int, 3>& array (x, y, z)
     * returns int hash value of the cell
     * in short: cell --> hash
     * */
    struct ArrayHash {
        std::size_t operator()(const std::array<int, 3> &array) const {
            std::size_t hash = 17;
            hash = hash * 31 + std::hash<int>()(array[0]);
            hash = hash * 31 + std::hash<int>()(array[1]);
            hash = hash * 31 + std::hash<int>()(array[2]);
            return hash;
        }
    };

private:

    /**
   * All the particles
   */
    std::vector<Particle> particleList;

    std::array<double, 3> domainSize;
    /**
     * cells size = r_c × r_c × 1 for 2D and r_c × r_c × r_c for 3D
     * r_c = cutoffRadius
     * */
    double cutoffRadius;

    /** cells maps the each cell of the grid to a list of particle indices that are currently
     * in that cell unordered_map is a hash table that stores key-value pairs
     * key : std::array<int, 3> represents a index of cell from the 3D grid
     * value : std::vector<int> stores indices of particles taken from particleList (before: ParticleContainer)
     * provides efficient hashing and lookup in the unordered map
     * */
    std::unordered_map<std::array<int, 3>, std::vector<int>, ArrayHash> cells;

public:
    LinkedCellContainer(ForceBase &model, std::array<double, 3> &domainSize, double &cutoffRadius,
                        std::vector<Particle> &particles,
                        std::unordered_map<std::array<int, 3>, std::vector<int>, ArrayHash> cellGrid)
            : ParticleContainerBase(model), domainSize(domainSize), cutoffRadius(cutoffRadius),
              particleList(particles), cells(cellGrid) {};

    LinkedCellContainer(ForceBase &model) : ParticleContainerBase(model), particleList(std::vector<Particle>()),
                                            domainSize({1, 1, 1}), cutoffRadius(1.0),
                                            cells(std::unordered_map<std::array<int, 3>, std::vector<int>, ArrayHash>()) {};

    LinkedCellContainer(ForceBase &model, std::vector<Particle> &particles)
            : ParticleContainerBase(model), particleList(particles), domainSize({1, 1, 1}), cutoffRadius(1.0),
              cells(std::unordered_map<std::array<int, 3>, std::vector<int>, ArrayHash>()) {};


    double &getCutoffRadius();

    std::vector<Particle> &getParticles();

    void addParticle(Particle &particle);

    std::size_t size() const;

    void resetF();

    void calculateF();

    void calculateX(double delta_t);

    void calculateV(double delta_t);

    void plotParticles(int iteration);

    bool isWithinDomain(const std::array<int, 3> &cellIndex);

    /**
     * calculates the cell index based on a particle's position
     * */
    std::array<int, 3> &calculateCellIndex(const Particle &particle);

    std::vector<int> &getParticlesInNeighborCells(const std::array<int, 3> &cellIndex);

    std::vector<int> &getParticlesInCell(const std::array<int, 3> &cellIndex);


    /**
     * calculates interaction forces between particles in the same or neighboring cells
     * */
    void calculateInteractions();

    /**
     * update the cell associations of all particles in the simulation at each time step
     * */
    void updateCells(double delta_t);
};

#endif //PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER_H
