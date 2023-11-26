//
// Created by Layla Zadina on 19.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_LINKEDCELL_H
#define PSEMOLDYN_GROUPE_LINKEDCELL_H

#include <vector>
#include <array>
#include <unordered_map>
#include "Particle.h"
#include "ParticleContainer.h"

/** ParticleContainer holds all particles in the simulation.
 * It does not need to know about the cells.
 * LinkedCell manages the cell grid and the association of particles to cells.
 * It needs to be updated every time step to reflect particle movements. */


class LinkedCell {
    /** ArrayHash is efficient hashing for each cell
     * takes the cell std::array<int, 3>& array
     * returns int hash value of the cell
     * */
    struct ArrayHash {
        std::size_t operator()(const std::array<int, 3>& array) const {
            std::size_t hash = 17;
            hash = hash * 31 + std::hash<int>()(array[0]);
            hash = hash * 31 + std::hash<int>()(array[1]);
            hash = hash * 31 + std::hash<int>()(array[2]);
            return hash;
        }
    };

private:
    ParticleContainer& container;
    std::array<double, 3> domainSize;

    /** cells size = r_c × r_c × 1 for 2D and r_c × r_c × r_c for 3D
     * r_c = cutoffRadius */
    double cutoffRadius;


    /** cell maps the each cell of the grid to a list of particle indices that are currently
     * in that cell unordered_map is a hash table that stores key-value pairs
     * <std::array<int, 3> is key, represents a 3D grid index
     * std::vector<int> is value, stores indices of particles taken from ParticleContainer
     * provides efficient hashing and lookup in the unordered map
     * */
    std::unordered_map<std::array<int, 3>, std::vector<int>, ArrayHash> cells;


public:
    /** The constructor initializes the Linked Cell structure
     * with the given domain size, cutoff radius,
     * and a reference to the ParticleContainer.
     * */
    LinkedCell(std::array<double, 3> domainSize, double cutoffRadius, ParticleContainer& container);

    /** update the cell associations of all particles in the simulation at each time step
     * */
    void updateCells();

    /** calculates the cell index based on a particle's position */
    std::array<int, 3> calculateCellIndex(const Particle& particle) const;

    /** calculates interaction forces between particles in the same or neighboring cells */
    void calculateInteractions();

    std::vector<int> getParticlesInCell(const std::array<int, 3>& cellIndex) const;
    std::vector<int> getParticlesinNeighborCells(const std::array<int, 3>& cellIndex) const;
    bool isWithinDomain(const std::array<int, 3>& cellIndex) const;

};


#endif //PSEMOLDYN_GROUPE_LINKEDCELL_H
