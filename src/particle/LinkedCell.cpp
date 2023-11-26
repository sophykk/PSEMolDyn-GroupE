//
// Created by Layla Zadina on 19.11.2023.

#include <vector>
#include <array>
#include <unordered_map>
#include "Particle.h"
#include "ParticleContainer.h"
#include "../Formulas.h"

#include "LinkedCell.h"

// Constructor
LinkedCell::LinkedCell(std::array<double, 3> domainSize, double cutoffRadius, ParticleContainer& container)
        : domainSize(domainSize), cutoffRadius(cutoffRadius), container(container) {

}




/** Calculate force interactions between particles*/
void LinkedCell::calculateInteractions() {
    // Constants for Lennard-Jones potential (could be parameters or class members)
    //todo: should they be assigned here?
    double sigma = 1.0;
    double eps = 5.0;

    // Iterate over all cells
    for (const auto& cell : cells) {
        const auto& cellIndex = cell.first;
        const std::vector<int>& particlesInCell = cell.second;

        // Get neighboring particles
        std::vector<int> neighboringParticles = getParticlesinNeighborCells(cellIndex);

        // Iterate over all particles in the current cell
        for (int particleIndex : particlesInCell) {
            const Particle& particle = container.getParticles()[particleIndex];
            std::array<double, 3> particlePosition = particle.getX();
            // Calculate force with neighboring particles
            for (int neighborIndex : neighboringParticles) {
                // Avoid double calculation and self-interaction
                if (neighborIndex <= particleIndex) continue;

                const Particle& neighborParticle = container.getParticles()[neighborIndex];
                std::array<double, 3> neighborPosition = neighborParticle.getX();
                // Calculate the Lennard-Jones force
                std::array<double, 3> force = Formulas::calculateLJForce(particlePosition, neighborPosition, sigma, eps);

                //todo: set/add force to the particles
            }
        }
    }
}


void LinkedCell::updateCells() {
    // Clear the current cell associations
    for (auto& cell : cells) {
        cell.second.clear();
    }

    // Iterate through all particles and update their cell associations
    const auto& allParticles = container.getParticles();
    for (std::size_t i = 0; i < allParticles.size(); ++i) {
        const auto& particle = allParticles[i];
        auto cellIndex = calculateCellIndex(particle);

        // Check if the cell index is within the domain
        if (isWithinDomain(cellIndex)) {
            cells[cellIndex].push_back(i); // Add particle index to the cell
        }
        // Todo: If the particle is outside the domain, handle accordingly
        //  (e.g., ignore, apply boundary conditions, etc.)
    }
}


/** calculates the cell index for a given particle based on its position */
std::array<int, 3> LinkedCell::calculateCellIndex(const Particle& particle) const {
    std::array<int, 3> index;
    auto pos = particle.getX(); // Get the position of the particle

    // Calculate the index in each dimension
    for (int i = 0; i < 3; ++i) {
        index[i] = static_cast<int>(pos[i] / cutoffRadius);
    }
    //Todo: isWithinDomain
    return index;
}


/**Get particles in a specific cell
 * uses the hash values of the cell
 * returns the vector of Particle indices (int) in the container*/

std::vector<int> LinkedCell::getParticlesInCell(const std::array<int, 3>& cellIndex) const {
    auto cellHash = cells.find(cellIndex);
    if (cellHash != cells.end()) {
        return cellHash->second; // Returns the vector of particle indices in the cell
    } else {
        return std::vector<int>(); // Returns an empty vector if the cell is not found
    }
}


/**Get particles (represented by indices) in neighboring cells of a specific cell */
std::vector<int> LinkedCell::getParticlesinNeighborCells(const std::array<int, 3>& cellIndex) const {
    std::vector<int> neighboringParticles;

    // Iterate over neighboring cells in 3D (3x3x3 grid)
    // offsets dx, dy, dz range from -1 to 1, covering all adjacent cells.
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                if (dx == 0 && dy == 0 && dz == 0) continue; // Skip the cell itself

                std::array<int, 3> neighborIndex = {cellIndex[0] + dx, cellIndex[1] + dy, cellIndex[2] + dz};

                // Check if neighbor index is within domain bounds
                if (isWithinDomain(neighborIndex)) {
                    auto found = cells.find(neighborIndex);
                    if (found != cells.end()) {
                        const std::vector<int>& particlesInNeighbor = found->second;
                        neighboringParticles.insert(neighboringParticles.end(), particlesInNeighbor.begin(), particlesInNeighbor.end());
                    }
                }
            }
        }
    }

    return neighboringParticles;
}

bool LinkedCell::isWithinDomain(const std::array<int, 3>& cellIndex) const {
    for (int i = 0; i < 3; ++i) {
        if (cellIndex[i] < 0 || cellIndex[i] >= domainSize[i]) {
            //TODO: in the cell what do I refer to as index? Left corner/center
            return false;
        }
    }
    return true;
}


