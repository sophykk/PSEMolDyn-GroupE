//
// Created by sophy on 26.11.2023.
//

#include "LinkedCellContainer.h"
#include <list>
#include <vector>
#include "Formulas.h"
#include "utils/ArrayUtils.h"
#include "inputOutput/outputWriter/VTKWriter.h"
#include "Forces/ForceBase.h"


double &LinkedCellContainer::getCutoffRadius() {
    return cutoffRadius;
}

std::vector<Particle> &LinkedCellContainer::getParticles() {
    return particleList;
}

void LinkedCellContainer::addParticle(Particle &particle) {
    particleList.push_back(particle);
}

std::size_t LinkedCellContainer::size() const {
    return particleList.size();
}

void LinkedCellContainer::resetF() {
    for (auto &p: getParticles()) {
        p.resetF();
    }
}

void LinkedCellContainer::calculateF() {
    // Iterate over all cells
    for (const auto &cell: cells) {
        const auto &cellIndex = cell.first;
        const std::vector<int> &particlesInCell = cell.second;

        // Get neighboring particles
        std::vector<int> neighboringParticles = getParticlesInNeighborCells(cellIndex);

        // Iterate over all particles in the current cell
        for (int particleIndex: particlesInCell) {
            // Calculate force with neighboring particles
            for (int neighborIndex: neighboringParticles) {
                // Avoid double calculation and self-interaction
                if (neighborIndex <= particleIndex) continue;
                // Calculate force
                std::array<double, 3> force = forceModel.calculateForce(getParticles()[particleIndex],
                                                                        getParticles()[neighborIndex]);
                getParticles()[particleIndex].addF(force);
                getParticles()[neighborIndex].addF(-1.0 * force);
                //set/add force to the particles done in calculateF? but do here, so it doesn't calculate for whole list
            }
        }
    }
    /*
    for (auto p1 = getParticles().begin(); p1 < getParticles().end(); p1++) {
        for (auto p2 = p1 + 1; p2 < getParticles().end(); p2++) {
            auto force = forceModel.calculateForce(*p1, *p2);
            p1->addF(force);
            p2->addF(-1.0 * force);
        }
    }*/
}

void LinkedCellContainer::calculateX(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: xi(tn+1) = xi(tn) + ∆t · vi(tn) + (∆t)^2 * (Fi(tn)/2mi)
        // Calculate xi(tn+1)
        auto xi_tn1 = Formulas::verletXStep(p.getX(), p.getV(), p.getF(), p.getM(), delta_t);
        p.setX(xi_tn1);  // Update the position
    }
}

void LinkedCellContainer::calculateV(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: vi(tn+1) = vi(tn) + ∆t * ((Fi(tn) + Fi(tn+1))/ 2mi)
        // Calculate the velocity at time tn+1
        auto vi_tn1 = Formulas::verletVStep(p.getV(), p.getOldF(), p.getF(), p.getM(), delta_t);
        p.setV(vi_tn1);
    }
}

void LinkedCellContainer::plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(size());
    for (auto &p: getParticles()) {
        writer.plotParticle(p);
    }
    writer.writeFile(out_name, iteration);
}

/**
 * @param cellIndex is a vector describing the position of that cell in the grid (x, y, z)
 * @return true if cell is in grid, else false
 */

bool LinkedCellContainer::isWithinDomain(const std::array<int, 3> &cellIndex) {
    for (int i = 0; i < 3; ++i) {
        if (cellIndex[i] < 0 || cellIndex[i] >= domainSize[i]) {
            //in the cell what do I refer to as index? Left corner/center, doesn't matter???? just (x, y, z)
            return false;
        }
    }
    return true;
}

/**
 * calculates the cell index based on a particle's position
 * attention: ref to local variable index returned
 * */
std::array<int, 3> &LinkedCellContainer::calculateCellIndex(const Particle &particle) {
    std::array<int, 3> index;
    auto pos = particle.getX(); // Get the position of the particle

    // Calculate the index in each dimension
    for (int i = 0; i < 3; ++i) {
        auto x = pos[i] / cutoffRadius;
        index[i] = static_cast<int>(pos[i]);
    }
    if (isWithinDomain(index)) {
        return index;
    }
    else{
        //if particle in outside cell
        std::array<int, 3> null{-1, -1, -1};
        return null;
    }

}

/**
 * @param cellIndex is a vector describing the position of that cell in the grid (x, y, z)
 * @return vector of indices of all Particles in that cell
 * Get particles in a specific cell and uses the hash values of the cell
 * returns the vector of Particle indices (int) in the container
 */
std::vector<int> &LinkedCellContainer::getParticlesInCell(const std::array<int, 3> &cellIndex) {

    auto cellHash = cells.find(cellIndex);
    if (cellHash != cells.end()) {
        return cellHash->second; // Returns the vector of particle indices in the cell
    } else {
        std::vector<int> pIndices{};
        return pIndices; // Returns an empty vector if the cell is not found
    }

}

/**
 * Get particles (represented by indices) in neighboring cells of a specific cell
 * @param cellIndex is a vector containing the indices of all particles in that cell
 * @return vector containing the indices of all particles from surrounding cells
 * attention here too: ref to local variable index returned
 * */

std::vector<int> &LinkedCellContainer::getParticlesInNeighborCells(const std::array<int, 3> &cellIndex) {
    std::vector<int> neighboringParticles;

    // Iterate over neighboring cells in 3D (3x3x3 grid)
    // offsets dx, dy, dz range from -1 to 1, covering all adjacent cells.
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                if (dx == 0 && dy == 0 && dz == 0) continue; // Skip the cell itself, is in the middle

                std::array<int, 3> neighborIndex = {cellIndex[0] + dx, cellIndex[1] + dy, cellIndex[2] + dz};

                // Check if neighbor index is within domain bounds
                if (isWithinDomain(neighborIndex)) {
                    auto found = cells.find(neighborIndex);
                    if (found != cells.end()) {
                        const std::vector<int> &particlesInNeighbor = found->second;
                        neighboringParticles.insert(neighboringParticles.end(), particlesInNeighbor.begin(),
                                                    particlesInNeighbor.end());
                    }
                }
            }
        }
    }
    return neighboringParticles;
}

/**
 * calculates interaction forces between particles in the same or neighboring cells
 * done in calculateF now
 * */
void LinkedCellContainer::calculateInteractions() {
// Constants for Lennard-Jones potential (could be parameters or class members)

    // Iterate over all cells
    for (const auto &cell: cells) {
        const auto &cellIndex = cell.first;
        const std::vector<int> &particlesInCell = cell.second;

        // Get neighboring particles
        std::vector<int> neighboringParticles = getParticlesInNeighborCells(cellIndex);

        // Iterate over all particles in the current cell
        for (int particleIndex: particlesInCell) {
            // Calculate force with neighboring particles
            for (int neighborIndex: neighboringParticles) {
                // Avoid double calculation and self-interaction
                if (neighborIndex <= particleIndex) continue;
                // Calculate force
                std::array<double, 3> force = forceModel.calculateForce(getParticles()[particleIndex],
                                                                        getParticles()[neighborIndex]);
                getParticles()[particleIndex].addF(force);
                getParticles()[neighborIndex].addF(-1.0 * force);
                //set/add force to the particles done in calculateF? but do here, so it doesn't calculate for whole list
            }
        }
    }
}

/**
  * update the cell associations of all particles in the simulation at each time step
  * actually done in MolSim by many steps like shown, but doesn't work yet, in MolSim better
* */
void updateCells(double delta_t) {
    /* calculateX(delta_t);

    // reset forces
    resetF();

    // calculate new f
    calculateF();

    // calculate new v
    calculateV(delta_t); */
}
