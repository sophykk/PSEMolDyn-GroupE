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


class LinkedCell {
    
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
    double cutoffRadius;
    std::unordered_map<std::array<int, 3>, std::vector<int>, ArrayHash> cells;


public:
    LinkedCell(std::array<double, 3> domainSize, double cutoffRadius, ParticleContainer& container);

    void populateCells();
    void calculateInteractions();
    std::vector<int> getParticlesInCell(const std::array<int, 3>& cellIndex) const;
    std::vector<int> getNeighboringCells(const std::array<int, 3>& cellIndex) const;

};


#endif //PSEMOLDYN_GROUPE_LINKEDCELL_H
