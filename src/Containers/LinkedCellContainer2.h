//
// Created by sophy on 30.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER2_H
#define PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER2_H
#pragma once

#include "Particle.h"
#include <vector>
#include <cstddef> // for std::size_t
#include "ParticleContainerBase.h"


class LinkedCellContainer2 : public ParticleContainerBase {
private:

    /**
   * All the particles
   */
    std::vector<Particle> particleList;
    std::vector<std::vector<std::vector<Particle>>> grid;
    //like {180, 90, 1}
    std::vector<double> domainSize;
    //like 3.0
    double cutoffRadius;
    char boundaryCon;

public:

    LinkedCellContainer2(ForceBase &model, std::vector<double> &dSize, double &cRadius, char bCon);

    LinkedCellContainer2(ForceBase &model, std::vector<Particle> &particles, std::vector<double> &dSize,
                         double &cRadius, char bCon);

    char &getBoundaryCon();

    void setBoundaryCon(char boundary);

    void addParticle(Particle &particle);

    std::vector<Particle> &getParticles();

    std::vector<std::vector<std::vector<Particle>>> &getGrid();

    //just 2D first
    //  std::vector<Particle> &getParticlesFromCell(int x, int y);

    std::size_t size() const;

    void resetF();

    void calculateF();

    //  bool withinCutoff(Particle &p1, Particle &p2) const;

    void calculateX(double delta_t);

    void calculateV(double delta_t);

    void plotParticles(int iteration);

    void initGrid();

    void applyReflecting(Particle &p);
};

#endif //PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER2_H
