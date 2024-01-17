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


class LinkedCellContainer : public ParticleContainerBase {
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
    std::array<char, 4> boundaryCon;
    std::vector<Particle> haloList;
    double gGrav;

public:

    LinkedCellContainer(ForceBase &model, std::vector<double> &dSize, double &cRadius, std::array<char, 4> bCon, double &gGrav);

    LinkedCellContainer(ForceBase &model, std::vector<Particle> &particles, std::vector<double> &dSize,
                        double &cRadius, std::array<char, 4> bCon, double &gGrav);

    char &getBoundaryCon(int index);

    void setBoundaryCon(std::array<char, 4> &boundary);

    bool checkBoundary(char b);

    void addParticle(Particle &particle);

    std::vector<Particle> &getParticles();

    std::vector<std::vector<std::vector<Particle>>> &getGrid();

    std::pair<int, int> &getCell(const std::array<double, 3> &pos);

    //just 2D first
    //std::vector<Particle> &getParticlesFromCell(int x, int y);

    std::size_t size() const;

    void resetF();

    void calculateF();

    void harmonicPotential(Particle &p1, Particle &p2, double k, double r0);

    //  bool withinCutoff(Particle &p1, Particle &p2) const;

    void calculateX(double delta_t);

    void calculateV(double delta_t);

    void plotParticles(int iteration);

    void initGrid();

    bool checkDistance(Particle &p, std::string border);

    void applyReflecting(Particle &p);

    void applyPeriodic(Particle &p, int x, int y);
};

#endif //PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER2_H
