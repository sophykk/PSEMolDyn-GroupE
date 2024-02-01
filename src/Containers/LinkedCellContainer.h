//
// Created by sophy on 30.11.2023.
//

#ifndef PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER_H
#define PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER_H
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
    std::vector<std::vector<std::vector<std::vector<Particle>>>> grid;
    //like {180, 90, 1} -> always z = 1 for 2D
    std::vector<double> domainSize;
    //like 3.0
    double cutoffRadius;
    std::array<char, 6> boundaryCon;
    //std::vector<Particle> haloList;
    double gGrav;
    //0 for gGrav along x-axis, 1 for y-axis, 2 for z-axis
    //i.e 1 for w5t3.xml & w5t4.xml and 2 for w5t1.xml
    int gGravPos = 1;
    bool isMembrane;
    int k = 300;
    double r0 = 2.2;
    double pullUpF = -0.8;


public:

    LinkedCellContainer(ForceBase &model, std::vector<double> &dSize, double &cRadius, std::array<char, 6> bCon, double &gGrav, bool &isMembrane);

    LinkedCellContainer(ForceBase &model, std::vector<Particle> &particles, std::vector<double> &dSize,
                        double &cRadius, std::array<char, 6> bCon, double &gGrav, bool &isMembrane);

    LinkedCellContainer(ForceBase &model, std::vector<double> &dSize, double &cRadius, std::array<char, 6> bCon, double &gGrav, bool &isMembrane,
                        int &k, double &r0, double &pullUpF);

    char &getBoundaryCon(int index);

    bool checkBoundary(char b);

    void addParticle(Particle &particle);

    std::vector<Particle> &getParticles();

    std::vector<std::vector<std::vector<std::vector<Particle>>>> &getGrid();

    int getDimension() const;

    std::vector<double> getDomainSize();

    std::size_t size() const;

    void resetF();

    void calculateF();

    void calcF();

    void calcNF();

    //  bool withinCutoff(Particle &p1, Particle &p2) const;

    void calculateX(double delta_t);

    void calculateV(double delta_t);

    void plotParticles(int iteration);

    void initGrid();

    bool checkDistance(Particle &p, std::string border);

    void applyReflecting(Particle &p);

    void applyPeriodic(Particle &p, int x, int y, int z);

    std::array<double, 3> applyNeighbouringForce(Particle &p1, Particle &p2, std::string type);
};

#endif //PSEMOLDYN_GROUPE_LINKEDCELLCONTAINER_H
