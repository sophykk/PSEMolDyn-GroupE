//
// Created by sophy on 30.11.2023.
//

#include "LinkedCellContainer2.h"
#include <list>
#include <vector>
#include "Formulas.h"
#include "utils/ArrayUtils.h"
#include "inputOutput/outputWriter/VTKWriter.h"
#include <spdlog/spdlog.h>


LinkedCellContainer2::LinkedCellContainer2(ForceBase &model, std::vector<int> &dSize, double &cRadius, int &cSize) :
        ParticleContainerBase(model), domainSize(dSize), cutoffRadius(cRadius), cellSize(cSize) {
    int x = domainSize[0] / cellSize;
    int y = domainSize[1] / cellSize;
    grid.resize(x, std::vector<std::vector<Particle>>(y, std::vector<Particle>(1)));
    initGrid();
}

LinkedCellContainer2::LinkedCellContainer2(ForceBase &model, std::vector<Particle> &particles, std::vector<int> &dSize,
                                           double &cRadius, int &cSize) :
        ParticleContainerBase(model), particleList(particles), domainSize(dSize), cutoffRadius(cRadius),
        cellSize(cSize) {

    int x = domainSize[0] / cellSize;
    int y = domainSize[1] / cellSize;

    grid.resize(x, std::vector<std::vector<Particle>>(y, std::vector<Particle>(1)));
    initGrid();
}

std::vector<Particle> &LinkedCellContainer2::getParticles() {
    return particleList;
}

std::vector<std::vector<std::vector<Particle>>> &LinkedCellContainer2::getGrid() {
    return grid;
}

void LinkedCellContainer2::addParticle(Particle &particle) {
    particleList.push_back(particle);
}

std::size_t LinkedCellContainer2::size() const {
    return particleList.size();
}

void LinkedCellContainer2::resetF() {
    for (auto &p: getParticles()) {
        p.resetF();
    }
}

/*
 * start in the left corner cell, only has  3 neighbours
 * always max 3 calcs for each cell because
 * e.g u start bottom left and calc above, right and right diagonal
 * now for right cell you don't need to calc for left neighbor anymore,
 * because it was just done before in other direction
 * now just calc for above, right and right diagonal if exists
 */
void LinkedCellContainer2::calculateF() {
    for (int i = 0; i < getGrid().size(); ++i) {
        for (int j = 0; j < getGrid()[i].size(); ++j) {
            //if above neigbour exists
            if (j + 1 < getGrid()[i].size()) {
                //calculate force between this and above cell for all particles
                for (auto p1 = getGrid()[i][j].begin(); p1 < getGrid()[i][j].end(); p1++) {
                    for (auto p2 = getGrid()[i][j + 1].begin(); p2 < getGrid()[i][j + 1].end(); p2++) {
                        if (withinCutoff(*p1, *p2)) {
                            auto force = forceModel.calculateForce(*p1, *p2);
                            p1->addF(force);
                            p2->addF(-1.0 * force);
                        }
                    }
                }
            }
            //if right hand neighbour exists
            if (i + 1 < getGrid().size()) {
                for (auto p1 = getGrid()[i][j].begin(); p1 < getGrid()[i][j].end(); p1++) {
                    for (auto p2 = getGrid()[i + 1][j].begin(); p2 < getGrid()[i + 1][j].end(); p2++) {
                        if (withinCutoff(*p1, *p2)) {
                            auto force = forceModel.calculateForce(*p1, *p2);
                            p1->addF(force);
                            p2->addF(-1.0 * force);
                        }
                    }
                }
            }
            //if right diagonal neighbour exists
            if (i + 1 < getGrid().size() && j + 1 < getGrid()[i].size()) {
                for (auto p1 = getGrid()[i][j].begin(); p1 < getGrid()[i][j].end(); p1++) {
                    for (auto p2 = getGrid()[i + 1][j + 1].begin(); p2 < getGrid()[i + 1][j + 1].end(); p2++) {
                        if (withinCutoff(*p1, *p2)) {
                            auto force = forceModel.calculateForce(*p1, *p2);
                            p1->addF(force);
                            p2->addF(-1.0 * force);
                        }
                    }
                }
            }
        }
    }

/* og calc F
for (auto p1 = getParticles().begin(); p1 < getParticles().end(); p1++) {
    for (auto p2 = p1 + 1; p2 < getParticles().end(); p2++) {
        auto force = forceModel.calculateForce(*p1, *p2);
        p1->addF(force);
        p2->addF(-1.0 * force);
    }
}
*/
}

bool LinkedCellContainer2::withinCutoff(Particle &p1, Particle &p2) const {
    std::array<double, 3> distance{};
    if (p1.getX() < p2.getX()) {
        distance = p1.getX() - p2.getX();
    } else {
        distance = p1.getX() - p2.getX();
    }
    if (distance[0] > cutoffRadius || distance[1] > cutoffRadius || distance[1] > cutoffRadius) {
        return false;
    }
    return true;
}

void LinkedCellContainer2::calculateX(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: xi(tn+1) = xi(tn) + ∆t · vi(tn) + (∆t)^2 * (Fi(tn)/2mi)
        // Calculate xi(tn+1)
        auto xi_tn1 = Formulas::verletXStep(p.getX(), p.getV(), p.getF(), p.getM(), delta_t);
        p.setX(xi_tn1);  // Update the position
    }
}

void LinkedCellContainer2::calculateV(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: vi(tn+1) = vi(tn) + ∆t * ((Fi(tn) + Fi(tn+1))/ 2mi)
        // Calculate the velocity at time tn+1
        auto vi_tn1 = Formulas::verletVStep(p.getV(), p.getOldF(), p.getF(), p.getM(), delta_t);
        p.setV(vi_tn1);
    }
}

void LinkedCellContainer2::plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(size());
    for (auto &p: getParticles()) {
        writer.plotParticle(p);
    }
    writer.writeFile(out_name, iteration);
}

/*
 * if cSize = 10, and dSize = {180, 90, 1}
 * then first cell from 0-10, next 11-20,
 */

void LinkedCellContainer2::initGrid() {

    int lowery = 0;
    int uppery = cellSize;
    //y --> 0,1 --> 1,1
    for (int i = 0; i < domainSize[1] / cellSize; ++i) {
        //2. dimension
        int lowerx = 0;
        int upperx = cellSize;
        //x, --> 0,0 --> 1,0 --> 2,0 --> 3,0
        for (int j = 0; j < domainSize[0] / cellSize; ++j) {
            //only for 2D domains {_, _, 1}, otherwise need other implementation (extra loop)
            //3. dimension
            for (auto p1 = getParticles().begin(); p1 < getParticles().end(); p1++) {
                spdlog::debug("x cellSize now {} - {}", lowerx, upperx);
                spdlog::debug("y cellSize now {} - {}", lowery, uppery);
                if (p1->getX()[0] <= upperx && p1->getX()[0] >= lowerx &&
                    p1->getX()[1] <= uppery && p1->getX()[1] >= lowery) {
                    spdlog::debug("we are in grid {}, {}", j, i);
                    spdlog::info("this is the particle before {}", p1->toString());
                    grid[j][i][0] = *p1;
                    spdlog::info("this is p1 after set to grid {}", grid[j][i][0].toString());
                }
            }
            lowerx = upperx + 1;
            upperx += cellSize;
        }
        lowery = uppery + 1;
        uppery += cellSize;
    }
}