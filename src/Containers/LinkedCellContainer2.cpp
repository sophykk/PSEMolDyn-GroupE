//
// Created by sophy on 30.11.2023.
//

#include "LinkedCellContainer2.h"
#include <vector>
#include "Formulas.h"
#include "utils/ArrayUtils.h"
#include "inputOutput/outputWriter/VTKWriter.h"
#include <cmath>
#include <limits>
#include <spdlog/spdlog.h>


LinkedCellContainer2::LinkedCellContainer2(ForceBase &model, std::vector<double> &dSize, double &cRadius, char bCon) :
        ParticleContainerBase(model), domainSize(dSize), cutoffRadius(cRadius), boundaryCon(bCon) {
    int x;
    int y;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? x = domainSize[0] / cutoffRadius : x = std::ceil(
            domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? y = domainSize[1] / cutoffRadius : y = std::ceil(
            domainSize[1] / cutoffRadius);
    grid.resize(x, std::vector<std::vector<Particle>>(y, std::vector<Particle>()));
    initGrid();
}

LinkedCellContainer2::LinkedCellContainer2(ForceBase &model, std::vector<Particle> &particles,
                                           std::vector<double> &dSize, double &cRadius, char bCon) :
        ParticleContainerBase(model), particleList(particles), domainSize(dSize), cutoffRadius(cRadius),
        boundaryCon(bCon) {
    int x;
    int y;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? x = domainSize[0] / cutoffRadius : x = std::ceil(
            domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? y = domainSize[1] / cutoffRadius : y = std::ceil(
            domainSize[1] / cutoffRadius);
    grid.resize(x, std::vector<std::vector<Particle>>(y, std::vector<Particle>()));
    initGrid();
}

std::vector<Particle> &LinkedCellContainer2::getParticles() {
    return particleList;
}

std::vector<std::vector<std::vector<Particle>>> &LinkedCellContainer2::getGrid() {
    return grid;
}

char &LinkedCellContainer2::getBoundaryCon() {
    return boundaryCon;
}

void LinkedCellContainer2::setBoundaryCon(char boundary) {
    if (boundary == 'o' || boundary == 'r') {
        spdlog::error("This is not a valid boundary condition!");
        exit(-1);
    } else {
        boundaryCon = boundary;
    }
}

void LinkedCellContainer2::addParticle(Particle &particle) {
    particleList.push_back(particle);
}

std::size_t LinkedCellContainer2::size() const {
    return particleList.size();
}


/*
 * if cutoffRadius = 10, and dSize = {180, 90, 1}
 * then first cell from 0-10, next 11-20,
 */
void LinkedCellContainer2::initGrid() {

    double minDouble = std::numeric_limits<double>::min();
    double lowery = 0.0;
    double uppery = cutoffRadius;
    //x0,y --> 0,0 --> 0,1
    for (int y = 0; y < std::ceil(domainSize[1] / cutoffRadius); ++y) {
        //2. dimension
        double lowerx = 0.0;
        double upperx = cutoffRadius;
        //x,y0 --> 0,0 --> 1,0 --> 2,0 --> 3,0
        for (int x = 0; x < std::ceil(domainSize[0] / cutoffRadius); ++x) {
            grid[x][y].clear();
            //3. dimension
            for (auto &p1: getParticles()) {
                if (p1.getX()[0] <= upperx && p1.getX()[0] >= lowerx &&
                    p1.getX()[1] <= uppery && p1.getX()[1] >= lowery) {
                    grid[x][y].push_back(p1);

                }
            }
            lowerx = upperx + minDouble;
            if (upperx + cutoffRadius > domainSize[0]) {
                upperx += std::fmod(domainSize[0], cutoffRadius);
            } else {
                upperx += cutoffRadius;
            }
        }
        lowery = uppery + minDouble;
        if (uppery + cutoffRadius > domainSize[1]) {
            uppery += std::fmod(domainSize[1], cutoffRadius);
        } else {
            uppery += cutoffRadius;
        }
    }
}

void LinkedCellContainer2::resetF() {
    for (auto &p: getParticles()) {
        p.resetF();
    }
}

/*
 * start in the left corner cell, only has  3 neighbours
 * always max 4 calcs for each cell because
 * e.g u start bottom left and calc above, right and right diagonal
 * now for right cell you don't need to calc for left neighbor anymore,
 * because it was just done before in other direction
 * now just calc for above, right and right diagonal and left diagonal if exists
 */
void LinkedCellContainer2::calculateF() {
    std::vector<Particle> nParticleList;
    initGrid();
    //down to up
    for (int y = 0; y < getGrid()[0].size(); ++y) {
        //left to right
        for (int x = 0; x < getGrid().size(); ++x) {
            for (auto p1 = grid[x][y].begin(); p1 < grid[x][y].end(); p1++) {
                for (auto p2 = p1 + 1; p2 < grid[x][y].end(); p2++) {
                    auto force = forceModel.calculateForce(*p1, *p2);
                    p1->addF(force);
                    p2->addF(-1.0 * force);
                }
                if (boundaryCon == 'r'){
                    applyReflecting(*p1);
                }
            }
            //if above neigbour exists
            if (y + 1 >= 0 && y + 1 < getGrid()[x].size()) {
                //calculate force between this and above cell for all particles
                for (auto &p1: getGrid()[x][y]) {
                    for (auto &p2: getGrid()[x][y + 1]) {
                        auto force = forceModel.calculateForce(p1, p2);
                        p1.addF(force);
                        p2.addF(-1.0 * force);

                    }
                }
                //if right hand neighbour exists
                if (x + 1.0 >= 0.0 && x + 1.0 < getGrid().size()) {
                    for (auto &p1: getGrid()[x][y]) {
                        for (auto &p2: getGrid()[x + 1][y]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF(force);
                            p2.addF(-1.0 * force);
                        }
                    }
                }
                //if upper right diagonal neighbour exists
                if (y + 1 >= 0 && y + 1 < getGrid()[x].size() && x + 1 >= 0 && x + 1 < getGrid().size()) {
                    for (auto &p1: getGrid()[x][y]) {
                        for (auto &p2: getGrid()[x + 1][y + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF(force);
                            p2.addF(-1.0 * force);
                        }
                    }

                }
                //if upper left diagonal neighbour exists
                if (y + 1 >= 0 && y + 1 < getGrid()[x].size() && x - 1 >= 0 && x - 1 < getGrid().size()) {
                    for (auto &p1: getGrid()[x][y]) {
                        for (auto &p2: getGrid()[x - 1][y + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF(force);
                            p2.addF(-1.0 * force);
                        }
                    }

                }
            }
        }
    }
    for (auto x: getGrid()) {
        for (auto y: x) {
            for (auto z: y) {
                nParticleList.push_back(z);
            }
        }
    }
    particleList = nParticleList;
}

/*
 * if cutOffRadius is not cellSize
bool LinkedCellContainer2::withinCutoff(Particle &p1, Particle &p2) const {
    std::array<double, 3> distance{};
    if (p1.getX() < p2.getX()) {
        distance = p2.getX() - p1.getX();
    } else {
        distance = p1.getX() - p2.getX();
    }
    if (distance[0] > cutoffRadius || distance[1] > cutoffRadius || distance[1] > cutoffRadius) {
        return false;
    }
    return true;
}
*/

void LinkedCellContainer2::calculateX(double delta_t) {
    for (auto &p: getParticles()) {
        // formula: xi(tn+1) = xi(tn) + ∆t · vi(tn) + (∆t)^2 * (Fi(tn)/2mi)
        // Calculate xi(tn+1)
        auto xi_tn1 = Formulas::verletXStep(p.getX(), p.getV(), p.getF(), p.getM(), delta_t);
        p.setX(xi_tn1);  // Update the position
    }
    initGrid();
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

void LinkedCellContainer2::applyReflecting(Particle &p) {
    //1.0 is sigma (not mentioned because doesn't change result)
    //divided by 2.0 because of border p____|____p
    auto disCheck = std::pow(2, 1.0 / 6) / 2.0;
    Particle halo(p);
    //near left border
    if (p.getX()[0] < disCheck && p.getX()[0] >= 0.0) {
        halo.setX({(-1) * p.getX()[0], p.getX()[1], p.getX()[2]});
        auto force = forceModel.calculateForce(p, halo);
        p.addF(force);
    }
    //near right border
    if (domainSize[0] - p.getX()[0] < disCheck && p.getX()[0] <= domainSize[0]) {
        halo.setX({domainSize[0] + (domainSize[0] - p.getX()[0]), p.getX()[1], p.getX()[2]});
        auto force = forceModel.calculateForce(p, halo);
        p.addF(force);
    }
    //near floor border
    if (p.getX()[1] < disCheck && p.getX()[1] >= 0.0) {
        halo.setX({p.getX()[0], (-1) * p.getX()[1], p.getX()[2]});
        auto force = forceModel.calculateForce(p, halo);
        p.addF(force);
    }
    //near upper border
    if (domainSize[1] - p.getX()[1] < disCheck && p.getX()[1] <= domainSize[1]) {
        halo.setX({p.getX()[0], domainSize[1] + (domainSize[1] - p.getX()[1]), p.getX()[2]});
        auto force = forceModel.calculateForce(p, halo);
        p.addF(force);
    }
}
