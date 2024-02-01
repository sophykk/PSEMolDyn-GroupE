//
// Created by sophy on 30.11.2023.
//

#include "LinkedCellContainer.h"
#include <vector>
#include "Formulas.h"
#include "utils/ArrayUtils.h"
#include "inputOutput/outputWriter/VTKWriter.h"
#include <cmath>
#include <limits>
#include <cstdlib>
#include <spdlog/spdlog.h>

//check
LinkedCellContainer::LinkedCellContainer(ForceBase &model, std::vector<double> &dSize, double &cRadius,
                                         std::array<char, 6> bCon, double &gGrav, int &task, bool &isMembrane) :
        ParticleContainerBase(model), domainSize(dSize), cutoffRadius(cRadius), boundaryCon(bCon), gGrav(gGrav),
        task(task), isMembrane(isMembrane) {
    int cellDimensionX;
    int cellDimensionY;
    int cellDimensionZ;
    if(task == 1){
        gGravPos = 2;
    } else {
        gGravPos = 1;
    }
    std::fmod(domainSize[0], cutoffRadius) == 0 ? cellDimensionX = domainSize[0] / cutoffRadius
                                                : cellDimensionX = std::ceil(domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? cellDimensionY = domainSize[1] / cutoffRadius
                                                : cellDimensionY = std::ceil(domainSize[1] / cutoffRadius);
    std::fmod(domainSize[2], cutoffRadius) == 0 ? cellDimensionZ = domainSize[2] / cutoffRadius
                                                : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius);
    grid.resize(cellDimensionX,
                std::vector<std::vector<std::vector<Particle>>>(
                        cellDimensionY,
                        std::vector<std::vector<Particle>>(
                                cellDimensionZ, std::vector<Particle>())));
    initGrid();
}

//check, if z grid = 1 => 2D
LinkedCellContainer::LinkedCellContainer(ForceBase &model, std::vector<Particle> &particles, std::vector<double> &dSize,
                                         double &cRadius, std::array<char, 6> bCon, double &gGrav, int &task, bool &isMembrane) :
        ParticleContainerBase(model), particleList(particles), domainSize(dSize), cutoffRadius(cRadius),
        boundaryCon(bCon), gGrav(gGrav), task(task), isMembrane(isMembrane) {
    int cellDimensionX;
    int cellDimensionY;
    int cellDimensionZ;
    if(task == 1){
        gGravPos = 2;
    } else {
        gGravPos = 1;
    }
    std::fmod(domainSize[0], cutoffRadius) == 0 ? cellDimensionX = domainSize[0] / cutoffRadius
                                                : cellDimensionX = std::ceil(domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? cellDimensionY = domainSize[1] / cutoffRadius
                                                : cellDimensionY = std::ceil(domainSize[1] / cutoffRadius);
    std::fmod(domainSize[2], cutoffRadius) == 0 ? cellDimensionZ = domainSize[2] / cutoffRadius
                                                : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius);
    grid.resize(cellDimensionX,
                std::vector<std::vector<std::vector<Particle>>>(
                        cellDimensionY,
                        std::vector<std::vector<Particle>>(
                                cellDimensionZ, std::vector<Particle>())));
    initGrid();
}

LinkedCellContainer::LinkedCellContainer(ForceBase &model, std::vector<double> &dSize, double &cRadius,
                                         std::array<char, 6> bCon, double &gGrav, int &task, bool &isMembrane, int &k,
                                         double &r0, double &pullUpF) :
        ParticleContainerBase(model), domainSize(dSize), cutoffRadius(cRadius), boundaryCon(bCon), gGrav(gGrav), task(task),
        isMembrane(isMembrane), k(k), r0(r0), pullUpF(pullUpF) {
    int cellDimensionX;
    int cellDimensionY;
    int cellDimensionZ;
    if(task == 1){
        gGravPos = 2;
    } else {
        gGravPos = 1;
    }
    std::fmod(domainSize[0], cutoffRadius) == 0 ? cellDimensionX = domainSize[0] / cutoffRadius
                                                : cellDimensionX = std::ceil(domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? cellDimensionY = domainSize[1] / cutoffRadius
                                                : cellDimensionY = std::ceil(domainSize[1] / cutoffRadius);
    std::fmod(domainSize[2], cutoffRadius) == 0 ? cellDimensionZ = domainSize[2] / cutoffRadius
                                                : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius);
    grid.resize(cellDimensionX,
                std::vector<std::vector<std::vector<Particle>>>(
                        cellDimensionY,
                        std::vector<std::vector<Particle>>(
                                cellDimensionZ, std::vector<Particle>())));
    initGrid();
}

std::vector<Particle> &LinkedCellContainer::getParticles() {
    return particleList;
}

std::vector<std::vector<std::vector<std::vector<Particle>>>> &LinkedCellContainer::getGrid() {
    return grid;
}

char &LinkedCellContainer::getBoundaryCon(int index) {
    return boundaryCon[index];
}

/**
 * if the z-Dimension vector is more than 1 it is 3D otherwise 2D
 * @return how many dimensions we are working with (2D/3D)
 */
int LinkedCellContainer::getDimension() const {
    return grid[0][0].size() == 1 ? 2 : 3;
}

bool LinkedCellContainer::checkBoundary(char b) {
    for (int i = 0; i < boundaryCon.size(); ++i) {
        if (getBoundaryCon(i) == b) {
            return true;
        }
    }
    return false;
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

void LinkedCellContainer::plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(size());
    for (auto &p: getParticles()) {
        writer.plotParticle(p);
    }
    writer.writeFile(out_name, iteration);
}

bool LinkedCellContainer::checkDistance(Particle &p, std::string border) {
    //divided by 2.0 because of border p____|____p
    auto disCheck = std::pow(2, (1 / 6)) * p.getSigma() / 2.0;

    if (border == "left") {
        if (p.getX()[0] < disCheck && p.getX()[0] >= 0.0) {
            return true;
        }
    }
    if (border == "right") {
        if (domainSize[0] - p.getX()[0] < disCheck && p.getX()[0] <= domainSize[0]) {
            return true;
        }
    }
    if (border == "floor") {
        if (p.getX()[1] < disCheck && p.getX()[1] >= 0.0) {
            return true;
        }
    }
    if (border == "upper") {
        if (domainSize[1] - p.getX()[1] < disCheck && p.getX()[1] <= domainSize[1]) {
            return true;
        }
    }
    if (border == "front") {
        if (p.getX()[2] < disCheck && p.getX()[2] >= 0.0) {
            return true;
        }
    }
    if (border == "back") {
        if (domainSize[2] - p.getX()[2] < disCheck && p.getX()[2] <= domainSize[2]) {
            return true;
        }
    }
    return false;
}

std::vector<double> LinkedCellContainer::getDomainSize() {
    return domainSize;
}

/**
 * @brief Calculate xi(tn+1)
 * formula: xi(tn+1) = xi(tn) + ∆t · vi(tn) + (∆t)^2 * (Fi(tn)/2mi)
 */

void LinkedCellContainer::calculateX(double delta_t) {

    for (auto &p: getParticles()) {
        auto xi_tn1 = Formulas::verletXStep(p.getX(), p.getV(), p.getF(), p.getM(), delta_t);
        //this is another periodic boundary same as in initGrid
        if (checkBoundary('p')) {
            if (boundaryCon[0] == 'p' && p.getX()[0] < 0.0) {
                xi_tn1[0] = p.getX()[0] + domainSize[0];
            }
            // upper
            if (boundaryCon[1] == 'p' && p.getX()[1] > domainSize[1]) {
                xi_tn1[1] = p.getX()[1] - domainSize[1];
            }
            // right
            if (boundaryCon[2] == 'p' && p.getX()[0] > domainSize[0]) {
                xi_tn1[0] = p.getX()[0] - domainSize[0];
            }
            // floor
            if (boundaryCon[3] == 'p' && p.getX()[1] < 0.0) {
                xi_tn1[1] = p.getX()[1] + domainSize[1];
            }
            // front
            if (boundaryCon[4] == 'p' && p.getX()[2] < 0.0) {
                xi_tn1[2] = p.getX()[2] + domainSize[2];
            }
            // back
            if (boundaryCon[5] == 'p' && p.getX()[2] > domainSize[0]) {
                xi_tn1[2] = p.getX()[2] - domainSize[2];
            }
        }
        if (!p.getIsWall()) {
            p.setX(xi_tn1); // Update the position
        }

    }
}

/**
 * @brief Calculate the velocity at time tn+1
 * formula: vi(tn+1) = vi(tn) + ∆t * ((Fi(tn) + Fi(tn+1))/ 2mi)
 */
void LinkedCellContainer::calculateV(double delta_t) {

    for (auto &p: getParticles()) {
        auto vi_tn1 = Formulas::verletVStep(p.getV(), p.getOldF(), p.getF(), p.getM(), delta_t);
        p.setV(vi_tn1);
    }
}

/**
 * @brief This methods puts all the particles in a cell in the grid according to their positions
 * if cutoffRadius = 10, and dSize = {180, 90, 1}
 * then first cell from 0-10, next 11-20,
 */
void LinkedCellContainer::initGrid() {

    //periodic con
    if (checkBoundary('p')) {
        for (auto &p: particleList) {
            //left
            if (boundaryCon[0] == 'p' && p.getX()[0] < 0.0) {
                p.setX({p.getX()[0] + domainSize[0], p.getX()[1], p.getX()[2]});
            }
            // upper
            if (boundaryCon[1] == 'p' && p.getX()[1] > domainSize[1]) {
                p.setX({p.getX()[0], p.getX()[1] - domainSize[1], p.getX()[2]});
            }
            // right
            if (boundaryCon[2] == 'p' && p.getX()[0] > domainSize[0]) {
                p.setX({p.getX()[0] - domainSize[0], p.getX()[1], p.getX()[2]});
            }
            // floor
            if (boundaryCon[3] == 'p' && p.getX()[1] < 0.0) {
                p.setX({p.getX()[0], p.getX()[1] + domainSize[1], p.getX()[2]});
                //spdlog::info("this is in initGrid periodic floor");
            }
            // front
            if (boundaryCon[4] == 'p' && p.getX()[2] < 0.0) {
                p.setX({p.getX()[0], p.getX()[1], p.getX()[2] + domainSize[2]});
            }
            // back
            if (boundaryCon[5] == 'p' && p.getX()[2] > domainSize[2]) {
                p.setX({p.getX()[0], p.getX()[1], p.getX()[2] - domainSize[2]});
            }
        }
    }
    double minDouble = std::numeric_limits<double>::min();
    double lowerZ = 0.0;
    double upperZ = cutoffRadius;
    //x0,y --> 0,0 --> 0,1
    for (int z = 0; z < std::ceil(domainSize[2] / cutoffRadius); ++z) {
        double lowerY = 0.0;
        double upperY = cutoffRadius;
        for (int y = 0; y < std::ceil(domainSize[1] / cutoffRadius); ++y) {
            double lowerX = 0.0;
            double upperX = cutoffRadius;
            //x,y0 --> 0,0 --> 1,0 --> 2,0 --> 3,0 then at the end x,y,0 -> x,y,1
            for (int x = 0; x < std::ceil(domainSize[0] / cutoffRadius); ++x) {
                grid[x][y][z].clear();
                for (auto p = getParticles().begin(); p < getParticles().end(); p++) {
                    if (p->getX()[0] <= upperX && p->getX()[0] >= lowerX &&
                        p->getX()[1] <= upperY && p->getX()[1] >= lowerY &&
                        p->getX()[2] <= upperZ && p->getX()[2] >= lowerZ) {
                        grid[x][y][z].push_back(*p);
                    }
                }
                /* before without pointers
                 * for (auto &p1: getParticles()) {
                    if (p1.getX()[0] <= upperX && p1.getX()[0] >= lowerX &&
                        p1.getX()[1] <= upperY && p1.getX()[1] >= lowerY &&
                        p1.getX()[2] <= upperZ && p1.getX()[2] >= lowerZ) {
                        grid[x][y][z].push_back(p1);
                }*/
                lowerX = upperX + minDouble;
                if (upperX + cutoffRadius > domainSize[0]) {
                    upperX += std::fmod(domainSize[0], cutoffRadius);
                } else {
                    upperX += cutoffRadius;
                }
            }
            lowerY = upperY + minDouble;
            if (upperY + cutoffRadius > domainSize[1]) {
                upperY += std::fmod(domainSize[1], cutoffRadius);
            } else {
                upperY += cutoffRadius;
            }
        }
        lowerZ = upperZ + minDouble;
        if (upperZ + cutoffRadius > domainSize[2]) {
            upperZ += std::fmod(domainSize[2], cutoffRadius);
        } else {
            upperZ += cutoffRadius;
        }
    }
    //spdlog::info("particleListSize in initGrid after {}", particleList.size());

}


void LinkedCellContainer::calculateF() {
    if (isMembrane) {
        calcNF();
    } else {
        calcF();
    }
}


/**
 * @brief calculates force between all particles in cell with itself and its neighbours
 * start in the left corner cell, only has  3(4) neighbours
 * always max 4 calculations for each cell because
 * e.g u start bottom left and calc above, right and left diagonal
 * now for right cell you don't need to calc for left neighbor anymore -> Newton 3. Law,
 * now just calc for above, right and right diagonal and left diagonal if exists
 * and same for one dimension behind if exists
 */
//if not membrane, then use
void LinkedCellContainer::calcF() {
    std::vector<Particle> nParticleList;
    initGrid();
    //down to up
    for (int z = 0; z < getGrid()[0][0].size(); ++z) {
        for (int y = 0; y < getGrid()[0].size(); ++y) {
            //left to right
            for (int x = 0; x < getGrid().size(); ++x) {
                for (auto p1 = grid[x][y][z].begin(); p1 < grid[x][y][z].end(); p1++) {
                    for (auto p2 = p1 + 1; p2 < grid[x][y][z].end(); p2++) {
                        auto force = forceModel.calculateForce(*p1, *p2);
                        if (!p1->getIsWall()) {
                            force[gGravPos] += (p1->getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            if (force[2] < -100) {
                                force[2] = -100;
                            }
                            p1->addF(force);
                        }
                        if (!p2->getIsWall()) {
                            force[gGravPos] += (p2->getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            if (force[2] < -100) {
                                force[2] = -100;
                            }
                            p2->addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                        }
                    }
                    // check done
                    if (checkBoundary('r')) {
                        applyReflecting(*p1);
                    }
                    // check done
                    if (checkBoundary('p')) {
                        applyPeriodic(*p1, x, y, z);
                    }
                }
                //if above neighbour exists (y+1)
                if (y + 1 > 0 && y + 1 < getGrid()[x].size()) {
                    //calculate force between this and above cell for all particles
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x][y + 1][z]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
                //if above back neighbour exists (y+1 & z+1)
                if (y + 1 > 0 && y + 1 < getGrid()[x].size() && z + 1 > 0 && z + 1 < getGrid()[x][y].size()) {
                    //calculate force between this and above back cell for all particles
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x][y + 1][z + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
                //if right hand neighbour exists (x+1)
                if (x + 1 > 0 && x + 1 < getGrid().size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x + 1][y][z]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
                //if right back neighbour exists (x+1 & z+1)
                if (x + 1 > 0 && x + 1 < getGrid().size() && z + 1 > 0 && z + 1 < getGrid()[x][y].size()) {
                    //calculate force between this and right back cell for all particles
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x + 1][y][z + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
                //if right upper diagonal neighbour exists (x+1 & y+1)
                if (x + 1 > 0 && x + 1 < getGrid().size() && y + 1 > 0 && y + 1 < getGrid()[x].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x + 1][y + 1][z]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
                //if right upper back diagonal neighbour exists (x+1 & y+1 & z+1)
                if (x + 1 > 0 && x + 1 < getGrid().size() && y + 1 > 0 && y + 1 < getGrid()[x].size()
                    && z + 1 > 0 && z + 1 < getGrid()[x][y].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x + 1][y + 1][z + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
                //if left upper diagonal neighbour exists (x-1 & y+1)
                if (x - 1 >= 0 && x - 1 < getGrid().size() && y + 1 > 0 && y + 1 < getGrid()[x].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x - 1][y + 1][z]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
                //if left upper back diagonal neighbour exists (x-1 & y+1 & z+1)
                if (x - 1 >= 0 && x - 1 < getGrid().size() && y + 1 > 0 && y + 1 < getGrid()[x].size() && z + 1 > 0 &&
                    z + 1 < getGrid()[x][y].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x - 1][y + 1][z + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
                //if back neighbour exists (z+1)
                if (z + 1 > 0 && z + 1 < getGrid()[x][y].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x][y][z + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            if (!p1.getIsWall()) {
                                force[gGravPos] += (p1.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p1.addF(force);
                            }
                            if (!p2.getIsWall()) {
                                force[gGravPos] += (p2.getM() * gGrav);
                                //capping
                                if (force[0] > 100) {
                                    force[0] = 100;
                                }
                                if (force[0] < -100) {
                                    force[0] = -100;
                                }
                                if (force[1] > 100) {
                                    force[1] = 100;
                                }
                                if (force[1] < -100) {
                                    force[1] = -100;
                                }
                                if (force[2] > 100) {
                                    force[2] = 100;
                                }
                                if (force[2] < -100) {
                                    force[2] = -100;
                                }
                                p2.addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
                            }
                        }
                    }
                }
            }
        }
    }
    for (auto &x: getGrid()) {
        for (auto &y: x) {
            for (auto &z: y) {
                for (auto &p: z) {
                    nParticleList.push_back(p);
                }
            }
        }
    }
    particleList = nParticleList;

}

//if membrane -> use neighboring forces
void LinkedCellContainer::calcNF() {
    initGrid();

    for (auto p1 = particleList.begin(); p1 < particleList.end(); p1++) {
        if (checkBoundary('r')) {
            applyReflecting(*p1);
        }
        //check done
        if (checkBoundary('p')) {
            applyPeriodic(*p1, ceil(p1->getX()[0] / cutoffRadius), ceil(p1->getX()[1] / cutoffRadius),
                          ceil(p1->getX()[2] / cutoffRadius));
        }
    }

    for (auto p1 = particleList.begin(); p1 < particleList.end() - 1; p1++) {
        for (auto p2 = p1 + 1; p2 < particleList.end(); p2++) {
            double sigma;
            if (p1->getType() != p2->getType()) {
                sigma = (p1->getSigma() + p2->getSigma()) / 2;
            } else {
                sigma = p1->getSigma();
            }
            auto dis = std::pow(2, (1 / 6)) * sigma;
            std::array<double, 3> force{0.0, 0.0, 0.0};

            if ((p1->getXIndex() == p2->getXIndex() && abs(p2->getYIndex() - p1->getYIndex()) == 1)
                || (p1->getYIndex() == p2->getYIndex() && abs(p2->getXIndex() - p1->getXIndex()) == 1)) {
                force = applyNeighbouringForce(*p1, *p2, "direct");

            } else if (abs(p2->getXIndex() - p1->getXIndex()) == 1 && abs(p2->getYIndex() - p1->getYIndex()) == 1) {
                force = applyNeighbouringForce(*p1, *p2, "diagonal");

            } else if (abs(p1->getX()[0] - p2->getX()[0]) < dis &&
                       abs(p1->getX()[1] - p2->getX()[1]) < dis &&
                       abs(p1->getX()[2] - p2->getX()[2]) < dis) {
                force = forceModel.calculateForce(*p1, *p2);
            }
            //The particles with x/y-indices (17/24), (17/25), (18/24) and (18/25) are pulled
            //by a constant force FZ-UP “upwards“ along the z-axis:
            /*if ((p1->getXIndex() == 17 && p1->getXIndex() == 24) ||
                (p1->getXIndex() == 17 && p1->getXIndex() == 25) ||
                (p1->getXIndex() == 18 && p1->getXIndex() == 24) ||
                (p1->getXIndex() == 18 && p1->getXIndex() == 25)) {*/
            if (p1->getType() == 5) {
                force[2] += pullUpF;
            }
            if (!p1->getIsWall()) {
                force[gGravPos] += (p1->getM() * gGrav);
                //capping
                if (force[0] > 100) {
                    force[0] = 100;
                }
                if (force[0] < -100) {
                    force[0] = -100;
                }
                if (force[1] > 100) {
                    force[1] = 100;
                }
                if (force[1] < -100) {
                    force[1] = -100;
                }
                if (force[2] > 100) {
                    force[2] = 100;
                }
                if (force[2] < -100) {
                    force[2] = -100;
                }
                p1->addF(force);
            }
            if (!p2->getIsWall()) {
                force[gGravPos] += (p2->getM() * gGrav);
                //capping
                if (force[0] > 100) {
                    force[0] = 100;
                }
                if (force[0] < -100) {
                    force[0] = -100;
                }
                if (force[1] > 100) {
                    force[1] = 100;
                }
                if (force[1] < -100) {
                    force[1] = -100;
                }
                if (force[2] > 100) {
                    force[2] = 100;
                }
                if (force[2] < -100) {
                    force[2] = -100;
                }
                p2->addF({-1.0 * force[0], -1.0 * force[1], -1.0 * force[2]});
            }
        }
    }

    /*for (auto &x: getGrid()) {
        for (auto &y: x) {
            for (auto &z: y) {
                for (auto &p: z) {
                    nParticleList.push_back(p);
                }
            }
        }
    }
    particleList = nParticleList;*/
}

/**
 * @brief Apply reflecting condition on particles
 * @param p Particle to check
 * bounces back borderline particle from border through
 * creating new halo particle which applies force on it
 */
void LinkedCellContainer::applyReflecting(Particle &p) {
    auto dis = std::pow(2, (1 / 6)) * p.getSigma();

    Particle halo(p);
    //near left border
    if (boundaryCon[0] == 'r') {
        if (checkDistance(p, "left")) {
            halo.setX({(-1) * p.getX()[0], p.getX()[1], p.getX()[2]});
            auto force = forceModel.calculateForce(p, halo);
            if (!p.getIsWall()) {
                force[gGravPos] += (p.getM() * gGrav);
                //capping
                if (force[0] > 100) {
                    force[0] = 100;
                }
                if (force[0] < -100) {
                    force[0] = -100;
                }
                if (force[1] > 100) {
                    force[1] = 100;
                }
                if (force[1] < -100) {
                    force[1] = -100;
                }
                if (force[2] > 100) {
                    force[2] = 100;
                }
                if (force[2] < -100) {
                    force[2] = -100;
                }
                p.addF(force);
            }
        }
    }
    //near upper border
    if (boundaryCon[1] == 'r') {
        if (checkDistance(p, "upper")) {
            halo.setX({p.getX()[0], domainSize[1] + (domainSize[1] - p.getX()[1]), p.getX()[2]});
            auto force = forceModel.calculateForce(p, halo);
            if (!p.getIsWall()) {
                force[gGravPos] += (p.getM() * gGrav);
                //capping
                if (force[0] > 100) {
                    force[0] = 100;
                }
                if (force[0] < -100) {
                    force[0] = -100;
                }
                if (force[1] > 100) {
                    force[1] = 100;
                }
                if (force[1] < -100) {
                    force[1] = -100;
                }
                if (force[2] > 100) {
                    force[2] = 100;
                }
                if (force[2] < -100) {
                    force[2] = -100;
                }
                p.addF(force);
            }
        }
    }
    //near right border
    if (boundaryCon[2] == 'r') {
        if (checkDistance(p, "right")) {
            halo.setX({domainSize[0] + (domainSize[0] - p.getX()[0]), p.getX()[1], p.getX()[2]});
            auto force = forceModel.calculateForce(p, halo);
            if (!p.getIsWall()) {
                force[gGravPos] += (p.getM() * gGrav);
                //capping
                if (force[0] > 100) {
                    force[0] = 100;
                }
                if (force[0] < -100) {
                    force[0] = -100;
                }
                if (force[1] > 100) {
                    force[1] = 100;
                }
                if (force[1] < -100) {
                    force[1] = -100;
                }
                if (force[2] > 100) {
                    force[2] = 100;
                }
                if (force[2] < -100) {
                    force[2] = -100;
                }
                p.addF(force);
            }
        }
    }
    //near floor border
    if (boundaryCon[3] == 'r') {
        if (checkDistance(p, "floor")) {
            halo.setX({p.getX()[0], (-1) * p.getX()[1], p.getX()[2]});
            auto force = forceModel.calculateForce(p, halo);
            if (!p.getIsWall()) {
                //capping
                force[gGravPos] += (p.getM() * gGrav);
                if (force[0] > 100) {
                    force[0] = 100;
                }
                if (force[0] < -100) {
                    force[0] = -100;
                }
                if (force[1] > 100) {
                    force[1] = 100;
                }
                if (force[1] < -100) {
                    force[1] = -100;
                }
                if (force[2] > 100) {
                    force[2] = 100;
                }
                if (force[2] < -100) {
                    force[2] = -100;
                }
                p.addF(force);
            }
        }
    }
    if (getDimension() != 2) {
        //near front border
        if (boundaryCon[4] == 'r') {
            if (checkDistance(p, "front")) {
                halo.setX({p.getX()[0], p.getX()[1], p.getX()[2] - dis});
                auto force = forceModel.calculateForce(p, halo);
                if (!p.getIsWall()) {
                    force[gGravPos] += (p.getM() * gGrav);
                    //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                    if (force[2] < -100) {
                        force[2] = -100;
                    }
                    p.addF(force);
                }
            }
        }
        //near back border
        if (boundaryCon[5] == 'r') {
            if (checkDistance(p, "back")) {
                halo.setX({p.getX()[0], p.getX()[1], p.getX()[2] + dis});
                auto force = forceModel.calculateForce(p, halo);
                if (!p.getIsWall()) {
                    force[gGravPos] += (p.getM() * gGrav);
                    //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                    if (force[2] < -100) {
                        force[2] = -100;
                    }
                    p.addF(force);
                }
            }
        }
    }
}

/**
 * @brief just halos, no inserting on the other side here
 * and creating halos if near border -> should affect particles there
 * @param p particle
 * @param x cell location x
 * @param y cell location y
 */
//check done
void LinkedCellContainer::applyPeriodic(Particle &p, int x, int y, int z) {
    Particle halo(p);

    // left
    if (boundaryCon[0] == 'p') {
        if (checkDistance(p, "left")) {
            // halo on right
            halo.setX({p.getX()[0] + domainSize[0], p.getX()[1], p.getX()[2]});
            for (auto &p1: getGrid()[getGrid().size() - 1][y][z]) {
                if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[gGravPos] += (p1.getM() * gGrav);
                    //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                    p1.addF(force);
                }
            }
            // left front
            if (boundaryCon[4] == 'p' && checkDistance(p, "front")) {
                // halo on right back
                halo.setX({p.getX()[0] + domainSize[0], p.getX()[1], p.getX()[2] + domainSize[2]});
                for (auto &p1: getGrid()[getGrid().size() - 1][y][getGrid()[x][y].size() - 1]) {
                    if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                        abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                        p1.addF(force);
                    }
                }
            }
            // left back
            if (boundaryCon[5] == 'p' && checkDistance(p, "back")) {
                // halo on right front
                halo.setX({p.getX()[0] + domainSize[0], p.getX()[1], p.getX()[2] - domainSize[2]});
                for (auto &p1: getGrid()[getGrid().size() - 1][y][0]) {
                    if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                        abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                        p1.addF(force);
                    }
                }
            }
        }
    }

    // right
    if (boundaryCon[2] == 'p') {
        if (checkDistance(p, "right")) {
            //halo on left
            halo.setX({p.getX()[0] - domainSize[0], p.getX()[1], p.getX()[2]});
            for (auto &p1: getGrid()[0][y][z]) {
                if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[gGravPos] += (p1.getM() * gGrav);
                    //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                    p1.addF(force);
                }
            }
            // right front
            if (boundaryCon[4] == 'p' && checkDistance(p, "front")) {
                // halo on left back
                halo.setX({p.getX()[0] - domainSize[0], p.getX()[1], p.getX()[2] + domainSize[2]});
                for (auto &p1: getGrid()[0][y][getGrid()[x][y].size() - 1]) {
                    if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                        abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                        p1.addF(force);
                    }
                }
            }
            // right back
            if (boundaryCon[5] == 'p' && checkDistance(p, "back")) {
                // halo on left front
                halo.setX({p.getX()[0] - domainSize[0], p.getX()[1], p.getX()[2] - domainSize[2]});
                for (auto &p1: getGrid()[0][y][0]) {
                    if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                        abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                        p1.addF(force);
                    }
                }
            }
        }
    }

    // upper
    if (boundaryCon[1] == 'p') {
        if (checkDistance(p, "upper")) {
            // halo on floor
            halo.setX({p.getX()[0], p.getX()[1] - domainSize[1], p.getX()[2]});
            for (auto &p1: getGrid()[x][0][z]) {
                if (abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[gGravPos] += (p1.getM() * gGrav);
                    //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                    p1.addF(force);
                }
            }
            // left (upper)
            if (boundaryCon[0] == 'p' && checkDistance(p, "left")) {
                // halo on right floor
                halo.setX({p.getX()[0] + domainSize[0], p.getX()[1] - domainSize[1], p.getX()[2]});
                for (auto &p1: getGrid()[getGrid().size() - 1][0][z]) {
                    if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                        abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                        if (force[0] > 100) {
                            force[0] = 100;
                        }
                        if (force[0] < -100) {
                            force[0] = -100;
                        }
                        if (force[1] > 100) {
                            force[1] = 100;
                        }
                        if (force[1] < -100) {
                            force[1] = -100;
                        }
                        if (force[2] > 100) {
                            force[2] = 100;
                        }
                        p1.addF(force);
                    }
                }
                // (left upper) front
                if (boundaryCon[4] == 'p' && checkDistance(p, "front")) {
                    // halo on right floor back
                    halo.setX({p.getX()[0] + domainSize[0], p.getX()[1] - domainSize[1], p.getX()[2] + domainSize[2]});
                    for (auto &p1: getGrid()[getGrid().size() - 1][0][getGrid()[x][y].size() - 1]) {
                        if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                            abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                            abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                            auto force = forceModel.calculateForce(p1, halo);
                            force[gGravPos] += (p1.getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            p1.addF(force);
                        }
                    }
                }
                // (left upper) back
                if (boundaryCon[5] == 'p' && checkDistance(p, "back")) {
                    //halo on right floor front
                    halo.setX({p.getX()[0] + domainSize[0], p.getX()[1] - domainSize[1], p.getX()[2] - domainSize[2]});
                    for (auto &p1: getGrid()[getGrid().size() - 1][0][0]) {
                        if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                            abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                            abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                            auto force = forceModel.calculateForce(p1, halo);
                            force[gGravPos] += (p1.getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            p1.addF(force);
                        }
                    }
                }
            }
            // right (upper)
            if (boundaryCon[2] == 'p' && checkDistance(p, "right")) {
                // halo on left floor
                halo.setX({p.getX()[0] - domainSize[0], p.getX()[1] - domainSize[1], p.getX()[2]});
                for (auto &p1: getGrid()[0][0][z]) {
                    if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                        abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                        if (force[0] > 100) {
                            force[0] = 100;
                        }
                        if (force[0] < -100) {
                            force[0] = -100;
                        }
                        if (force[1] > 100) {
                            force[1] = 100;
                        }
                        if (force[1] < -100) {
                            force[1] = -100;
                        }
                        if (force[2] > 100) {
                            force[2] = 100;
                        }
                        p1.addF(force);
                    }
                }
                // (right upper) front
                if (boundaryCon[4] == 'p' && checkDistance(p, "front")) {
                    // halo on left floor back
                    halo.setX({p.getX()[0] - domainSize[0], p.getX()[1] - domainSize[1], p.getX()[2] + domainSize[2]});
                    for (auto &p1: getGrid()[0][0][getGrid()[x][y].size() - 1]) {
                        if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                            abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                            abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                            auto force = forceModel.calculateForce(p1, halo);
                            force[gGravPos] += (p1.getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            p1.addF(force);
                        }
                    }
                }
                // (right upper) back
                if (boundaryCon[5] == 'p' && checkDistance(p, "back")) {
                    //halo on left floor front
                    halo.setX({p.getX()[0] - domainSize[0], p.getX()[1] - domainSize[1], p.getX()[2] - domainSize[2]});
                    for (auto &p1: getGrid()[0][0][0]) {
                        if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                            abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                            abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                            auto force = forceModel.calculateForce(p1, halo);
                            force[gGravPos] += (p1.getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            p1.addF(force);
                        }
                    }
                }
            }
            // (upper) front
            if (boundaryCon[4] == 'p' && checkDistance(p, "front")) {
                // halo on floor back
                halo.setX({p.getX()[0], p.getX()[1] - domainSize[1], p.getX()[2] + domainSize[2]});
                for (auto &p1: getGrid()[x][0][getGrid()[x][y].size() - 1]) {
                    if (abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                        abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                        if (force[0] > 100) {
                            force[0] = 100;
                        }
                        if (force[0] < -100) {
                            force[0] = -100;
                        }
                        if (force[1] > 100) {
                            force[1] = 100;
                        }
                        if (force[1] < -100) {
                            force[1] = -100;
                        }
                        if (force[2] > 100) {
                            force[2] = 100;
                        }
                        p1.addF(force);
                    }
                }
            }
            // (upper) back
            if (boundaryCon[5] == 'p' && checkDistance(p, "back")) {
                // halo on floor front
                halo.setX({p.getX()[0], p.getX()[1] - domainSize[1], p.getX()[2] - domainSize[2]});
                for (auto &p1: getGrid()[x][0][0]) {
                    if (abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                        abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                        if (force[0] > 100) {
                            force[0] = 100;
                        }
                        if (force[0] < -100) {
                            force[0] = -100;
                        }
                        if (force[1] > 100) {
                            force[1] = 100;
                        }
                        if (force[1] < -100) {
                            force[1] = -100;
                        }
                        if (force[2] > 100) {
                            force[2] = 100;
                        }
                        p1.addF(force);
                    }
                }
            }
        }
    }

    // floor
    if (boundaryCon[3] == 'p') {
        if (checkDistance(p, "floor")) {
            // halo on upper
            halo.setX({p.getX()[0], p.getX()[1] + domainSize[1], p.getX()[2]});
            for (auto &p1: getGrid()[x][getGrid()[0].size() - 1][z]) {
                if (abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[gGravPos] += (p1.getM() * gGrav);
                    //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                    p1.addF(force);
                }
            }
            // left (floor)
            if (boundaryCon[0] == 'p' && checkDistance(p, "left")) {
                // halo on right upper
                halo.setX({p.getX()[0] + domainSize[0], p.getX()[1] + domainSize[1], p.getX()[2]});
                for (auto &p1: getGrid()[getGrid().size() - 1][getGrid()[x].size() - 1][z]) {
                    if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                        abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                        if (force[0] > 100) {
                            force[0] = 100;
                        }
                        if (force[0] < -100) {
                            force[0] = -100;
                        }
                        if (force[1] > 100) {
                            force[1] = 100;
                        }
                        if (force[1] < -100) {
                            force[1] = -100;
                        }
                        if (force[2] > 100) {
                            force[2] = 100;
                        }
                        p1.addF(force);
                    }
                }
                // (left floor) front
                if (boundaryCon[4] == 'p' && checkDistance(p, "front")) {
                    // halo on right upper back
                    halo.setX({p.getX()[0] + domainSize[0], p.getX()[1] + domainSize[1], p.getX()[2] + domainSize[2]});
                    for (auto &p1: getGrid()[getGrid().size() - 1][getGrid()[x].size() - 1][getGrid()[x][y].size() -
                                                                                            1]) {
                        if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                            abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                            abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                            auto force = forceModel.calculateForce(p1, halo);
                            force[gGravPos] += (p1.getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            p1.addF(force);
                        }
                    }
                }
                // (left floor) back
                if (boundaryCon[5] == 'p' && checkDistance(p, "back")) {
                    //halo on right upper front
                    halo.setX({p.getX()[0] + domainSize[0], p.getX()[1] + domainSize[1], p.getX()[2] - domainSize[2]});
                    for (auto &p1: getGrid()[getGrid().size() - 1][getGrid()[x].size() - 1][0]) {
                        if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                            abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                            abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                            auto force = forceModel.calculateForce(p1, halo);
                            force[gGravPos] += (p1.getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            p1.addF(force);
                        }
                    }
                }
            }
            // right (floor)
            if (boundaryCon[2] == 'p' && checkDistance(p, "right")) {
                // halo on left upper
                halo.setX({p.getX()[0] - domainSize[0], p.getX()[1] + domainSize[1], p.getX()[2]});
                for (auto &p1: getGrid()[0][getGrid()[x].size() - 1][z]) {
                    if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                        abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                        if (force[0] > 100) {
                            force[0] = 100;
                        }
                        if (force[0] < -100) {
                            force[0] = -100;
                        }
                        if (force[1] > 100) {
                            force[1] = 100;
                        }
                        if (force[1] < -100) {
                            force[1] = -100;
                        }
                        if (force[2] > 100) {
                            force[2] = 100;
                        }
                        p1.addF(force);
                    }
                }
                // (right floor) front
                if (boundaryCon[4] == 'p' && checkDistance(p, "front")) {
                    // halo on left upper back
                    halo.setX({p.getX()[0] - domainSize[0], p.getX()[1] + domainSize[1], p.getX()[2] + domainSize[2]});
                    for (auto &p1: getGrid()[0][getGrid()[x].size() - 1][getGrid()[x][y].size() - 1]) {
                        if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                            abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                            abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                            auto force = forceModel.calculateForce(p1, halo);
                            force[gGravPos] += (p1.getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            p1.addF(force);
                        }
                    }
                }
                // (right floor) back
                if (boundaryCon[5] == 'p' && checkDistance(p, "back")) {
                    //halo on left upper front
                    halo.setX({p.getX()[0] - domainSize[0], p.getX()[1] + domainSize[1], p.getX()[2] - domainSize[2]});
                    for (auto &p1: getGrid()[0][getGrid()[x].size() - 1][0]) {
                        if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius &&
                            abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                            abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                            auto force = forceModel.calculateForce(p1, halo);
                            force[gGravPos] += (p1.getM() * gGrav);
                            //capping
                            if (force[0] > 100) {
                                force[0] = 100;
                            }
                            if (force[0] < -100) {
                                force[0] = -100;
                            }
                            if (force[1] > 100) {
                                force[1] = 100;
                            }
                            if (force[1] < -100) {
                                force[1] = -100;
                            }
                            if (force[2] > 100) {
                                force[2] = 100;
                            }
                            p1.addF(force);
                        }
                    }
                }
            }
            // (floor) front
            if (boundaryCon[4] == 'p' && checkDistance(p, "front")) {
                // halo on upper back
                halo.setX({p.getX()[0], p.getX()[1] + domainSize[1], p.getX()[2] + domainSize[2]});
                for (auto &p1: getGrid()[x][getGrid()[x].size() - 1][getGrid()[x][y].size() - 1]) {
                    if (abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                        abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                        if (force[0] > 100) {
                            force[0] = 100;
                        }
                        if (force[0] < -100) {
                            force[0] = -100;
                        }
                        if (force[1] > 100) {
                            force[1] = 100;
                        }
                        if (force[1] < -100) {
                            force[1] = -100;
                        }
                        if (force[2] > 100) {
                            force[2] = 100;
                        }
                        p1.addF(force);
                    }
                }
            }
            // (floor) back
            if (boundaryCon[5] == 'p' && checkDistance(p, "back")) {
                // halo on upper front
                halo.setX({p.getX()[0], p.getX()[1] + domainSize[1], p.getX()[2] - domainSize[2]});
                for (auto &p1: getGrid()[x][getGrid()[x].size() - 1][0]) {
                    if (abs(p1.getX()[1] - halo.getX()[1]) <= cutoffRadius &&
                        abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                        auto force = forceModel.calculateForce(p1, halo);
                        force[gGravPos] += (p1.getM() * gGrav);
                        //capping
                        if (force[0] > 100) {
                            force[0] = 100;
                        }
                        if (force[0] < -100) {
                            force[0] = -100;
                        }
                        if (force[1] > 100) {
                            force[1] = 100;
                        }
                        if (force[1] < -100) {
                            force[1] = -100;
                        }
                        if (force[2] > 100) {
                            force[2] = 100;
                        }
                        p1.addF(force);
                    }
                }
            }
        }
    }

    // front
    if (boundaryCon[4] == 'p') {
        if (checkDistance(p, "front")) {
            // halo on back
            halo.setX({p.getX()[0], p.getX()[1], p.getX()[2] + domainSize[2]});
            for (auto &p1: getGrid()[x][y][getGrid()[x][y].size() - 1]) {
                if (abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[gGravPos] += (p1.getM() * gGrav);
                    //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                    p1.addF(force);
                }
            }
        }
    }

    // back
    if (boundaryCon[5] == 'p') {
        if (checkDistance(p, "back")) {
            // halo on front
            halo.setX({p.getX()[0], p.getX()[1], p.getX()[2] - domainSize[2]});
            for (auto &p1: getGrid()[x][y][0]) {
                if (abs(p1.getX()[2] - halo.getX()[2]) <= cutoffRadius) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[gGravPos] += (p1.getM() * gGrav);
                    //capping
                    if (force[0] > 100) {
                        force[0] = 100;
                    }
                    if (force[0] < -100) {
                        force[0] = -100;
                    }
                    if (force[1] > 100) {
                        force[1] = 100;
                    }
                    if (force[1] < -100) {
                        force[1] = -100;
                    }
                    if (force[2] > 100) {
                        force[2] = 100;
                    }
                    p1.addF(force);
                }
            }
        }
    }
}


std::array<double, 3> LinkedCellContainer::applyNeighbouringForce(Particle &p1, Particle &p2, std::string type) {
    auto p = p2.getX() - p1.getX();
    std::array<double, 3> x = {p[0] / Formulas::secondNorm(p1.getX() - p2.getX()),
                               p[1] / Formulas::secondNorm(p1.getX() - p2.getX()),
                               p[2] / Formulas::secondNorm(p1.getX() - p2.getX())};

    // direct neighbors
    if (type == "direct") {
        return k * ((Formulas::secondNorm(p1.getX() - p2.getX())) - r0) * x;
        // diagonal neighbors
    } else if (type == "diagonal") {
        return k * (Formulas::secondNorm(p1.getX() - p2.getX()) - sqrt(2 * r0)) * x;
    }
    return std::array<double, 3>();
}

/**
 * @brief if cutOffRadius is not cellSize, check distance then
 * @code
 * bool LinkedCellContainer::withinCutoff(Particle &p1, Particle &p2) const {
 *     std::array<double, 3> distance{};
 *     if (p1.getX() < p2.getX()) {
 *         distance = p2.getX() - p1.getX();
 *     } else {
 *         distance = p1.getX() - p2.getX();
 *     }
 *     if (distance[0] > cutoffRadius || distance[1] > cutoffRadius || distance[1] > cutoffRadius) {
 *         return false;
 *     }
 *     return true;
 * }
 * @endcode
 */


