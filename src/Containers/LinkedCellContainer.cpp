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
                                         std::array<char, 6> bCon, double &gGrav, bool &isMembrane) :
        ParticleContainerBase(model), domainSize(dSize), cutoffRadius(cRadius), boundaryCon(bCon), gGrav(gGrav), isMembrane(isMembrane) {
    int cellDimensionX;
    int cellDimensionY;
    int cellDimensionZ;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? cellDimensionX = domainSize[0] / cutoffRadius
                                                : cellDimensionX = std::ceil(domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? cellDimensionY = domainSize[1] / cutoffRadius
                                                : cellDimensionY = std::ceil(domainSize[1] / cutoffRadius);
    std::fmod(domainSize[2], cutoffRadius) == 0 ? cellDimensionZ = domainSize[2] / cutoffRadius
                                                : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius);
    //domainSize[2] == 1 ? cellDimensionZ = 1 : (std::fmod(domainSize[2], cutoffRadius) == 0 ?
    //        cellDimensionZ = domainSize[2] / cutoffRadius : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius));
    grid.resize(cellDimensionX,
                std::vector<std::vector<std::vector<Particle>>>(
                        cellDimensionY,
                        std::vector<std::vector<Particle>>(
                                cellDimensionZ, std::vector<Particle>())));
    initGrid();
}

//check, maybe if z grid = 1 => 2D
LinkedCellContainer::LinkedCellContainer(ForceBase &model, std::vector<Particle> &particles, std::vector<double> &dSize,
                                         double &cRadius, std::array<char, 6> bCon, double &gGrav, bool &isMembrane) :
        ParticleContainerBase(model), particleList(particles), domainSize(dSize), cutoffRadius(cRadius),
        boundaryCon(bCon), gGrav(gGrav), isMembrane(isMembrane) {
    int cellDimensionX;
    int cellDimensionY;
    int cellDimensionZ;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? cellDimensionX = domainSize[0] / cutoffRadius
                                                : cellDimensionX = std::ceil(domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? cellDimensionY = domainSize[1] / cutoffRadius
                                                : cellDimensionY = std::ceil(domainSize[1] / cutoffRadius);
    std::fmod(domainSize[2], cutoffRadius) == 0 ? cellDimensionZ = domainSize[2] / cutoffRadius
                                                : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius);
    //domainSize[2] == 1 ? cellDimensionZ = 1 : (std::fmod(domainSize[2], cutoffRadius) == 0 ?
    //        cellDimensionZ = domainSize[2] / cutoffRadius : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius));
    grid.resize(cellDimensionX,
                std::vector<std::vector<std::vector<Particle>>>(
                        cellDimensionY,
                        std::vector<std::vector<Particle>>(
                                cellDimensionZ, std::vector<Particle>())));
    initGrid();
}

LinkedCellContainer::LinkedCellContainer(ForceBase &model, std::vector<double> &dSize, double &cRadius,
                                         std::array<char, 6> bCon, double &gGrav, bool &isMembrane, int &k, double &r0, double &pullUpF) :
        ParticleContainerBase(model), domainSize(dSize), cutoffRadius(cRadius), boundaryCon(bCon), gGrav(gGrav), isMembrane(isMembrane), k(k), r0(r0), pullUpF(pullUpF) {
    int cellDimensionX;
    int cellDimensionY;
    int cellDimensionZ;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? cellDimensionX = domainSize[0] / cutoffRadius
                                                : cellDimensionX = std::ceil(domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? cellDimensionY = domainSize[1] / cutoffRadius
                                                : cellDimensionY = std::ceil(domainSize[1] / cutoffRadius);
    std::fmod(domainSize[2], cutoffRadius) == 0 ? cellDimensionZ = domainSize[2] / cutoffRadius
                                                : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius);
    //domainSize[2] == 1 ? cellDimensionZ = 1 : (std::fmod(domainSize[2], cutoffRadius) == 0 ?
    //        cellDimensionZ = domainSize[2] / cutoffRadius : cellDimensionZ = std::ceil(domainSize[2] / cutoffRadius));
    grid.resize(cellDimensionX,
                std::vector<std::vector<std::vector<Particle>>>(
            cellDimensionY,
                    std::vector<std::vector<Particle>>(
                            cellDimensionZ, std::vector<Particle>())));
    initGrid();
}

//check done
std::vector<Particle> &LinkedCellContainer::getParticles() {
    return particleList;
}

//check done
std::vector<std::vector<std::vector<std::vector<Particle>>>> &LinkedCellContainer::getGrid() {
    return grid;
}

//check done
char &LinkedCellContainer::getBoundaryCon(int index) {
    return boundaryCon[index];
}

//check done
/**
 * if the z-Dimension vector is more than 1 it is 3D otherwise 2D
 * @return how many dimensions we are working with (2D/3D)
 */
int LinkedCellContainer::getDimension() const {
    return grid[0][0].size() == 1 ? 2 : 3;
}

//check done
/**
 * @brief Sets the boundary to outflow / reflecting / periodic
 * @param boundary a char, should be 'o' or 'r' or 'p'
 * @throws error message invalid_argument if boundary is not 'o' or 'r' or 'p'.
 */
void LinkedCellContainer::setBoundaryCon(std::array<char, 6> &nBoundary) {
    for (int i = 0; i < nBoundary.size(); i++) {
        if (nBoundary[i] == 'o' || nBoundary[i] == 'p' || nBoundary[i] == 'r') {
            boundaryCon[i] = nBoundary[i];
        } else {
            spdlog::error("Index {} is not a valid boundary condition!", i);
        }
    }
}

//check done
bool LinkedCellContainer::checkBoundary(char b) {
    for (int i = 0; i < boundaryCon.size(); ++i) {
        if (getBoundaryCon(i) == b) {
            return true;
        }
    }
    return false;
}

//check done
void LinkedCellContainer::addParticle(Particle &particle) {
    particleList.push_back(particle);
}

//check done
std::size_t LinkedCellContainer::size() const {
    return particleList.size();
}


//check done
/**
 * @brief This methods puts all the particles in a cell in the grid according to their positions
 * if cutoffRadius = 10, and dSize = {180, 90, 1}
 * then first cell from 0-10, next 11-20,
 */
void LinkedCellContainer::initGrid() {
    //spdlog::info("this is initGrid");

    //periodic con
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
            spdlog::info("this is in initGrid periodic floor");
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
                for (auto &p1: getParticles()) {
                    if (p1.getX()[0] <= upperX && p1.getX()[0] >= lowerX &&
                        p1.getX()[1] <= upperY && p1.getX()[1] >= lowerY &&
                        p1.getX()[2] <= upperZ && p1.getX()[2] >= lowerZ) {
                        grid[x][y][z].push_back(p1);
                    }
                }
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
        if (upperZ + cutoffRadius > domainSize[2]) {
            upperZ += std::fmod(domainSize[2], cutoffRadius);
        } else {
            upperZ += cutoffRadius;
        }
    }
}

void LinkedCellContainer::resetF() {
    for (auto &p: getParticles()) {
        p.resetF();
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
//check done something wrong here -> in reflecting something wrong
void LinkedCellContainer::calculateF() {

    //spdlog::info("this is calcF");
    std::vector<Particle> nParticleList;
    initGrid();
    //down to up
    //fixed here rather than getGrid[0].size to getGrid[0][0].size -> runs now at least
    for (int z = 0; z < getGrid()[0][0].size(); ++z) {
        for (int y = 0; y < getGrid()[0].size(); ++y) {
            //left to right
            for (int x = 0; x < getGrid().size(); ++x) {
                for (auto p1 = grid[x][y][z].begin(); p1 < grid[x][y][z].end(); p1++) {
                    for (auto p2 = p1 + 1; p2 < grid[x][y][z].end(); p2++) {
                        double sigma;
                        if (p1->getType() != p2->getType()) {
                            sigma = ((p1->getSigma() + p2->getSigma()) / 2);
                        } else {
                            sigma = p1->getSigma();
                        }
                        auto dis = std::pow(2, (1 / 6)) * sigma;
                        std::array<double, 3> force{0.0, 0.0, 0.0};
                        if (isMembrane) {
                            if ((abs(p1->getX()[0] - p2->getX()[0]) && p1->getX()[1] == p2->getX()[1] &&
                                 p1->getX()[2] == p2->getX()[2])
                                || (abs(p1->getX()[1] - p2->getX()[1]) && p1->getX()[0] == p2->getX()[0] &&
                                    p1->getX()[2] == p2->getX()[2])
                                || (abs(p1->getX()[2] - p2->getX()[2]) && p1->getX()[1] == p2->getX()[1] &&
                                    p1->getX()[0] == p2->getX()[0])) {
                                force = applyNeighbouringForce(*p1, *p2, "direct");
                            } else if ((abs(p1->getX()[0] - p2->getX()[0]) && abs(p1->getX()[1] - p2->getX()[1]))
                                       || (abs(p1->getX()[0] - p2->getX()[0]) && abs(p1->getX()[2] - p2->getX()[2]))
                                       || (abs(p1->getX()[1] - p2->getX()[1]) && abs(p1->getX()[2] - p2->getX()[2]))) {
                                force = applyNeighbouringForce(*p1, *p2, "diagonal");
                            }
                        } else if (abs(p1->getX()[0] - p2->getX()[0]) < dis &&
                                       abs(p1->getX()[1] - p2->getX()[1]) < dis &&
                                       abs(p1->getX()[2] - p2->getX()[2]) < dis) {
                                force = forceModel.calculateForce(*p1, *p2);
                        }

                        spdlog::info("this is calcX ps position {}, {}, {}", p1->getX()[0], p1->getX()[1],
                                     p1->getX()[2]);
                        spdlog::info("this is calcX ps velocity{}, {}, {}", p1->getV()[0], p1->getV()[1],
                                     p1->getV()[2]);
                        spdlog::info("in calcF this is force of p1 before {}, {}, {}", p1->getF()[0], p1->getF()[1],
                                     p1->getF()[2]);
                        spdlog::info("in calcF this is calculated force {}, {}, {}", force[0], force[1], force[2]);
                        p1->addF({force[0], force[1] + (p1->getM() * gGrav), force[2]});
                        p2->addF({-1.0 * force[0], -1.0 * (force[1] + (p2->getM() * gGrav)), -1.0 * force[2]});
                        spdlog::info("in calcF this is force of p1 after {}, {}, {}", p1->getOldF()[0], p1->getOldF()[1],
                                     p1->getOldF()[2]);
                    }
                    // check done
                    if (checkBoundary('r')) {
                        applyReflecting(*p1);
                        spdlog::info("in calcF reflect this is force of p1 {}, {}, {}", p1->getF()[0], p1->getF()[1],
                                     p1->getF()[2]);
                    }
                    // check done
                    if (checkBoundary('p')) {
                        applyPeriodic(*p1, x, y, z);
                        spdlog::info("in calcF periodic this is force of p1 {}, {}, {}", p1->getF()[0], p1->getF()[1],
                                     p1->getF()[2]);
                    }
                }
                //if above neighbour exists (y+1)
                if (y + 1 > 0 && y + 1 < getGrid()[x].size()) {
                    //calculate force between this and above cell for all particles
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x][y + 1][z]) {
                            // spdlog::info("calcF, above");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
                //if above back neighbour exists (y+1 & z+1)
                if (y + 1 > 0 && y + 1 < getGrid()[x].size() && z + 1 > 0 && z + 1 < getGrid()[x][y].size()) {
                    //calculate force between this and above back cell for all particles
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x][y + 1][z + 1]) {
                            //  spdlog::info("calcF, above back");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
                //if right hand neighbour exists (x+1)
                if (x + 1 > 0 && x + 1 < getGrid().size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x + 1][y][z]) {
                            //   spdlog::info("calcF, right");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
                //if right back neighbour exists (x+1 & z+1)
                if (x + 1 > 0 && x + 1 < getGrid().size() && z + 1 > 0 && z + 1 < getGrid()[x][y].size()) {
                    //calculate force between this and above back cell for all particles
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x + 1][y][z + 1]) {
                            //   spdlog::info("calcF, right back");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
                //if right upper diagonal neighbour exists (x+1 & y+1)
                if (x + 1 > 0 && x + 1 < getGrid().size() && y + 1 > 0 && y + 1 < getGrid()[x].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x + 1][y + 1][z]) {
                            //spdlog::info("calcF, upper right");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
                //if right upper back diagonal neighbour exists (x+1 & y+1 & z+1)
                if (x + 1 > 0 && x + 1 < getGrid().size() && y + 1 > 0 && y + 1 < getGrid()[x].size()
                    && z + 1 > 0 && z + 1 < getGrid()[x][y].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x + 1][y + 1][z + 1]) {
                            //   spdlog::info("calcF, upper right diagonal back");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
                //if left upper diagonal neighbour exists (x-1 & y+1)
                if (x - 1 >= 0 && x - 1 < getGrid().size() && y + 1 > 0 && y + 1 < getGrid()[x].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x - 1][y + 1][z]) {
                            //   spdlog::info("calcF, upper left diagonal");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
                //if left upper back diagonal neighbour exists (x-1 & y+1 & z+1)
                if (x - 1 >= 0 && x - 1 < getGrid().size() && y + 1 > 0 && y + 1 < getGrid()[x].size() && z + 1 > 0 &&
                    z + 1 < getGrid()[x][y].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x - 1][y + 1][z + 1]) {
                            //  spdlog::info("calcF, upper left diagonal back");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
                //if back neighbour exists (z+1)
                if (z + 1 > 0 && z + 1 < getGrid()[x][y].size()) {
                    for (auto &p1: getGrid()[x][y][z]) {
                        for (auto &p2: getGrid()[x][y][z + 1]) {
                            //   spdlog::info("calcF, back");
                            auto force = forceModel.calculateForce(p1, p2);
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
                            p2.addF({-1.0 * force[0], -1.0 * (force[1] + (p2.getM() * gGrav)), -1.0 * force[2]});
                        }
                    }
                }
            }
        }
    }
    for (auto x: getGrid()) {
        for (auto y: x) {
            for (auto z: y) {
                for (auto p: z) {
                    nParticleList.push_back(p);
                }
            }
        }
    }
    particleList = nParticleList;
}

/**
 * @brief Calculate xi(tn+1)
 * formula: xi(tn+1) = xi(tn) + ∆t · vi(tn) + (∆t)^2 * (Fi(tn)/2mi)
 */
//check done, +front and back periodic
void LinkedCellContainer::calculateX(double delta_t) {
    // spdlog::info("this is calcX");

    for (auto &p: getParticles()) {
        auto xi_tn1 = Formulas::verletXStep(p.getX(), p.getV(), p.getF(), p.getM(), delta_t);
        //xi_tn1[2] = 0.0;
        //this is another periodic boundary same as in initGrid
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
            //spdlog::info("this is p[2] in initGrid periodic con floor {}", p.getX()[2]);
        }
        // front
        if (boundaryCon[4] == 'p' && p.getX()[2] < 0.0) {
            xi_tn1[2] = p.getX()[2] + domainSize[2];
        }
        // back
        if (boundaryCon[5] == 'p' && p.getX()[2] > domainSize[0]) {
            xi_tn1[2] = p.getX()[2] - domainSize[2];
        }
        p.setX(xi_tn1);  // Update the position
        spdlog::info("this is calcX ps position {}, {}, {}", p.getX()[0], p.getX()[1], p.getX()[2]);
        spdlog::info("this is calcX ps force {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
        spdlog::info("this is calcX ps velocity{}, {}, {}", p.getV()[0], p.getV()[1], p.getV()[2]);

    }
    //initGrid();
}

/**
 * @brief Calculate the velocity at time tn+1
 * formula: vi(tn+1) = vi(tn) + ∆t * ((Fi(tn) + Fi(tn+1))/ 2mi)
 */
//check done
void LinkedCellContainer::calculateV(double delta_t) {
    // spdlog::info("this is calcV");

    for (auto &p: getParticles()) {
        auto vi_tn1 = Formulas::verletVStep(p.getV(), p.getOldF(), p.getF(), p.getM(), delta_t);
        //vi_tn1[2] = 0.0;
        p.setV(vi_tn1);
        spdlog::info("this is calcX ps position {}, {}, {}", p.getX()[0], p.getX()[1], p.getX()[2]);
        spdlog::info("this is calcX ps force {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
        spdlog::info("this is calcX ps velocity{}, {}, {}", p.getV()[0], p.getV()[1], p.getV()[2]);
        // spdlog::info("in calcV this is velocity of p1 {}, {}, {}", p.getV()[0], p.getV()[1], p.getV()[2]);
    }
}

//check done
void LinkedCellContainer::plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::VTKWriter writer;
    writer.initializeOutput(size());
    for (auto &p: getParticles()) {
        writer.plotParticle(p);
    }
    writer.writeFile(out_name, iteration);
}

//check done
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

/**
 * @brief Apply reflecting condition on particles
 * @param p Particle to check
 * bounces back borderline particle from border through
 * creating new halo particle which applies force on it
 */
//check done
void LinkedCellContainer::applyReflecting(Particle &p) {
    auto dis = std::pow(2, (1 / 6)) * p.getSigma();

    Particle halo(p);
    //near left border
    if (boundaryCon[0] == 'r') {
        if (checkDistance(p, "left")) {
            spdlog::info("this is appReflect, left");
            halo.setX({(-1) * p.getX()[0], p.getX()[1], p.getX()[2]});
            spdlog::info("this is appReflect left ps pos {}, {}, {}", p.getX()[0], p.getX()[1], p.getX()[2]);
            spdlog::info("this is appReflect left halos pos {}, {}, {}", halo.getX()[0], halo.getX()[1],
                         halo.getX()[2]);
            spdlog::info("this is appReflect left ps force before {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
            auto force = forceModel.calculateForce(p, halo);
            spdlog::info("this is appReflect left force {}, {}, {}", force[0], force[1], force[2]);
            p.addF({force[0], force[1] + (p.getM() * gGrav), force[2]});
            spdlog::info("this is appReflect left ps force after {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
        }
    }
    //near upper border
    if (boundaryCon[1] == 'r') {
        if (checkDistance(p, "upper")) {
            spdlog::info("this is appReflect, upper");
            halo.setX({p.getX()[0], domainSize[1] + (domainSize[1] - p.getX()[1]), p.getX()[2]});
            spdlog::info("this is appReflect upper ps pos {}, {}, {}", p.getX()[0], p.getX()[1], p.getX()[2]);
            spdlog::info("this is appReflect upper halos pos {}, {}, {}", halo.getX()[0], halo.getX()[1],
                         halo.getX()[2]);
            spdlog::info("this is appReflect upper ps force before {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
            auto force = forceModel.calculateForce(p, halo);
            spdlog::info("this is appReflect upper force {}, {}, {}", force[0], force[1], force[2]);
            p.addF({force[0], force[1] + (p.getM() * gGrav), force[2]});
            spdlog::info("this is appReflect upper ps force after {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);

        }
    }
    //near right border
    if (boundaryCon[2] == 'r') {
        if (checkDistance(p, "right")) {
            spdlog::info("this is appReflect, right");
            halo.setX({domainSize[0] + (domainSize[0] - p.getX()[0]), p.getX()[1], p.getX()[2]});
            spdlog::info("this is appReflect right ps pos {}, {}, {}", p.getX()[0], p.getX()[1], p.getX()[2]);
            spdlog::info("this is appReflect right halos pos {}, {}, {}", halo.getX()[0], halo.getX()[1],
                         halo.getX()[2]);
            spdlog::info("this is appReflect right ps force before{}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
            auto force = forceModel.calculateForce(p, halo);
            spdlog::info("this is appReflect right force {}, {}, {}", force[0], force[1], force[2]);
            p.addF({force[0], force[1] + (p.getM() * gGrav), force[2]});
            spdlog::info("this is appReflect right ps force after {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
        }
    }
    //near floor border
    if (boundaryCon[3] == 'r') {
        if (checkDistance(p, "floor")) {
            spdlog::info("this is appReflect, floor");
            halo.setX({p.getX()[0], (-1) * p.getX()[1], p.getX()[2]});
            spdlog::info("this is appReflect floor ps pos {}, {}, {}", p.getX()[0], p.getX()[1], p.getX()[2]);
            spdlog::info("this is appReflect floor halos pos {}, {}, {}", halo.getX()[0], halo.getX()[1],
                         halo.getX()[2]);
            spdlog::info("this is appReflect floor ps force before {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
            auto force = forceModel.calculateForce(p, halo);
            spdlog::info("this is appReflect floor force {}, {}, {}", force[0], force[1], force[2]);
            p.addF({force[0], force[1] + (p.getM() * gGrav), force[2]});
            spdlog::info("this is appReflect floor ps force after {}, {}, {}", p.getF()[0], p.getF()[1], p.getF()[2]);
        }
    }
    if (getDimension() != 2) {
        //near front border
        if (boundaryCon[4] == 'r') {
            if (checkDistance(p, "front")) {
                spdlog::info("this is appReflect, front");
                //halo.setX({p.getX()[0], p.getX()[1], (-1) * p.getX()[2]});
                halo.setX({p.getX()[0], p.getX()[1], p.getX()[2] - dis});
                spdlog::info("this is appReflect front ps pos {}, {}, {}", p.getX()[0], p.getX()[1], p.getX()[2]);
                spdlog::info("this is appReflect front halos pos {}, {}, {}", halo.getX()[0], halo.getX()[1],
                             halo.getX()[2]);
                spdlog::info("this is appReflect front ps force before {}, {}, {}", p.getF()[0], p.getF()[1],
                             p.getF()[2]);
                auto force = forceModel.calculateForce(p, halo);
                spdlog::info("this is appReflect front force {}, {}, {}", force[0], force[1], force[2]);
                p.addF({force[0], force[1] + (p.getM() * gGrav), force[2]});
                spdlog::info("this is appReflect front ps force after {}, {}, {}", p.getF()[0], p.getF()[1],
                             p.getF()[2]);
            }
        }
        //near back border
        if (boundaryCon[5] == 'r') {
            if (checkDistance(p, "back")) {
                spdlog::info("this is appReflect, back");
                //changed from p.getX()[2] + (domainSize[2] - p.getX()[2]) -> domainSize[2] + (domainSize[2] - p.getX()[2])
                //halo.setX({p.getX()[0], p.getX()[1], domainSize[2] + (domainSize[2] - p.getX()[2])});
                halo.setX({p.getX()[0], p.getX()[1], p.getX()[2] + dis});
                spdlog::info("this is appReflect left ps pos {}, {}, {}", p.getX()[0], p.getX()[1], p.getX()[2]);
                spdlog::info("this is appReflect left halos pos {}, {}, {}", halo.getX()[0], halo.getX()[1],
                             halo.getX()[2]);
                spdlog::info("this is appReflect back ps force before {}, {}, {}", p.getF()[0], p.getF()[1],
                             p.getF()[2]);
                auto force = forceModel.calculateForce(p, halo);
                spdlog::info("this is appReflect back force {}, {}, {}", force[0], force[1], force[2]);
                p.addF({force[0], force[1] + (p.getM() * gGrav), force[2]});
                spdlog::info("this is appReflect back ps force after {}, {}, {}", p.getF()[0], p.getF()[1],
                             p.getF()[2]);
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
    //spdlog::info("this is appPeriod");
    Particle halo(p);

    // left
    if (boundaryCon[0] == 'p') {
        if (checkDistance(p, "left")) {
            // halo on right
            halo.setX({p.getX()[0] + domainSize[0], p.getX()[1], p.getX()[2]});
            for (auto &p1: getGrid()[getGrid().size() - 1][y][z]) {
                if (abs(p1.getX()[0] - halo.getX()[0]) <= cutoffRadius) {
                    auto force = forceModel.calculateForce(p1, halo);
                    p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                    p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                    p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                    p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                            p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                        p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                    p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
                    p1.addF({force[0], force[1] + (p1.getM() * gGrav), force[2]});
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
    //if ((abs(p1.getX()[0] - p2.getX()[0]) && p1.getX()[1] == p2.getX()[1] && p1.getX()[2] == p2.getX()[2])
    //    || (abs(p1.getX()[1] - p2.getX()[1]) && p1.getX()[0] == p2.getX()[0] && p1.getX()[2] == p2.getX()[2])
    //    || (abs(p1.getX()[2] - p2.getX()[2]) && p1.getX()[1] == p2.getX()[1] && p1.getX()[0] == p2.getX()[0])) {
    if (type == "direct") {
        return k * ((Formulas::secondNorm(p1.getX() - p2.getX())) - r0) * x;
    } // diagonal neighbors
        //else if ((abs(p1.getX()[0] - p2.getX()[0]) && abs(p1.getX()[1] - p2.getX()[1]))
        //         || (abs(p1.getX()[0] - p2.getX()[0]) && abs(p1.getX()[2] - p2.getX()[2]))
        //         || (abs(p1.getX()[1] - p2.getX()[1]) && abs(p1.getX()[2] - p2.getX()[2]))) {
    else if (type == "diagonal") {
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


