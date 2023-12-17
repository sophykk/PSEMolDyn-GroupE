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


LinkedCellContainer2::LinkedCellContainer2(ForceBase &model, std::vector<double> &dSize, double &cRadius,
                                           std::array<char, 4> bCon) :
        ParticleContainerBase(model), domainSize(dSize), cutoffRadius(cRadius), boundaryCon(bCon) {
    int cellDimensionX;
    int cellDimensionY;
    double gGrav = -12.44;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? cellDimensionX = domainSize[0] / cutoffRadius
                                                : cellDimensionX = std::ceil(domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? cellDimensionY = domainSize[1] / cutoffRadius
                                                : cellDimensionY = std::ceil(domainSize[1] / cutoffRadius);
    grid.resize(cellDimensionX, std::vector<std::vector<Particle>>(cellDimensionY, std::vector<Particle>()));
    initGrid();
}

LinkedCellContainer2::LinkedCellContainer2(ForceBase &model, std::vector<Particle> &particles,
                                           std::vector<double> &dSize, double &cRadius, std::array<char, 4> bCon) :
        ParticleContainerBase(model), particleList(particles), domainSize(dSize), cutoffRadius(cRadius),
        boundaryCon(bCon) {
    int cellDimensionX;
    int cellDimensionY;
    double gGrav = -12.44;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? cellDimensionX = domainSize[0] / cutoffRadius
                                                : cellDimensionX = std::ceil(
            domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? cellDimensionY = domainSize[1] / cutoffRadius
                                                : cellDimensionY = std::ceil(
            domainSize[1] / cutoffRadius);
    grid.resize(cellDimensionX, std::vector<std::vector<Particle>>(cellDimensionY, std::vector<Particle>()));
    initGrid();
}

std::vector<Particle> &LinkedCellContainer2::getParticles() {
    return particleList;
}

std::vector<std::vector<std::vector<Particle>>> &LinkedCellContainer2::getGrid() {
    return grid;
}

std::vector<Particle> &LinkedCellContainer2::getCell(std::array<double, 3> &pos) {
    return getGrid()[std::ceil(pos[0] / cutoffRadius)][std::ceil(pos[1] / cutoffRadius)];
}

char &LinkedCellContainer2::getBoundaryCon(int index) {
    return boundaryCon[index];
}

/**
 * @brief Sets the boundary
 * @param boundary a char, should be 'o' or 'r'
 * @throws error message invalid_argument if boundary is not 'o' or 'r'.
 */
void LinkedCellContainer2::setBoundaryCon(std::array<char, 4> &nBoundary) {
    for (int i = 0; i < nBoundary.size(); i++) {
        if (nBoundary[i] == 'o' || nBoundary[i] == 'p' || nBoundary[i] == 'r') {
            boundaryCon[i] = nBoundary[i];
        } else {
            spdlog::error("Index {} is not a valid boundary condition!", i);
        }
    }
}

bool LinkedCellContainer2::checkBoundary(char b) {
    for (int i = 0; i < boundaryCon.size(); ++i) {
        if (getBoundaryCon(i) == b) {
            return true;
        }
    }
    return false;
}

void LinkedCellContainer2::addParticle(Particle &particle) {
    particleList.push_back(particle);
}

std::size_t LinkedCellContainer2::size() const {
    return particleList.size();
}

/**
 * @brief This methods puts all the particles in a cell in the grid according to their positions
 * if cutoffRadius = 10, and dSize = {180, 90, 1}
 * then first cell from 0-10, next 11-20,
 * */
void LinkedCellContainer2::initGrid() {

    double minDouble = std::numeric_limits<double>::min();
    double lowerY = 0.0;
    double upperY = cutoffRadius;
    //x0,y --> 0,0 --> 0,1
    for (int y = 0; y < std::ceil(domainSize[1] / cutoffRadius); ++y) {
        double lowerX = 0.0;
        double upperX = cutoffRadius;
        //x,y0 --> 0,0 --> 1,0 --> 2,0 --> 3,0
        for (int x = 0; x < std::ceil(domainSize[0] / cutoffRadius); ++x) {
            grid[x][y].clear();
            for (auto &p1: getParticles()) {
                if (p1.getX()[0] <= upperX && p1.getX()[0] >= lowerX &&
                    p1.getX()[1] <= upperY && p1.getX()[1] >= lowerY) {
                    grid[x][y].push_back(p1);
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
}

void LinkedCellContainer2::resetF() {
    for (auto &p: getParticles()) {
        p.resetF();
    }
}

/**
 * @brief calculates force between all particles in cell with itself and its neighbours
 * start in the left corner cell, only has  3 neighbours
 * always max 4 calcs for each cell because
 * e.g u start bottom left and calc above, right and left diagonal
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
                    std::array<double, 3> force2 = {force[0], force[1] * p2->getM(), force[2]};
                    force[1] += gGrav * p1->getM();
                    p1->addF(force);
                    p2->addF(-1.0 * force2);
                }
                if (checkBoundary('r')) {
                    applyReflecting(*p1);
                }
                if (checkBoundary('p')) {
                    applyPeriodic(*p1);
                }
            }
            //if above neighbour exists
            if (y + 1 >= 0 && y + 1 < getGrid()[x].size()) {
                //calculate force between this and above cell for all particles
                for (auto &p1: getGrid()[x][y]) {
                    for (auto &p2: getGrid()[x][y + 1]) {
                        auto force = forceModel.calculateForce(p1, p2);
                        std::array<double, 3> force2 = {force[0], gGrav * p2.getM(), force[2]};
                        force[1] += gGrav * p1.getM();
                        p1.addF(force);
                        p2.addF(-1.0 * force2);
                    }
                }
                //if right hand neighbour exists
                if (x + 1.0 >= 0.0 && x + 1.0 < getGrid().size()) {
                    for (auto &p1: getGrid()[x][y]) {
                        for (auto &p2: getGrid()[x + 1][y]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            std::array<double, 3> force2 = {force[0], gGrav * p2.getM(), force[2]};
                            force[1] += gGrav * p1.getM();
                            p1.addF(force);
                            p2.addF(-1.0 * force2);
                        }
                    }
                }
                //if upper right diagonal neighbour exists
                if (y + 1 >= 0 && y + 1 < getGrid()[x].size() && x + 1 >= 0 && x + 1 < getGrid().size()) {
                    for (auto &p1: getGrid()[x][y]) {
                        for (auto &p2: getGrid()[x + 1][y + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            std::array<double, 3> force2 = {force[0], gGrav * p2.getM(), force[2]};
                            force[1] += gGrav * p1.getM();
                            p1.addF(force);
                            p2.addF(-1.0 * force2);
                        }
                    }
                }
                //if upper left diagonal neighbour exists
                if (y + 1 >= 0 && y + 1 < getGrid()[x].size() && x - 1 >= 0 && x - 1 < getGrid().size()) {
                    for (auto &p1: getGrid()[x][y]) {
                        for (auto &p2: getGrid()[x - 1][y + 1]) {
                            auto force = forceModel.calculateForce(p1, p2);
                            std::array<double, 3> force2 = {force[0], gGrav * p2.getM(), force[2]};
                            force[1] += gGrav * p1.getM();
                            p1.addF(force);
                            p2.addF(-1.0 * force2);
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

/**
 * @brief if cutOffRadius is not cellSize, check distance then
 * @code
 * bool LinkedCellContainer2::withinCutoff(Particle &p1, Particle &p2) const {
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

/**
 * @brief Calculate xi(tn+1)
 * formula: xi(tn+1) = xi(tn) + ∆t · vi(tn) + (∆t)^2 * (Fi(tn)/2mi)
 */
void LinkedCellContainer2::calculateX(double delta_t) {
    for (auto &p: getParticles()) {
        auto xi_tn1 = Formulas::verletXStep(p.getX(), p.getV(), p.getF(), p.getM(), delta_t);

        if (p.getX()[0] < 0.0 && boundaryCon[0] == 'p') {
            xi_tn1[0] = xi_tn1[0] + domainSize[0];
        }
        if (p.getX()[1] > domainSize[1] && boundaryCon[1] == 'p') {
            xi_tn1[1] = xi_tn1[1] - domainSize[1];
        }
        if (p.getX()[0] > domainSize[0] && boundaryCon[2] == 'p') {
            xi_tn1[0] = xi_tn1[0] - domainSize[0];
        }
        if (p.getX()[1] < 0.0 && boundaryCon[3] == 'p') {
            xi_tn1[1] = xi_tn1[1] + domainSize[1];
        }

        p.setX(xi_tn1);  // Update the position
    }
    initGrid();
}

/**
 * @brief Calculate the velocity at time tn+1
 * formula: vi(tn+1) = vi(tn) + ∆t * ((Fi(tn) + Fi(tn+1))/ 2mi)
 */
void LinkedCellContainer2::calculateV(double delta_t) {
    for (auto &p: getParticles()) {
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

bool LinkedCellContainer2::checkDistance(Particle &p, std::string border) {
    //1.0 is sigma (not mentioned because doesn't change result)
    //divided by 2.0 because of border p____|____p
    auto disCheck = std::pow(2, 1.0 / 6) / 2.0;

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
    return false;
}

/**
 * @brief Apply reflecting condition on particles
 * @param p Particle to check
 * bounces back borderline particle from border through
 * creating new halo particle which applies force on it
 */
void LinkedCellContainer2::applyReflecting(Particle &p) {

    Particle halo(p);
    //near left border
    if (boundaryCon[0] == 'r') {
        if (checkDistance(p, "left")) {
            halo.setX({(-1) * p.getX()[0], p.getX()[1], p.getX()[2]});
            auto force = forceModel.calculateForce(p, halo);
            force[1] += gGrav * p.getM();
            p.addF(force);
        }
    }
    //near right border
    if (boundaryCon[2] == 'r') {
        if (checkDistance(p, "right")) {
            halo.setX({domainSize[0] + (domainSize[0] - p.getX()[0]), p.getX()[1], p.getX()[2]});
            auto force = forceModel.calculateForce(p, halo);
            force[1] += gGrav * p.getM();
            p.addF(force);
        }
    }
    //near floor border
    if (boundaryCon[3] == 'r') {
        if (checkDistance(p, "floor")) {
            halo.setX({p.getX()[0], (-1) * p.getX()[1], p.getX()[2]});
            auto force = forceModel.calculateForce(p, halo);
            force[1] += gGrav * p.getM();
            p.addF(force);
        }
    }
    //near upper border
    if (boundaryCon[1] == 'r') {
        if (checkDistance(p, "upper")) {
            halo.setX({p.getX()[0], domainSize[1] + (domainSize[1] - p.getX()[1]), p.getX()[2]});
            auto force = forceModel.calculateForce(p, halo);
            force[1] += gGrav * p.getM();
            p.addF(force);
        }
    }
}

void LinkedCellContainer2::applyPeriodic(Particle &p) {
    auto distance = std::pow(2, 1.0 / 6);
    Particle halo(p);

    // left
    if (boundaryCon[0] == 'p') {
        if (p.getX()[0] < 0.0) {
            //insert particle if outside boundary on the opposite side
            p.setX({p.getX()[0] + domainSize[0], p.getX()[1], p.getX()[2]});
        }
        /*
        //creating halos if near border
        if (checkDistance(p, "left")) {
            halo.setX({p.getX()[0] + domainSize[0], p.getX()[1], p.getX()[2]});
            for (auto &p1: getGrid()[getGrid().size()][std::ceil(p.getX()[1] / cutoffRadius)]) {
                if (distance >= abs(p1.getX()[0] - halo.getX()[0])) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[1] += gGrav + p.getM();
                    p1.addF(force);
                }
            }
        }*/
    }

    // upper
    if (boundaryCon[1] == 'p') {
        if (p.getX()[1] > domainSize[1]) {
            std::array<double, 3> x_new = {p.getX()[0], p.getX()[1] - domainSize[1], p.getX()[2]};
            p.setX(x_new);
        }
        /*
        if (checkDistance(p, "upper")) {
            halo.setX({p.getX()[0], p.getX()[1] - domainSize[1], p.getX()[2]});
            for (auto &p1: getGrid()[p.getX()[0] / cutoffRadius][std::ceil(p.getX()[1])]) {
                if (distance >= abs(p1.getX()[1] - halo.getX()[1])) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[1] += gGrav + p.getM();
                    p1.addF(force);
                }
            }
        }*/
    }
    // right
    if (boundaryCon[2] == 'p') {
        if (p.getX()[0] > domainSize[0]) {
            std::array<double, 3> x_new = {p.getX()[0] - domainSize[0], p.getX()[1], p.getX()[2]};
            p.setX(x_new);
        }/*
        if (checkDistance(p, "right")) {
            halo.setX({p.getX()[0] - domainSize[0], p.getX()[1], p.getX()[2]});
            for (auto &p1: getGrid()[0][p.getX()[1] / cutoffRadius]) {
                if (distance >= abs(p1.getX()[0] - halo.getX()[0])) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[1] += gGrav + p.getM();
                    p1.addF(force);
                }
            }
        }*/
    }
    // floor
    if (boundaryCon[3] == 'p') {
        if (p.getX()[1] < 0.0) {
            std::array<double, 3> x_new = {p.getX()[0], p.getX()[1] + domainSize[1], p.getX()[2]};
            p.setX(x_new);
        }/*
        if (checkDistance(p, "floor")) {
            halo.setX({p.getX()[0] + domainSize[0], p.getX()[1], p.getX()[2]});
            for (auto &p1: getGrid()[p.getX()[0] / cutoffRadius][getGrid()[0].size()]) {
                if (distance >= abs(p1.getX()[1] - halo.getX()[1])) {
                    auto force = forceModel.calculateForce(p1, halo);
                    force[1] += gGrav + p.getM();
                    p1.addF(force);
                }
            }
        }*/
    }
}

