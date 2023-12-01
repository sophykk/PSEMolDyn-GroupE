//
// Created by sophy on 30.11.2023.
//

#include "LinkedCellContainer2.h"
#include <vector>
#include "Formulas.h"
#include "utils/ArrayUtils.h"
#include "inputOutput/outputWriter/VTKWriter.h"
#include <spdlog/spdlog.h>
#include <cmath>
#include <limits>


LinkedCellContainer2::LinkedCellContainer2(ForceBase &model, std::vector<double> &dSize, double &cRadius) :
        ParticleContainerBase(model), domainSize(dSize), cutoffRadius(cRadius) {
    int x;
    int y;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? x = domainSize[0] / cutoffRadius : x = std::ceil(
            domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? y = domainSize[1] / cutoffRadius : y = std::ceil(
            domainSize[1] / cutoffRadius);

    grid.resize(x, std::vector<std::vector<Particle>>(y, std::vector<Particle>()));
    initGrid();
   // spdlog::info("x is {}", getGrid().size());
   // spdlog::info("y is {}", getGrid()[0].size());
}

LinkedCellContainer2::LinkedCellContainer2(ForceBase &model, std::vector<Particle> &particles,
                                           std::vector<double> &dSize, double &cRadius) :
        ParticleContainerBase(model), particleList(particles), domainSize(dSize), cutoffRadius(cRadius) {
    int x;
    int y;
    std::fmod(domainSize[0], cutoffRadius) == 0 ? x = domainSize[0] / cutoffRadius : x = std::ceil(
            domainSize[0] / cutoffRadius);
    std::fmod(domainSize[1], cutoffRadius) == 0 ? y = domainSize[1] / cutoffRadius : y = std::ceil(
            domainSize[1] / cutoffRadius);

    grid.resize(x, std::vector<std::vector<Particle>>(y, std::vector<Particle>()));
    initGrid();
  //  spdlog::info("x is {}", getGrid().size());
  //  spdlog::info("y is {}", getGrid()[0].size());
}

std::vector<Particle> &LinkedCellContainer2::getParticles() {
    return particleList;
}

std::vector<std::vector<std::vector<Particle>>> &LinkedCellContainer2::getGrid() {
    return grid;
}

/*
//just 2D first
std::vector<Particle> &LinkedCellContainer2::getParticlesFromCell(int x, int y) {
    return getGrid()[x][y];
}
*/

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
    spdlog::info("We are now in calcF");
    //down to up
    for (int y = 0; y < getGrid()[0].size(); ++y) {
        //    spdlog::info("first for loop");
        //left to right
        for (int x = 0; x < getGrid().size(); ++x) {
            //      spdlog::info("second for loop");
            //if above neigbour exists
            if (y + 1 >= 0 && y + 1 < getGrid()[x].size()) {
                //     spdlog::info("first if stmt");
                //calculate force between this and above cell for all particles
                for (auto &p1: getGrid()[x][y]) {
                    //       spdlog::info("first for loop in first if");
                    for (auto &p2: getGrid()[x][y + 1]) {
                        //         spdlog::info("second for loop in first if");
                        auto force = forceModel.calculateForce(p1, p2);
                        p1.addF(force);
                        p2.addF(-1.0 * force);
                        //       spdlog::debug("above force for Particle {} is {}, {}, {}", p1.toString(), p1.getF()[0],
                        //                   p1.getF()[1], p1.getF()[2]);
                        //
                    }
                }
            }
            //if right hand neighbour exists
            if (x + 1 >= 0 && x + 1 < getGrid().size()) {
                //     spdlog::info("second if stmt");

                for (auto &p1: getGrid()[x][y]) {
                    //       spdlog::info("first for loop in second if");

                    for (auto &p2: getGrid()[x + 1][y]) {
                        //         spdlog::info("second for loop in second if");

                        auto force = forceModel.calculateForce(p1, p2);
                        p1.addF(force);
                        p2.addF(-1.0 * force);
                        //        spdlog::debug("right force for Particle {} is {}, {}, {}", p1.toString(), p1.getF()[0],
                        //                      p1.getF()[1], p1.getF()[2]);
                    }
                }
            }
            //if upper right diagonal neighbour exists
            if (y + 1 >= 0 && y + 1 < getGrid()[x].size() && x + 1 >= 0 && x + 1 < getGrid().size()) {
                //    spdlog::info("third if stmt");

                for (auto &p1: getGrid()[x][y]) {
                    //    spdlog::info("first for loop in third if");

                    for (auto &p2: getGrid()[x + 1][y + 1]) {
                        //        spdlog::info("second for loop in third if");

                        auto force = forceModel.calculateForce(p1, p2);
                        p1.addF(force);
                        p2.addF(-1.0 * force);
                        //        spdlog::debug("right diagonal force for Particle {} is {}, {}, {}", p1.toString(), p1.getF()[0],
                        //                      p1.getF()[1], p1.getF()[2]);
                    }
                }
            }
            //if upper left diagonal neighbour exists
            if (y + 1 >= 0 && y + 1 < getGrid()[x].size() && x - 1 >= 0 && x - 1 < getGrid().size()) {
                //    spdlog::info("fourth if stmt");

                for (auto &p1: getGrid()[x][y]) {
                    //    spdlog::info("first for loop in fourth if");

                    for (auto &p2: getGrid()[x - 1][y]) {
                        //    spdlog::info("second for loop in fourth if");
                        auto force = forceModel.calculateForce(p1, p2);
                        p1.addF(force);
                        p2.addF(-1.0 * force);
                        //    spdlog::debug("left diagonal force for Particle {} is {}, {}, {}", p1.toString(), p1.getF()[0],
                        //                   p1.getF()[1], p1.getF()[2]);
                    }
                }
            }
        }
    }
    spdlog::info("plist has {} particles", particleList.size());
    for (auto x : getGrid()) {
        for (auto y : x) {
            for (auto z : y) {
                spdlog::info("This is {}", z.toString());
            }}}
}

/*
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
    spdlog::info("We are now in calcX");
    initGrid();
    spdlog::info("in cell 0, 0 are {} p", getGrid()[0][0].size());
    spdlog::info("in cell 1, 0 are {} p", getGrid()[1][0].size());
    spdlog::info("in cell 2, 0 are {} p", getGrid()[2][0].size());
    spdlog::info("in cell 3, 0 are {} p", getGrid()[3][0].size());
    spdlog::info("in cell 0, 1 are {} p", getGrid()[0][1].size());
    spdlog::info("in cell 1, 1 are {} p", getGrid()[1][1].size());
    spdlog::info("in cell 2, 1 are {} p", getGrid()[2][1].size());
    spdlog::info("in cell 3, 1 are {} p", getGrid()[3][1].size());
    spdlog::info("plist has {} particles", getParticles().size());
    for (auto x : getParticles()) {
        spdlog::info("This is {}", x.toString());
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

    spdlog::info("We are now in initGrid");

    double minDouble = 1.0;
            //std::numeric_limits<double>::min();
    double lowery = 0.0;
    double uppery = cutoffRadius;
    //y --> 0,1 --> 1,1
    for (int i = 0; i < std::ceil(domainSize[1] / cutoffRadius); ++i) {
        //2. dimension
        double lowerx = 0.0;
        double upperx = cutoffRadius;
        //x, --> 0,0 --> 1,0 --> 2,0 --> 3,0
        for (int j = 0; j < std::ceil(domainSize[0] / cutoffRadius); ++j) {
            //only for 2D domains {_, _, 1}, otherwise need other implementation (extra loop)
            //3. dimension
            int u = 0;
            for (auto &p1: getParticles()) {
                spdlog::debug("x cellSize now {} - {}", lowerx, upperx);
                spdlog::debug("y cellSize now {} - {}", lowery, uppery);
                if (p1.getX()[0] <= upperx && p1.getX()[0] >= lowerx &&
                    p1.getX()[1] <= uppery && p1.getX()[1] >= lowery) {
                    spdlog::debug("we are in grid {}, {}", j, i);
                    spdlog::info("this is the particle before {}", p1.toString());
                    grid[j][i].push_back(p1);
                    spdlog::info("this is p1 after set to grid {}", grid[j][i][u].toString());
                    u++;
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
    spdlog::info("in cell 0, 0 are {} p", getGrid()[0][0].size());
    spdlog::info("in cell 1, 0 are {} p", getGrid()[1][0].size());
    spdlog::info("in cell 2, 0 are {} p", getGrid()[2][0].size());
    spdlog::info("in cell 3, 0 are {} p", getGrid()[3][0].size());
    spdlog::info("in cell 0, 1 are {} p", getGrid()[0][1].size());
    spdlog::info("in cell 1, 1 are {} p", getGrid()[1][1].size());
    spdlog::info("in cell 2, 1 are {} p", getGrid()[2][1].size());
    spdlog::info("in cell 3, 1 are {} p", getGrid()[3][1].size());
    spdlog::info("plist has {} particles", particleList.size());
    for (auto x : getGrid()) {
        for (auto y : x) {
            for (auto z : y) {
            spdlog::info("This is {}", z.toString());
    }}}
}