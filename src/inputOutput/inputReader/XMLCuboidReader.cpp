//
// Created by markus on 11/26/23.
//

#include "XMLCuboidReader.h"
#include "xsd/Simulation.hxx"

XMLCuboidReader::XMLCuboidReader() = default;

XMLCuboidReader::~XMLCuboidReader() = default;

std::vector<CuboidParticleGenerator> XMLCuboidReader::readFile(const char *filename) {

    std::vector<CuboidParticleGenerator> generators;
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    std::array<int, 3> N;
    double spacing;
    int type;

    std::unique_ptr<simulation> parameters = simulation_(filename);

    auto cuboids = parameters->cuboid();

    for(auto c: cuboids){

        x = {c.position().x(), c.position().y(), c.position().z()};
        v = {c.velocity().v(), c.velocity().w(), c.velocity().z()};
        m = c.mass();
        N = {c.grid().Nx(), c.grid().Ny(), c.grid().Nz()};
        spacing = c.spacing();
        type = c.type();

        generators.emplace_back(N, spacing, m, v, x, type);
    }

    return generators;
}
