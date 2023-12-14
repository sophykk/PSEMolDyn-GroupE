//
// Created by markus on 12/1/23.
//
#include "XMLFileReader.h"
#include "../../xsd/Simulation.hxx"

XMLFileReader::XMLFileReader() = default;

XMLFileReader::~XMLFileReader() = default;

void XMLFileReader::readSimulationParams(const char *filename, double &endTime, double &deltaT, std::string &modelType, std::string &containerType, std::string &objectType, int &plotInterval){

    std::unique_ptr<simulation> parameters = simulation_(filename);

    endTime = parameters->simulationParams().endTime();
    deltaT = parameters->simulationParams().deltaT();
    modelType = parameters->simulationParams().modelType();
    containerType = parameters->simulationParams().containerType();
    objectType = parameters->simulationParams().objectType();
    plotInterval = parameters->simulationParams().plotInterval();
}

void XMLFileReader::readLennardJonesForceParams(const char *filename, double &sigma, double &epsilon) {

    std::unique_ptr<simulation> parameters = simulation_(filename);

    sigma = parameters->lennardJonesForceParams()->sigma();
    epsilon = parameters->lennardJonesForceParams()->epsilon();
}

void XMLFileReader::readLinkedCellParams(const char *filename, std::vector<double> &x, double &cutOffR, std::array<char, 4> &boundaryC) {

    std::unique_ptr<simulation> parameters = simulation_(filename);

    x = {parameters->linkedCellParams()->domainSize().Lx(), parameters->linkedCellParams()->domainSize().Ly(), parameters->linkedCellParams()->domainSize().Lz()};
    cutOffR = parameters->linkedCellParams()->cutoffRadius();
    boundaryC = {parameters->linkedCellParams()->boundaryConditions().left()[0], parameters->linkedCellParams()->boundaryConditions().up()[0], parameters->linkedCellParams()->boundaryConditions().right()[0], parameters->linkedCellParams()->boundaryConditions().down()[0]};
}

std::vector<CuboidParticleGenerator> XMLFileReader::readCuboids(const char *filename) {

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

std::vector<SphereParticleGenerator> XMLFileReader::readSpheres(const char *filename) {

    std::vector<SphereParticleGenerator> generators;
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    double spacing;
    double r;
    int type;

    std::unique_ptr<simulation> parameters = simulation_(filename);

    auto spheres = parameters->sphere();

    for(auto s: spheres){

        x = {s.position().x(), s.position().y(), s.position().z()};
        v = {s.velocity().v(), s.velocity().w(), s.velocity().z()};
        m = s.mass();
        spacing = s.spacing();
        r = s.radius();
        type = s.type();

        generators.emplace_back(x,spacing,m,v,r,type);
    }

    return generators;
}