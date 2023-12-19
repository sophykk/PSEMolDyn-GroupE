//
// Created by markus on 12/1/23.
//
#include "XMLFileReader.h"
#include "../../xsd/Simulation.hxx"
#include <iostream>

XMLFileReader::XMLFileReader() = default;

XMLFileReader::~XMLFileReader() = default;

void XMLFileReader::readSimulationParams(const char *filename, double &endTime, double &deltaT, std::string &modelType, std::string &containerType, std::string &objectType, int &plotInterval, bool &checkpointing){

    std::unique_ptr<simulation> parameters = simulation_(filename);

    endTime = parameters->simulationParams().endTime();
    deltaT = parameters->simulationParams().deltaT();
    modelType = parameters->simulationParams().modelType();
    containerType = parameters->simulationParams().containerType();
    objectType = parameters->simulationParams().objectType();
    plotInterval = parameters->simulationParams().plotInterval();
    checkpointing = parameters->simulationParams().checkpointing();
}

void XMLFileReader::readLinkedCellParams(const char *filename, std::vector<double> &x, double &cutOffR, std::array<char, 4> &boundaryC, double &gGrav) {

    std::unique_ptr<simulation> parameters = simulation_(filename);

    x = {parameters->linkedCellParams()->domainSize().Lx(), parameters->linkedCellParams()->domainSize().Ly(), parameters->linkedCellParams()->domainSize().Lz()};
    cutOffR = parameters->linkedCellParams()->cutoffRadius();
    boundaryC = {parameters->linkedCellParams()->boundaryConditions().left()[0], parameters->linkedCellParams()->boundaryConditions().up()[0], parameters->linkedCellParams()->boundaryConditions().right()[0], parameters->linkedCellParams()->boundaryConditions().down()[0]};
    gGrav = parameters->linkedCellParams()->gravitationalAcceleration();
}

void XMLFileReader::readThermostatParams(const char *filename, double &initT, int &thermostatInterval){

    std::unique_ptr<simulation> parameters = simulation_(filename);

    initT = parameters->thermostat()->initialTemperature();
    thermostatInterval = parameters->thermostat()->thermostatInterval();
}

std::vector<CuboidParticleGenerator> XMLFileReader::readCuboids(const char *filename) {


    std::vector<CuboidParticleGenerator> generators;
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    std::array<int, 3> N;
    double spacing;
    double sigma;
    double epsilon;
    double gGrav;
    int type;

    std::unique_ptr<simulation> parameters = simulation_(filename);

    auto cuboids = parameters->cuboid();

    for(auto c : cuboids){

        x = {c.position().x(), c.position().y(), c.position().z()};
        v = {c.velocity().v(), c.velocity().w(), c.velocity().z()};
        m = c.mass();
        N = {c.grid().Nx(), c.grid().Ny(), c.grid().Nz()};
        spacing = c.spacing();
        sigma = c.sigma();
        epsilon = c.epsilon();
        gGrav = c.gravitationalAcceleration();
        type = c.type();

        generators.emplace_back(N, spacing, m, v, x, sigma, epsilon, gGrav, type);
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
    double sigma;
    double epsilon;
    double gGrav;
    int type;

    std::unique_ptr<simulation> parameters = simulation_(filename);

    auto spheres = parameters->sphere();

    for(auto s: spheres){

        x = {s.position().x(), s.position().y(), s.position().z()};
        v = {s.velocity().v(), s.velocity().w(), s.velocity().z()};
        m = s.mass();
        spacing = s.spacing();
        r = s.radius();
        sigma = s.sigma();
        epsilon = s.epsilon();
        gGrav = s.gravitationalAcceleration();
        type = s.type();

        generators.emplace_back(x,spacing,m,v,r,sigma,epsilon,gGrav,type);
    }

    return generators;
}