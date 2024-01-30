//
// Created by markus on 12/1/23.
//
#include "XMLFileReader.h"
#include "../../xsd/Simulation.hxx"

XMLFileReader::XMLFileReader() = default;

XMLFileReader::~XMLFileReader() = default;

/** Function that reads the general parameters used for the simulation out of an xml file */
void XMLFileReader::readSimulationParams(const char *filename, double &endTime, double &deltaT, std::string &modelType,
                                         std::string &containerType, std::string &objectType, int &plotInterval,
                                         bool &checkpointing, bool &useParallelization) {

    std::unique_ptr<simulation> parameters = simulation_(filename);

    endTime = parameters->simulationParams().endTime();
    deltaT = parameters->simulationParams().deltaT();
    modelType = parameters->simulationParams().modelType();
    containerType = parameters->simulationParams().containerType();
    objectType = parameters->simulationParams().objectType();
    plotInterval = parameters->simulationParams().plotInterval();
    checkpointing = parameters->simulationParams().checkpointing();
    useParallelization = parameters->simulationParams().useParallelization();
}

/** Function that reads the parameters for the linkedCellContainer out of an xml file */
void XMLFileReader::readLinkedCellParams(const char *filename, std::vector<double> &x, double &cutOffR,
                                         std::array<char, 6> &boundaryC, double &gGrav, bool &isMembrane) {

    std::unique_ptr<simulation> parameters = simulation_(filename);

    x = {parameters->linkedCellParams()->domainSize().Lx(), parameters->linkedCellParams()->domainSize().Ly(),
         parameters->linkedCellParams()->domainSize().Lz()};
    cutOffR = parameters->linkedCellParams()->cutoffRadius();
    boundaryC = {parameters->linkedCellParams()->boundaryConditions().left()[0],
                 parameters->linkedCellParams()->boundaryConditions().up()[0],
                 parameters->linkedCellParams()->boundaryConditions().right()[0],
                 parameters->linkedCellParams()->boundaryConditions().down()[0],
                 parameters->linkedCellParams()->boundaryConditions().front()[0],
                 parameters->linkedCellParams()->boundaryConditions().back()[0]
    };
    gGrav = parameters->linkedCellParams()->gravitationalAcceleration();
    isMembrane = parameters->linkedCellParams()->isMembrane();
}

void XMLFileReader::readMembraneParams(const char *filename, int &k, double &r0, double &pullUpF) {

    std::unique_ptr<simulation> parameters = simulation_(filename);

    k = parameters->membraneParams()->stiffness();
    r0 = parameters->membraneParams()->averageBond();
    pullUpF = parameters->membraneParams()->pullUpForce();
}

/** Function that reads the parameters for the thermostat initialization used for the simulation out of an xml file */
void XMLFileReader::readThermostatParams(const char *filename, double &initT, int &thermostatInterval) {

    std::unique_ptr<simulation> parameters = simulation_(filename);

    initT = parameters->thermostat()->initialTemperature();
    thermostatInterval = parameters->thermostat()->thermostatInterval();
}

/** Function that reads the parameters for the cuboid generators out of an xml file and initializes the generators */
std::vector<CuboidParticleGenerator> XMLFileReader::readCuboids(const char *filename, bool &isMembrane) {


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
    bool isWall;

    std::unique_ptr<simulation> parameters = simulation_(filename);

    auto cuboids = parameters->cuboid();

    for (auto c: cuboids) {

        x = {c.position().x(), c.position().y(), c.position().z()};
        v = {c.velocity().v(), c.velocity().w(), c.velocity().z()};
        m = c.mass();
        N = {c.grid().Nx(), c.grid().Ny(), c.grid().Nz()};
        spacing = c.spacing();
        sigma = c.sigma();
        epsilon = c.epsilon();
        gGrav = c.gravitationalAcceleration();
        type = c.type();
        isWall = c.isWall();

        generators.emplace_back(N, spacing, m, v, x, sigma, epsilon, gGrav, type, isMembrane, isWall);
    }

    return generators;
}

/** Function that reads the parameters for the sphere generators out of an xml file and initializes the generators */
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
    bool isWall;

    std::unique_ptr<simulation> parameters = simulation_(filename);

    auto spheres = parameters->sphere();

    for (auto s: spheres) {

        x = {s.position().x(), s.position().y(), s.position().z()};
        v = {s.velocity().v(), s.velocity().w(), s.velocity().z()};
        m = s.mass();
        spacing = s.spacing();
        r = s.radius();
        sigma = s.sigma();
        epsilon = s.epsilon();
        gGrav = s.gravitationalAcceleration();
        type = s.type();
        isWall = s.isWall();

        generators.emplace_back(x, spacing, m, v, r, sigma, epsilon, gGrav, type, isWall);
    }

    return generators;
}