//
// Created by markus on 12/1/23.
//

#pragma once

#include "Particle.h"
#include "Generators/CuboidParticleGenerator.h"
#include "Generators/SphereParticleGenerator.h"

#include <list>
#include <vector>
#include <array>

class XMLFileReader {

public:
    XMLFileReader();

    virtual ~XMLFileReader();

    void readSimulationParams(const char *filename, double &endTime, double &deltaT, std::string &modelType,
                              std::string &containerType, std::string &objectType, int &plotInterval,
                              bool &checkpointing, int &task);

    void readLinkedCellParams(const char *filename, std::vector<double> &x, double &cutOffR, std::array<char, 6> &boundaryC,
                         double &gGrav, bool &isMembrane);

    void readMembraneParams(const char *filename, int &k, double &r0, double &pullUpF);

    void readThermostatParams(const char *filename, double &initT, int &thermostatInterval);

    std::vector<CuboidParticleGenerator> readCuboids(const char *filename, bool &isMembrane);

    std::vector<SphereParticleGenerator> readSpheres(const char *filename);
};
