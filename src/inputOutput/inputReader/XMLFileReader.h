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

    void readSimulationParams(const char *filename, double &endTime, double &deltaT, std::string &modelType, std::string &containerType, std::string &objectType, int &plotInterval);

    void readLennardJonesForceParams(const char *filename, double &sigma, double &epsilon);

    void readLinkedCellParams(const char *filename, std::vector<double> &x, double &cutOffR, std::array<char, 4> &boundaryC);

    std::vector<CuboidParticleGenerator> readCuboids(const char *filename);

    std::vector<SphereParticleGenerator> readSpheres(const char *filename);
};
