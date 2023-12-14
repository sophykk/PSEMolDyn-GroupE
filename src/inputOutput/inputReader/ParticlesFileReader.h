//
// Created by markus on 12/14/23.
//

#ifndef PSEMOLDYN_PARTICLESFILEREADER_H
#define PSEMOLDYN_PARTICLESFILEREADER_H

#pragma once

#include "Particle.h"

#include <vector>

class ParticlesFileReader {

public:
    ParticlesFileReader();
    virtual ~ParticlesFileReader();

    void readFile(std::vector<Particle> &particles, std::string &filename);
};

#endif //PSEMOLDYN_PARTICLESFILEREADER_H
