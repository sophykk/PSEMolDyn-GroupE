//
// Created by markus on 11/26/23.
//
/*
 * CuboidFileReader.h
 */

#pragma once

#include "Particle.h"
#include "Generators/CuboidParticleGenerator.h"

#include <list>
#include <vector>

class XMLCuboidReader {

public:
    XMLCuboidReader();
    virtual ~XMLCuboidReader();

    std::vector<CuboidParticleGenerator> readFile(const char *filename);
};




