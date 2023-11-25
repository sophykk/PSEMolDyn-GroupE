/*
 * CuboidFileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "Particle.h"
#include "Generators/CuboidParticleGenerator.h"

#include <list>
#include <vector>


class CuboidFileReader {

public:
  CuboidFileReader();
  virtual ~CuboidFileReader();

  std::vector<CuboidParticleGenerator> readFile(char *filename);
};
