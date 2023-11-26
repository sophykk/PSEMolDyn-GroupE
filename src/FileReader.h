/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "particle/Particle.h"
#include "particle/ParticleGenerator.h"

#include <list>
#include <vector>


class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  void readFile(std::vector<ParticleGenerator> &gen, char *filename);
};
