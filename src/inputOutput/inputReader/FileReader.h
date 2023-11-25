/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "Particle.h"
#include "ParticleGenerator.h"

#include <list>
#include <vector>


class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  std::vector<ParticleGenerator> readFile(char *filename);
};
