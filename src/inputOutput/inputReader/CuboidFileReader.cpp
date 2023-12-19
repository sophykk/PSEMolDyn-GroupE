/*
 * CuboidFileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "CuboidFileReader.h"
#include "Generators/CuboidParticleGenerator.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include "spdlog/spdlog.h"


CuboidFileReader::CuboidFileReader() = default;

CuboidFileReader::~CuboidFileReader() = default;

std::vector<CuboidParticleGenerator> CuboidFileReader::readFile(char *filename) {
    std::vector<CuboidParticleGenerator> generators;
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    double spacing;
    double sigma;
    double epsilon;
    double gGrav;
    int type;
    int num_cuboids = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {

        getline(input_file, tmp_string);
        spdlog::debug("Read line: {}", tmp_string);

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            spdlog::debug("Read line: {}", tmp_string);
        }

        std::istringstream numstream(tmp_string);
        numstream >> num_cuboids;
        spdlog::debug("Reading {}", num_cuboids);
        /**
       * start reading dispositions, velocities and mass of particles
       */
        getline(input_file, tmp_string);
        spdlog::debug("Read line: {}", tmp_string);

        for (int i = 0; i < num_cuboids; i++) {
            std::istringstream datastream(tmp_string);

            for (auto &xj : x) {
                datastream >> xj;
            }
            for (auto &vj : v) {
                datastream >> vj;
            }

            if (datastream.eof()) {
                spdlog::error("Error reading file: eof reached unexpectedly reading from line {}", i);
                exit(-1);
            }

            datastream >> m;
            /**If there is N (Grid), generate the Grid first and then pass it to particles*/
            std::array<int, 3> N; // Array to store N values

            for (auto &Nj : N) {
                datastream >> Nj;
            }

            datastream >> spacing;
            datastream >> sigma;
            datastream >> epsilon;
            datastream >> gGrav;
            datastream >> type;

            generators.emplace_back(N, spacing, m, v, x, sigma, epsilon, gGrav,type);

            getline(input_file, tmp_string);
            spdlog::debug("Read line: {}", tmp_string);
        }
        return generators;
    } else {
        spdlog::error("Error: could not open file {}", filename);
        exit(-1);
    }
}