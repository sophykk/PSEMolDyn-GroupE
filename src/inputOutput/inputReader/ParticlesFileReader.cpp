//
// Created by markus on 12/14/23.
//
#include "ParticlesFileReader.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <array>
#include <spdlog/spdlog.h>

ParticlesFileReader::ParticlesFileReader() = default;

ParticlesFileReader::~ParticlesFileReader() = default;

void ParticlesFileReader::readFile(std::vector<Particle> &particles, std::string &filename) {
    std::array<double, 3> x;
    std::array<double, 3> v;
    std::array<double, 3> f;
    std::array<double, 3> old_f;
    double m;
    double sigma;
    double epsilon;
    double gGrav;
    double membrane;
    int type;
    int num_particles = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {

        particles.clear();

        getline(input_file, tmp_string);

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
        }

        std::istringstream numstream(tmp_string);
        numstream >> num_particles;
        getline(input_file, tmp_string);

        for (int i = 0; i < num_particles; i++) {
            std::istringstream datastream(tmp_string);

            for (auto &xj : x) {
                datastream >> xj;
            }
            for (auto &vj : v) {
                datastream >> vj;
            }
            for (auto &fj : f) {
                datastream >> fj;
            }
            for (auto &old_fj : old_f) {
                datastream >> old_fj;
            }
            if (datastream.eof()) {
                spdlog::error("Error reading file: eof reached unexpectedly reading from line {}", i);
                exit(-1);
            }
            datastream >> m;
            datastream >> sigma;
            datastream >> epsilon;
            datastream >> gGrav;
            datastream >> membrane;
            datastream >> type;

            particles.emplace_back(x, v, f, old_f, m, gGrav, sigma, epsilon, membrane, type);

            getline(input_file, tmp_string);
        }
    } else {
        spdlog::error("Error: could not open file {}", filename);
        exit(-1);
    }
}