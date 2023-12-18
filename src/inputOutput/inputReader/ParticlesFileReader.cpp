//
// Created by markus on 12/14/23.
//
#include "ParticlesFileReader.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <array>

ParticlesFileReader::ParticlesFileReader() = default;

ParticlesFileReader::~ParticlesFileReader() = default;

void ParticlesFileReader::readFile(std::vector<Particle> &particles, std::string &filename) {
    std::array<double, 3> x;
    std::array<double, 3> v;
    std::array<double, 3> f;
    std::array<double, 3> old_f;
    double m;
    int type;
    int num_particles = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {

        particles.clear();

        getline(input_file, tmp_string);
        std::cout << "Read line: " << tmp_string << std::endl;

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            std::cout << "Read line: " << tmp_string << std::endl;
        }

        std::istringstream numstream(tmp_string);
        numstream >> num_particles;
        std::cout << "Reading " << num_particles << "." << std::endl;
        getline(input_file, tmp_string);
        std::cout << "Read line: " << tmp_string << std::endl;

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
                std::cout
                        << "Error reading file: eof reached unexpectedly reading from line "
                        << i << std::endl;
                exit(-1);
            }
            datastream >> m;
            datastream >> type;
            //todo read sigma and eps and change here
            particles.emplace_back(x, v, f, old_f, m, 1.0, 5.0, type);

            getline(input_file, tmp_string);
            std::cout << "Read line: " << tmp_string << std::endl;
        }
    } else {
        std::cout << "Error: could not open file " << filename << std::endl;
        exit(-1);
    }
}