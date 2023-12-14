//
// Created by markus on 12/14/23.
//
#include "TXTWriter.h"
#include <iomanip>
#include <sstream>
#include <iostream>

namespace outputWriter {

    TXTWriter::TXTWriter() = default;

    TXTWriter::~TXTWriter() = default;

    void TXTWriter::plotParticles(std::vector<Particle> &particles, const std::string &fileName) {
        std::ofstream outputFile(fileName);

        if (!outputFile.is_open()) {
            std::cerr << "Error opening file: " << fileName << std::endl;
            return;
        }

        // Write the size of the particle list in the first line
        outputFile << particles.size() << std::endl;

        // Write each particle on a new line
        for (auto &particle: particles) {
            outputFile << particle.getX()[0] << ' ' << particle.getX()[1] << ' ' << particle.getX()[2] << ' '
                       << particle.getV()[0] << ' ' << particle.getV()[1] << ' ' << particle.getV()[2] << ' '
                       << particle.getF()[0] << ' ' << particle.getF()[1] << ' ' << particle.getF()[2] << ' '
                       << particle.getOldF()[0] << ' ' << particle.getOldF()[1] << ' ' << particle.getOldF()[2] << ' '
                       << particle.getM() << ' ' << particle.getType() << std::endl;
        }

        outputFile.close();
    }
}