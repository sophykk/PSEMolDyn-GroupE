//
// Created by markus on 1/30/24.
//
#include "CSVFileExporter.h"
#include <string>

CSVFileExporter::CSVFileExporter(ParticleContainerBase& container) {

    numBins = 50;
    fileNumber = 1;
    binSize = container.getDomainSize()[0] / numBins;

    // initialize density and averageVelocity vectors
    density.resize(numBins, 0.0);
    averageVelocity.resize(numBins, 0.0);
}

CSVFileExporter::~CSVFileExporter() = default;

void CSVFileExporter::updateAndWriteStatistics(ParticleContainerBase& container) {

    // reset statistics
    std::fill(density.begin(), density.end(), 0.0);
    std::fill(averageVelocity.begin(), averageVelocity.end(), 0.0);

    // iterate over particles and update statistics
    for (const Particle& particle : container.getParticles()) {
        // Determine the bin index based on the particles x-coordinate
        int binIndex = static_cast<int>(particle.getX()[0] / binSize);

        // update density and velocity for the corresponding bin
        density[binIndex] += 1.0;
        averageVelocity[binIndex] += particle.getV()[0];
    }

    // normalize averageVelocity by the number of particles in each bin
    for (int i = 0; i < averageVelocity.size(); ++i) {
        if (density[i] > 0) {
            averageVelocity[i] /= density[i];
        }
    }

    // write statistics to CSV file
    writeStatisticsToCSV();
}

void CSVFileExporter::writeStatisticsToCSV() {
    // open CSV file
    std::ofstream outputFile("statistics_" + std::to_string(fileNumber) + ".csv");

    // write header
    outputFile << "Bin Index, Density, Average Velocity\n";

    // write data for each bin
    for (int i = 0; i < numBins; ++i) {
        outputFile << i << "," << density[i] << "," << averageVelocity[i] << "\n";
    }

    // close the file
    outputFile.close();

    // increase the file number to avoid overwriting
    fileNumber ++;
}
