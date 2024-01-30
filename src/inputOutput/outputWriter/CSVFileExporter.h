//
// Created by markus on 1/30/24.
//

#ifndef BACHELORPRAKTIKUM_CSVFILEEXPORTER_H
#define BACHELORPRAKTIKUM_CSVFILEEXPORTER_H

#include "../../Containers/LinkedCellContainer.h"
#include <fstream>
#include <vector>

class CSVFileExporter {
private:

    std::vector<double> density;
    std::vector<double> averageVelocity;
    int numBins;
    int fileNumber;
    double binSize;

public:
    CSVFileExporter(ParticleContainerBase& container);
    virtual ~CSVFileExporter();

    void updateAndWriteStatistics(ParticleContainerBase& container);
    void writeStatisticsToCSV();
};

#endif //BACHELORPRAKTIKUM_CSVFILEEXPORTER_H
