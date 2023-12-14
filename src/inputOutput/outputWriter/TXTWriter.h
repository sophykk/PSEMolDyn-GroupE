//
// Created by markus on 12/14/23.
//

#ifndef PSEMOLDYN_TXTWRITER_H
#define PSEMOLDYN_TXTWRITER_H

#pragma once

#include "Particle.h"

#include <fstream>
#include <vector>

namespace outputWriter {

    class TXTWriter {

    public:
        TXTWriter();

        virtual ~TXTWriter();

        void plotParticles(std::vector<Particle> &particles, const std::string &fileName);
    };

}

#endif //PSEMOLDYN_TXTWRITER_H
