//
// Created by markus on 11/25/23.
//

#ifndef PSEMOLDYN_GROUPE_FORCEBASE_H
#define PSEMOLDYN_GROUPE_FORCEBASE_H

#include <array>
#include "Particle.h"

class ForceBase {
public:

    virtual ~ForceBase() = default;

    virtual std::array<double, 3> calculateForce(Particle& p1, Particle& p2) = 0;
};
#endif //PSEMOLDYN_GROUPE_FORCEBASE_H
