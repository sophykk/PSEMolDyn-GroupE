//
// Created by markus on 11/25/23.
//

#ifndef PSEMOLDYN_GROUPE_GRAVITATIONALFORCE_H
#define PSEMOLDYN_GROUPE_GRAVITATIONALFORCE_H

#include "ForceBase.h"

class GravitationalForce : public ForceBase {
public:

    std::array<double, 3> calculateForce(Particle& p1, Particle& p2) override;
};

#endif //PSEMOLDYN_GROUPE_GRAVITATIONALFORCE_H
