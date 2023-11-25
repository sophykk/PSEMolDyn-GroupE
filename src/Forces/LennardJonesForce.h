//
// Created by markus on 11/25/23.
//

#ifndef PSEMOLDYN_GROUPE_LENNARDJONESFORCE_H
#define PSEMOLDYN_GROUPE_LENNARDJONESFORCE_H

#include "ForceBase.h"

class LennardJonesForce : public ForceBase {
public:

    LennardJonesForce(double sigma, double eps);

    std::array<double, 3> calculateForce(Particle& p1, Particle& p2) override;

private:
    double sigma;
    double eps;
};


#endif //PSEMOLDYN_GROUPE_LENNARDJONESFORCE_H
