//
// Created by markus on 11/25/23.
//

#ifndef PSEMOLDYN_GROUPE_LENNARDJONESFORCE_H
#define PSEMOLDYN_GROUPE_LENNARDJONESFORCE_H

#include "ForceBase.h"

class LennardJonesForce : public ForceBase {
public:

    LennardJonesForce();

    std::array<double, 3> calculateForce(Particle& p1, Particle& p2) override;

    void setSigma(const double nSigma);

    void setEps(const double nEps);

private:
    double sigma;
    double eps;
};


#endif //PSEMOLDYN_GROUPE_LENNARDJONESFORCE_H
