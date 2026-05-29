// =====================================================================
//  model.h  --  declaration of the Model class for Section 4.1
//  (Simple Airframe/Autopilot Dynamics Model)
// =====================================================================

#ifndef EX1_MODEL_H
#define EX1_MODEL_H

#include "../../osk/osk.h"

class Model : public osk::Block {
public:
    // state, derivative, command, gain, velocity, intermediate signal,
    // and the constructor-supplied initial flight path angle.
    double gamma, gammad;
    double gamma_cmd;
    double a_cmd, k, v;
    double gamma0;

    Model(double gamma_);

    void init()   override;
    void update() override;
    void rpt()    override;
};

#endif
