// =====================================================================
//  model.h  --  declaration of the Model class for Section 4.2.1
//  (Discrete-Sampler variant of the Section 4.1 first-order servo)
// =====================================================================

#ifndef EX2_MODEL_H
#define EX2_MODEL_H

#include "../../osk/osk.h"

class Model : public osk::Block {
public:
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
