// =====================================================================
//  model.h  --  Section 4.6.3.2 (PDF page 66)
//  Use sample to generate finer time-steps than the nominal dt.
// =====================================================================

#ifndef EX6_2_MODEL_H
#define EX6_2_MODEL_H

#include "../../osk/osk.h"

class Model : public osk::Block {
public:
    double gamma, gammad;
    double gamma_cmd;
    double a_cmd, k, v;
    double gamma0;
    double te;            // rolling event time -- declared on each fire

    Model(double gamma_);

    void init()   override;
    void update() override;
    void rpt()    override;
};

#endif
