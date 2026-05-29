// =====================================================================
//  model.h  --  Section 4.6.3.3 (PDF page 68)
//  Use sample to control time-step over a window (0.55-0.6 s).
// =====================================================================

#ifndef EX6_3_MODEL_H
#define EX6_3_MODEL_H

#include "../../osk/osk.h"

class Model : public osk::Block {
public:
    double gamma, gammad;
    double gamma_cmd;
    double a_cmd, k, v;
    double gamma0;
    double te;

    Model(double gamma_);

    void init()   override;
    void update() override;
    void rpt()    override;
};

#endif
