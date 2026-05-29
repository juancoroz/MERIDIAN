// =====================================================================
//  model.h  --  Section 4.6.3.1 (PDF page 65)
//  Capturing asynchronous events: three one-shot events fire at times
//  0.15, 0.185, 0.1855 s -- none of which are even multiples of the
//  0.1 s integrating step.
// =====================================================================

#ifndef EX6_1_MODEL_H
#define EX6_1_MODEL_H

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
