// =====================================================================
//  model.h  --  declaration of the Model class for Section 4.3
//  (Scheduled-Event variant: gain scheduler triggered at t = 1 s)
//
//  Structurally identical to the Section 4.2.1 model -- the only
//  change is in update() where a one-shot State::sample(State::EVENT,
//  1.0) block flips k from 1000 to 2000 at t = 1 s.  No new fields
//  are needed because k is already a member of Model.
// =====================================================================

#ifndef EX3_MODEL_H
#define EX3_MODEL_H

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
