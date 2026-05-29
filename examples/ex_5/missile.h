// =====================================================================
//  missile.h  --  Section 4.5 (PDF page 58-59)
//
//  Same Missile class as Section 4.4 but with two additions for
//  multi-stage operation:
//    - an istage member variable that the model uses to decide which
//      derivative law to apply (0 = open-loop ballistic, 1 = closed-
//      loop under autopilot control)
//    - init() uses the kernel's initCount to set istage according to
//      which stage entry is currently being initialised
//    - update() polls for the staging event with State::sample(EVENT,
//      ...) and sets Sim::stop=1 to advance the train-of-objects
// =====================================================================

#ifndef EX5_MISSILE_H
#define EX5_MISSILE_H

#include "../../osk/osk.h"

class Autopilot;

class Missile : public osk::Block {
public:
    Missile(double gamma0);

    void getsFrom(Autopilot* obj) { this->autopilot = obj; }

    ACCESS_FN(double, gamma);

    void init()   override;
    void update() override;
    void rpt()    override;

protected:
    double v;
    double gamma0, gamma, gammad;
    int    istage;            // which stage we're currently in
    Autopilot* autopilot;
};

#endif
