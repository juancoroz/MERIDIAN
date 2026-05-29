// =====================================================================
//  autopilot.h  --  Section 4.6.3.4 (Comprehensive Sample Demo)
//
//  Same Autopilot class as Section 4.4/4.5 but with two additions
//  for the comprehensive sample demo:
//    - te member: rolling event time for the async fine-stepping
//      window [1.55, 1.63] in stage 2
//    - update() declares two extra "throwaway" EVENT samples at
//      1.6001 and 1.6002 just to demonstrate that the async scheduler
//      can capture closely spaced one-time events
// =====================================================================

#ifndef EX64_AUTOPILOT_H
#define EX64_AUTOPILOT_H

#include "../../osk/osk.h"

class Missile;

class Autopilot : public osk::Block {
public:
    Autopilot(double k);

    void getsFrom(Missile* obj) { this->missile = obj; }

    ACCESS_FN(double, a_cmd);

    void init()   override;
    void update() override;
    void rpt()    override;

protected:
    double gamma_cmd;
    double k;
    double a_cmd;
    double te;          // event-time accumulator for the fine-step window
    Missile* missile;
};

#endif
