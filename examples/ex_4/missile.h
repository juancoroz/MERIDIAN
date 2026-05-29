// =====================================================================
//  missile.h  --  Section 4.4 (PDF page 53)
//
//  The missile carries the airframe state -- a single scalar gamma
//  whose derivative is a_cmd / v.  It reads a_cmd from the autopilot
//  via its access function.
// =====================================================================

#ifndef EX4_MISSILE_H
#define EX4_MISSILE_H

#include "../../osk/osk.h"

class Autopilot;   // fwd; full definition lives in autopilot.h

class Missile : public osk::Block {
public:
    Missile(double gamma0);

    void getsFrom(Autopilot* obj) { this->autopilot = obj; }

    // Expose gamma to other blocks.  Generates "double gamma_() const".
    ACCESS_FN(double, gamma);

    void init()   override;
    void update() override;
    void rpt()    override;

protected:
    double v;
    double gamma0, gamma, gammad;
    Autopilot* autopilot;
};

#endif
