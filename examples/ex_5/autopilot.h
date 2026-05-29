// =====================================================================
//  autopilot.h  --  Section 4.4 (PDF page 52)
//
//  The autopilot is one half of the decomposed missile-steering
//  system.  It reads gamma from the missile model (via the
//  ACCESS_FN-generated gamma_() function on Missile) and produces an
//  acceleration command a_cmd that the missile then consumes.  The
//  connection between the two objects is established in main.cpp by
//  calling getsFrom() on each.
// =====================================================================

#ifndef EX4_AUTOPILOT_H
#define EX4_AUTOPILOT_H

#include "../../osk/osk.h"

// Forward declaration: the Missile class is defined in missile.h.
// We only need a pointer to it here, so a forward declaration is
// enough -- this avoids a circular #include between autopilot.h and
// missile.h.
class Missile;

class Autopilot : public osk::Block {
public:
    Autopilot(double k);

    // Cross-block wiring API.  main.cpp will call this with a Missile*
    // so the autopilot knows which missile to read gamma from.
    void getsFrom(Missile* obj) { this->missile = obj; }

    // Expose a_cmd to other blocks.  Generates "double a_cmd_() const".
    ACCESS_FN(double, a_cmd);

    void init()   override;
    void update() override;
    void rpt()    override;

protected:
    double gamma_cmd;
    double k;
    double a_cmd;
    Missile* missile;
};

#endif
