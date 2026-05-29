// =====================================================================
//  autopilot.cpp  --  Section 4.4 (PDF page 53)
// =====================================================================

#include "autopilot.h"
#include "missile.h"            // need the full Missile class for gamma_()
#include <cstdio>

using namespace osk;

Autopilot::Autopilot(double k_) {
    k = k_;
    gamma_cmd = 1.0;
    a_cmd     = 0.0;
    missile   = nullptr;
}

void Autopilot::init() {
    // No state of our own to initialise -- a_cmd starts at 0 and is
    // overwritten on the very first 10 Hz sample.
}

void Autopilot::update() {
    if (State::sample(0.1)) {
        // Read gamma from the missile via its access function.  The
        // access function returns by value, so we cannot accidentally
        // mutate the missile's state from here.
        a_cmd = k * (gamma_cmd - missile->gamma_());
    }
}

void Autopilot::rpt() {
    if (State::sample(0.1)) {
        std::printf("Autopilot %.3f %.6f\n", State::t, this->a_cmd);
    }
}
