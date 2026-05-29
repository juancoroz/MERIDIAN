// =====================================================================
//  missile.cpp  --  Section 4.5 (PDF pages 58-59)
// =====================================================================

#include "missile.h"
#include "autopilot.h"
#include <cstdio>

using namespace osk;

Missile::Missile(double gamma0_) {
    gamma0 = gamma0_;
    gamma  = 0.0;
    gammad = 0.0;
    v      = 1000.0;
    istage = 0;
    autopilot = nullptr;
    addIntegrator(gamma, gammad);
}

// Section 4.5.3: initCount lets a model tailor its setup per stage
// entry.  This missile appears in BOTH stages of the train-of-objects,
// so init() runs twice -- once before stage 0 (initCount=0) and once
// before stage 1 (initCount=1).  We set istage from initCount so the
// model knows which derivative law to apply in update().
void Missile::init() {
    if (initCount == 0) {
        gamma  = gamma0;
        v      = 1000.0;
        istage = 0;
    } else if (initCount == 1) {
        istage = 1;
    }
}

// Section 4.5.4: poll for the staging event and select the derivative
// law for the current stage.
void Missile::update() {
    // One-shot staging event at t = 1.0 s.  When this fires we ask
    // the kernel to advance to stage index 1 -- the closed-loop
    // configuration with the autopilot in the loop.
    if (State::sample(State::EVENT, 1.00)) {
        Sim::stop = 1;
    }

    // Apply the derivative law for whichever stage we're in.
    if (istage == 0) {
        gammad = 0.0;                          // open-loop ballistic
    } else if (istage == 1) {
        gammad = autopilot->a_cmd_() / v;      // closed-loop
    }
}

void Missile::rpt() {
    if (State::sample(0.1)) {
        std::printf("Missile %.3f %.6f\n", State::t, this->gamma);
    }
}
