// =====================================================================
//  missile.cpp  --  Section 4.6.3.4 (PDF page 74)
//
//  Same staging logic as Section 4.5 but rpt() uses the no-arg form
//  of sample() so the missile reports at the start of EVERY
//  integrating step.  This is what makes the async stepping visible
//  in the output -- we see the kernel land precisely on 1.55, 1.56,
//  ..., 1.5999, 1.6, 1.6001, 1.6002 etc.
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

void Missile::init() {
    if (initCount == 0) {
        gamma  = gamma0;
        v      = 1000.0;
        istage = 0;
    } else if (initCount == 1) {
        istage = 1;
    }
}

void Missile::update() {
    if (State::sample(State::EVENT, 1.00)) {
        Sim::stop = 1;
    }
    if (istage == 0) {
        gammad = 0.0;
    } else if (istage == 1) {
        gammad = autopilot->a_cmd_() / v;
    }
}

void Missile::rpt() {
    if (State::sample()) {       // no-arg = every integrating step
        std::printf("Missile %8.4f %8.4f %8.6f\n",
                    State::t, State::dt, this->gamma);
    }
}
