// =====================================================================
//  autopilot.cpp  --  Section 4.6.3.4 (PDF pages 73-74)
// =====================================================================

#include "autopilot.h"
#include "missile.h"
#include <cstdio>

using namespace osk;

static const double EPS = 1.0e-9;

Autopilot::Autopilot(double k_) {
    k         = k_;
    gamma_cmd = 1.0;
    a_cmd     = 0.0;
    te        = 1.55;
    missile   = nullptr;
}

void Autopilot::init() {
    if (initCount == 0) {
        gamma_cmd = 1.0;
        te        = 1.55;
    }
}

void Autopilot::update() {
    // Note: unlike Section 4.4, a_cmd is recomputed every step (not
    // gated by sample(0.1)).  The autopilot is now a continuous
    // proportional controller; only the rpt method throttles output.
    a_cmd = k * (gamma_cmd - missile->gamma_());

    // Two "throwaway" events at finer resolution than any nominal dt,
    // declared as statement-form to demonstrate they get captured by
    // the async scheduler.  These appear in the missile rpt output
    // because the missile uses no-arg sample() to print every step.
    State::sample(State::EVENT, 1.6002);
    State::sample(State::EVENT, 1.6001);

    // Rolling-event fine-stepping window [1.55, 1.63] -- forces dt=0.01
    // inside the window.  Outside the window, the nominal dt=0.1 from
    // main.cpp prevails.
    if (State::sample(State::EVENT, te)) {
        if (te <= 1.63 - EPS) {
            te += 0.01;
        }
        State::sample(State::EVENT, te);
    }
}

void Autopilot::rpt() {
    if (State::sample(0.2)) {
        std::printf("Autopilot %6.4f %8.6f\n", State::t, this->a_cmd);
    }
}
