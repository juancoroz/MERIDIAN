// =====================================================================
//  missile.cpp  --  Section 4.4 (PDF pages 53-54)
// =====================================================================

#include "missile.h"
#include "autopilot.h"          // need the full Autopilot class for a_cmd_()
#include <cstdio>

using namespace osk;

Missile::Missile(double gamma0_) {
    gamma0 = gamma0_;
    v      = 1000.0;
    gamma  = 0.0;
    gammad = 0.0;
    autopilot = nullptr;
    addIntegrator(gamma, gammad);
}

void Missile::init() {
    if (initCount == 0) {
        gamma = gamma0;
    }
}

void Missile::update() {
    gammad = autopilot->a_cmd_() / v;
}

void Missile::rpt() {
    if (State::sample(0.1)) {
        std::printf("Missile %.3f %.6f\n", State::t, this->gamma);
    }
}
