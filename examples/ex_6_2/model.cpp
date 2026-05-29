// =====================================================================
//  model.cpp  --  Section 4.6.3.2 (PDF pages 66-67)
//
//  main passes dt=1.0 to the kernel, but the model declares an event
//  at "te" via sample(EVENT,...), then -- whenever that event fires
//  -- bumps te by 0.1 and immediately re-declares the next event.
//  The asynchronous scheduler shrinks each step from the nominal 1.0
//  down to 0.1 to land on the rolling event time.
// =====================================================================

#include "model.h"
#include <cstdio>
#include <iostream>

using namespace osk;
using std::cout;
using std::endl;

Model::Model(double gamma_) {
    gamma0 = gamma_;
    addIntegrator(gamma, gammad);
}

void Model::init() {
    cout << "starting Model...\n";
    if (initCount == 0) {
        gamma = gamma0;
        k     = 1000.0;
        te    = 0.0;
    }
}

void Model::update() {
    if (State::sample(State::EVENT, te)) {
        te += 0.1;
        // Statement form: immediately tell the kernel about the next
        // event time so this step's actual_dt can include it.  Note
        // that this is a function CALL (return value discarded), not
        // a special syntax -- C++ allows discarding any function's
        // return value.
        State::sample(State::EVENT, te);
        cout << "**** sample " << State::t << " " << State::dt << endl;
    }
    gamma_cmd = 1.0;
    a_cmd     = k * (gamma_cmd - gamma);
    v         = 1000.0;
    gammad    = a_cmd / v;
}

void Model::rpt() {
    if (State::sample()) {
        std::printf("%8.3f %8.4f %8.6f\n",
                    State::t, State::dt, this->gamma);
    }
}
