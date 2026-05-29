// =====================================================================
//  model.cpp  --  Section 4.6.3.3 (PDF page 68)
//
//  main passes the nominal dt=0.1 to the kernel.  This model declares
//  a rolling event time te that starts at 0.55, advances by 0.01 each
//  time the event fires, and stops advancing once te crosses 0.6.
//  The result: the kernel uses 0.1 s steps everywhere EXCEPT the
//  0.55-0.60 window where it uses 0.01 s steps.
// =====================================================================

#include "model.h"
#include <cstdio>
#include <iostream>

using namespace osk;
using std::cout;
using std::endl;

// EPS absorbs floating-point noise in the (te <= 0.6 + EPS) bound
// test, matching the report's EPS shorthand on PDF page 68.
static const double EPS = 1.0e-9;

Model::Model(double gamma_) {
    gamma0 = gamma_;
    addIntegrator(gamma, gammad);
}

void Model::init() {
    cout << "starting Model...\n";
    if (initCount == 0) {
        gamma = gamma0;
        k     = 1000.0;
        te    = 0.55;
    }
}

void Model::update() {
    if (State::sample(State::EVENT, te)) {
        te += 0.01;
        if (te <= 0.6 + EPS) {
            State::sample(State::EVENT, te);   // continue fine stepping
        }
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
