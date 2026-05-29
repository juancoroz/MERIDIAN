// =====================================================================
//  model.cpp  --  Section 4.3 model implementation (PDF page 49)
//
//  Same as the Section 4.2.1 model except for the addition of a
//  one-shot scheduled event at t = 1 s that flips the gain k from
//  1000 to 2000.  Two State::sample() forms appear in the same
//  update() method:
//
//      State::sample(State::EVENT, 1.0)  one-shot, fires at t = 1
//      State::sample(0.1)                periodic, fires every 0.1 s
//
//  Expected output is on PDF page 49 -- the "EVENT 1 2000" log line
//  appears between the t=0.9 and t=1.0 sample rows, and gamma's
//  trajectory after t=1 reflects the doubled gain (final gamma
//  reaches 0.962561 instead of 0.878423).
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
        k = 1000.0;
    }
}

void Model::update() {
    gamma_cmd = 1.0;
    v         = 1000.0;

    // One-shot gain switch at t = 1 s.
    if (State::sample(State::EVENT, 1.0)) {
        k = 2000;
        cout << " EVENT " << State::t << " " << k << endl;
    }

    // 10 Hz sample-and-hold of the acceleration command.
    if (State::sample(0.1)) {
        a_cmd = k * (gamma_cmd - gamma);
        cout << " sample " << State::t << " " << a_cmd << endl;
    }

    gammad = a_cmd / v;
}

void Model::rpt() {
    if (State::sample(0.1)) {
        std::printf("%8.3f %8.6f\n", State::t, this->gamma);
    }
}
