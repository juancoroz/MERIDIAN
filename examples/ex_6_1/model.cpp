// =====================================================================
//  model.cpp  --  Section 4.6.3.1 (PDF page 65)
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

    // Three one-shot events at irregular sub-step times.  None of
    // 0.15, 0.185, 0.1855 is an even increment of dt=0.1 -- the
    // asynchronous scheduler in the OSK is what makes the kernel
    // actually land on these times so the events can fire.
    if (State::sample(State::EVENT, 0.1855)) {
        k = 1000;
        cout << " EVENT " << State::t << " " << k << endl;
    }
    if (State::sample(State::EVENT, 0.15)) {
        k = 1000;
        cout << " EVENT " << State::t << " " << k << endl;
    }
    if (State::sample(State::EVENT, 0.185)) {
        k = 1000;
        cout << " EVENT " << State::t << " " << k << endl;
    }

    a_cmd  = k * (gamma_cmd - gamma);
    v      = 1000.0;
    gammad = a_cmd / v;
}

void Model::rpt() {
    if (State::sample(0.1)) {
        std::printf("%8.3f %8.4f %8.6f\n",
                    State::t, State::dt, this->gamma);
    }
}
