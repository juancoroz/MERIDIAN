// =====================================================================
//  model.cpp  --  Section 4.2.1 model implementation (PDF page 44)
//
//  Identical to the Section 4.1 model except a 10 Hz sample-and-hold
//  has been added for the lateral acceleration command (a_cmd).  The
//  command is recomputed only every 0.1 second using State::sample()
//  inside update(); between samples, a_cmd retains its previous
//  value (sample-and-hold behavior).
//
//  Expected output is documented on PDF page 45 and shows the
//  interleaving of "sample <t> <a_cmd>" log lines from update()
//  with "<t> <gamma>" rows from rpt().  Both fire at exactly the
//  same set of times (0, 0.1, 0.2, ..., 2.0).
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
