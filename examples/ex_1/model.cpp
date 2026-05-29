// =====================================================================
//  model.cpp  --  Section 4.1 model implementation (PDF page 42)
//
//  Closed-loop first-order system:
//      e        = gamma_cmd - gamma
//      a_cmd    = K * e
//      gammad   = a_cmd / V
//      gamma    = integral( gammad ) dt
//
//  Gain K=1000, V=1000 -> tau = V/K = 1 s.  Step input gamma_cmd=1
//  starting from gamma(0)=0 gives the analytical response
//      gamma(t) = 1 - exp(-t)
//  and in particular gamma(1) = 1 - 1/e = 0.6321...
// =====================================================================

#include "model.h"
#include <cstdio>
#include <iostream>

using namespace osk;
using std::cout;

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
    a_cmd     = k * (gamma_cmd - gamma);
    gammad    = a_cmd / v;
}

void Model::rpt() {
    if (State::sample(0.1)) {
        std::printf("%8.3f %8.6f\n", State::t, this->gamma);
    }
}
