// =====================================================================
//  model.cpp  --  Section 4.1 model, parameter-driven via Filer
//
//  Same equations as the original (gammad = k*(gamma_cmd - gamma)/v),
//  but values for gamma0, gamma_cmd, k, and v are read from the
//  "Model" block of input.txt at construction time.  The output of
//  this program should match PDF page 43 byte-for-byte, demonstrating
//  that the io utility cleanly satisfies Appendix 4.4.
// =====================================================================

#include "model.h"
#include <cstdio>
#include <iostream>

using namespace osk;

Model::Model(Filer& ff) {
    ff.setLine0("Model");
    gamma0    = ff.getDouble("gamma0");
    gamma_cmd = ff.getDouble("gamma_cmd");
    k         = ff.getDouble("k");
    v         = ff.getDouble("v");
    addIntegrator(gamma, gammad);
}

void Model::init() {
    std::cout << "starting Model...\n";
    if (initCount == 0) {
        gamma = gamma0;
    }
}

void Model::update() {
    a_cmd  = k * (gamma_cmd - gamma);
    gammad = a_cmd / v;
}

void Model::rpt() {
    if (State::sample(0.1)) {
        std::printf("%8.3f %8.6f\n", State::t, gamma);
    }
}
