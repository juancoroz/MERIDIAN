// =====================================================================
//  model.h  --  Appendix 2 test problem
//
//    d²y/dt² + t*y = 0,  y(0) = 1, y'(0) = 0
//
//  State-space:
//    x1 = y         x1' = x2
//    x2 = y'        x2' = -t * x1
//
//  Expected at t = 4.0: x1 ~= 0.219970 (exact), 0.219972 (RK2 dt=0.001)
// =====================================================================

#ifndef EX_APP2_MODEL_H
#define EX_APP2_MODEL_H

#include "../../osk/osk.h"
#include "state_rk2.h"

class Model : public osk::Block {
public:
    double x1, x1d;
    double x2, x2d;
    double x10, x20;

    Model(double x10_, double x20_) : x10(x10_), x20(x20_) {
        x1 = x1d = x2 = x2d = 0.0;
        // Use the RK2 integrator for both states via the templated
        // addIntegrator form.  No change to the kernel was needed --
        // State_rk2 is just a State subclass that overrides
        // propagate() and stages().
        addIntegrator<State_rk2>(x1, x1d);
        addIntegrator<State_rk2>(x2, x2d);
    }

    void init() override {
        if (initCount == 0) {
            std::printf("starting Model...\n");
            x1 = x10;
            x2 = x20;
        }
    }

    void update() override {
        x1d =  x2;
        x2d = -osk::State::t * x1;
    }

    void rpt() override {
        if (osk::State::sample(0.8) || osk::State::ticklast) {
            std::printf("%8.3f %8.6f\n", osk::State::t, this->x1);
        }
    }
};

#endif
