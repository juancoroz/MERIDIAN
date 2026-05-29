// =====================================================================
//  state_rk2.h  --  Midpoint Runge-Kutta (RK2) integrator
//
//  Demonstrates Appendix 2's "plug in your own integrator" pattern
//  using my kernel's existing extensibility hooks:
//
//    1. Inherit from osk::State.
//    2. Override stages() to declare how many derivative evaluations.
//    3. Override propagate(stage, dt) for the algorithm itself.
//    4. Register with Block::addIntegrator<State_rk2>(state, deriv).
//
//  No changes to State, Block, or Sim are required.  The report's
//  Appendix 2 design uses a factory/prototype pattern and exposes
//  more State internals (kpass, ready, t1, etc.) -- this version
//  uses virtual dispatch instead, which is more idiomatic C++14.
//
//  Algorithm (per Appendix 2 page A-16):
//
//    y* = y_n + (dt/2) * f(y_n, t_n)         <-- evaluate at start
//    y_{n+1} = y_n + dt * f(y*, t_n + dt/2)  <-- evaluate at midpoint
// =====================================================================

#ifndef EX_APP2_STATE_RK2_H
#define EX_APP2_STATE_RK2_H

#include "../../osk/state.h"

class State_rk2 : public osk::State {
public:
    State_rk2(double& state, double& deriv) : osk::State(state, deriv) {
        // Two-stage method needs two k-slots.  Base ctor already
        // sized k_ for RK4 (4 slots); we keep the extras (they
        // cost almost nothing) and only use [0] and [1].
    }

    int  stages() const override { return 2; }

    // RK2 midpoint method evaluates derivatives at:
    //    k = 0 -> t0          (initial slope)
    //    k = 1 -> t0 + dt/2   (midpoint slope)
    double t_at_stage(int stage, double t0, double dt) const override {
        switch (stage) {
            case 0:  return t0;
            case 1:  return t0 + 0.5 * dt;
            default: return t0 + dt;
        }
    }

    void propagate(int stage, double dt) override {
        switch (stage) {
            case 0:
                // Save y_n and k1 = f(y_n, t_n).  Step to midpoint
                // for the next derivative evaluation.
                save_  = state_;
                k_[0]  = deriv_;
                state_ = save_ + 0.5 * dt * k_[0];
                break;
            case 1:
                // k2 = f(y*, t + dt/2).  Apply full step with k2.
                k_[1]  = deriv_;
                state_ = save_ + dt * k_[1];
                break;
            default:
                break;
        }
    }
};

#endif
