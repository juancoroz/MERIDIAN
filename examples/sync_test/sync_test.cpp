// =====================================================================
//  sync_test.cpp  --  end-to-end verification of Section 4.2.2 claim
//
//  Section 4.2.2 argues that putting State::sample() blocks anywhere
//  inside a single update() method lets the modeler precisely control
//  the order in which discrete and continuous operations execute,
//  avoiding the "artificial lag" that two-method architectures
//  (Figure 12) impose.
//
//  This test verifies the claim by implementing the two systems
//  shown on PDF page 47:
//
//    Figure 11:  u -> [sample dt] -> w -> [K] -> xdot -> [INT] -> y=x
//                The discrete sampler runs BEFORE the continuous
//                derivative, so xdot must read the freshly-sampled w.
//
//    Figure 13:  u -> xdot -> [INT] -> x -> [K] -> w -> [sample] -> y
//                The continuous derivative runs BEFORE the discrete
//                sampler, so the sampler must read the freshly-
//                computed w.
//
//  For each system we run TWO models:
//    *_correct  -- the CMD-recommended source ordering
//    *_lagged   -- the opposite ordering (the bad case 4.2.2 warns
//                  against)
//
//  The two should differ by exactly one sample period.  If my kernel
//  respects source-code order inside update() -- as it must to support
//  the 4.2.2 architecture -- the correct/lagged outputs will diverge
//  by one sample at every sample boundary, with the kernel showing
//  zero ambiguity about which version is which.
// =====================================================================

#include "../../osk/osk.h"
#include <cstdio>
#include <cmath>

using namespace osk;

// Sampling period for both systems.
static const double DT_SAMPLE = 0.1;

// =====================================================================
//  Figure 11: discrete-then-continuous
// =====================================================================

// u(t) = t  (a ramp input)
static double u_of(double t) { return t; }

// Correct ordering: sample w from u FIRST, then use w in xdot.
class Fig11_correct : public Block {
public:
    double x, xdot;
    double w;
    double K;
    Fig11_correct() {
        x = 0.0; xdot = 0.0;
        w = 0.0;
        K = 1.0;
        addIntegrator(x, xdot);
    }
    void update() override {
        if (State::sample(DT_SAMPLE)) {
            w = u_of(State::t);          // discrete update FIRST
        }
        xdot = K * w;                    // continuous reads fresh w
    }
};

// Lagged ordering: compute xdot from the OLD w, then sample.  The
// derivative is one sample behind the input.
class Fig11_lagged : public Block {
public:
    double x, xdot;
    double w;
    double K;
    Fig11_lagged() {
        x = 0.0; xdot = 0.0;
        w = 0.0;
        K = 1.0;
        addIntegrator(x, xdot);
    }
    void update() override {
        xdot = K * w;                    // continuous reads STALE w
        if (State::sample(DT_SAMPLE)) {
            w = u_of(State::t);          // discrete update LAST
        }
    }
};

// =====================================================================
//  Figure 13: continuous-then-discrete
// =====================================================================

// u(t) = 1  (constant unity input -> x = t)
static double u13() { return 1.0; }

// Correct ordering: update xdot and w FIRST, then sample y from w.
class Fig13_correct : public Block {
public:
    double x, xdot;
    double w, y;
    double K;
    Fig13_correct() {
        x = 0.0; xdot = 0.0;
        w = 0.0; y = 0.0;
        K = 1.0;
        addIntegrator(x, xdot);
    }
    void update() override {
        xdot = u13();                    // continuous first
        w    = K * x;
        if (State::sample(DT_SAMPLE)) {
            y = w;                       // sampler reads fresh w
        }
    }
};

// Lagged ordering: sample y from w FIRST (using the PREVIOUS step's
// value of w), then update xdot and w.  The output y is one sample
// behind the integrated state.
class Fig13_lagged : public Block {
public:
    double x, xdot;
    double w, y;
    double K;
    Fig13_lagged() {
        x = 0.0; xdot = 0.0;
        w = 0.0; y = 0.0;
        K = 1.0;
        addIntegrator(x, xdot);
    }
    void update() override {
        if (State::sample(DT_SAMPLE)) {
            y = w;                       // sampler reads STALE w
        }
        xdot = u13();                    // continuous AFTER sampler
        w    = K * x;
    }
};

// =====================================================================
//  Helper: run a Block for tmax seconds and return the value of a
//  named scalar at the end, plus a copy of its value just after t=1.
//  We pass references so the runner doesn't have to know the model's
//  internals.
// =====================================================================
template <class M>
void run(M* m, double tmax, double dt) {
    Sim::stop = 0;
    State::t  = 0.0;
    std::vector<Block*> stage = { m };
    std::vector< std::vector<Block*> > stages = { stage };
    double dts[] = { dt };
    Sim sim(dts, tmax, stages);
    sim.run();
}

// =====================================================================
//  Test driver
// =====================================================================
int main() {
    const double tmax = 2.0;
    const double dt   = 0.01;

    std::printf("=== Figure 11: u -> [sample] -> w -> K -> xdot -> INT -> x ===\n");
    std::printf("Input u(t)=t, K=1, sample period=%.2f s.\n", DT_SAMPLE);
    std::printf("Analytical for CORRECT ordering: w is a staircase of u\n");
    std::printf("  starting at u(0)=0, so x is the integral of that staircase.\n");
    std::printf("  At t=2.0:  w=u(2.0)=2.0, and x = sum over n=0..19 of\n");
    std::printf("             0.1 * u(n*0.1) = 0.1 * (0+0.1+...+1.9) = 0.1*19 = 1.9\n");
    std::printf("Analytical for LAGGED ordering: w stays at 0 for the first\n");
    std::printf("  0.1 s (because xdot reads w BEFORE the t=0 sample updates w),\n");
    std::printf("  so x is shifted by one sample period.  Final x = 1.8.\n");
    std::printf("\n");

    {
        Fig11_correct* m = new Fig11_correct();
        run(m, tmax, dt);
        std::printf("  CORRECT: x(tmax) = %.6f  w(tmax) = %.6f\n", m->x, m->w);
        delete m;
    }
    {
        Fig11_lagged* m = new Fig11_lagged();
        run(m, tmax, dt);
        std::printf("  LAGGED : x(tmax) = %.6f  w(tmax) = %.6f\n", m->x, m->w);
        delete m;
    }

    std::printf("\n=== Figure 13: u -> xdot -> INT -> x -> K -> w -> [sample] -> y ===\n");
    std::printf("Input u(t)=1, K=1, sample period=%.2f s.\n", DT_SAMPLE);
    std::printf("So x = t, w = t, and y should sample w at the boundaries.\n");
    std::printf("Analytical for CORRECT ordering: y(tmax) = w(tmax) = 2.0\n");
    std::printf("Analytical for LAGGED ordering: y is the PREVIOUS w sample,\n");
    std::printf("  so at t=tmax we have y = w(tmax - dt_sample) = 1.9\n");
    std::printf("\n");

    {
        Fig13_correct* m = new Fig13_correct();
        run(m, tmax, dt);
        std::printf("  CORRECT: x(tmax) = %.6f  y(tmax) = %.6f\n", m->x, m->y);
        delete m;
    }
    {
        Fig13_lagged* m = new Fig13_lagged();
        run(m, tmax, dt);
        std::printf("  LAGGED : x(tmax) = %.6f  y(tmax) = %.6f\n", m->x, m->y);
        delete m;
    }

    return 0;
}
