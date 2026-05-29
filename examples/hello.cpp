// =====================================================================
//  hello.cpp  --  the Section 2 missile / target / radar example
//
//  Exercises:
//    - Block derivation
//    - addIntegrator() for constant-derivative (target) and updated-
//      derivative (missile pursuit, radar 2nd-order servo) models
//    - State::sample() / State::tickfirst / State::ticklast for
//      throttled output
//    - Sim::stop = +1 to advance from stage 0 (target+radar) to
//      stage 1 (target+missile+radar) when the radar locks on
//    - Sim::stop = -1 to terminate at intercept
// =====================================================================

#include "../osk/osk.h"
#include <cmath>
#include <cstdio>

using namespace osk;

static const double R  = 180.0 / PI;

// ---------------------------------------------------------------------
//  Target: constant velocity in two channels.
// ---------------------------------------------------------------------
class Target : public Block {
public:
    double x, y, vx, vy;

    Target(double x0, double y0, double vx0, double vy0) {
        x = x0; y = y0; vx = vx0; vy = vy0;
        addIntegrator(x, vx);
        addIntegrator(y, vy);
    }
    void rpt() override {
        if (State::sample(1.0) || State::tickfirst || State::ticklast) {
            std::printf("%12s %8.3f %8.3f %8.3f\n",
                        "Target", State::t, x, y);
        }
    }
};

// ---------------------------------------------------------------------
//  Missile: pursues the target at constant speed.
// ---------------------------------------------------------------------
class Missile : public Block {
public:
    double x, y, vx, vy, vel, d;
    Target* target;

    Missile(Target* tgt, double x0, double y0, double v) {
        target = tgt;
        x = x0; y = y0; vel = v;
        vx = 0.0; vy = 0.0; d = 0.0;
        addIntegrator(x, vx);
        addIntegrator(y, vy);
    }
    void update() override {
        double dx = target->x - x;
        double dy = target->y - y;
        d  = std::sqrt(dx * dx + dy * dy);
        if (d > 0.0) {
            vx = vel * dx / d;
            vy = vel * dy / d;
        }
        if (d <= 0.1) Sim::stop = -1;     // intercept -> terminate
    }
    void rpt() override {
        if (State::sample(1.0) || State::tickfirst || State::ticklast) {
            std::printf("%12s %8.3f %8.3f %8.3f %8.3f\n",
                        "Missile", State::t, x, y, d);
        }
    }
};

// ---------------------------------------------------------------------
//  Radar: second-order servo pointing at the target.  When its
//  pointing error drops below 1 degree it fires Sim::stop = 1 to
//  hand control over to stage 1 (which adds the missile).
// ---------------------------------------------------------------------
class Radar : public Block {
public:
    double x1, x1d, x2, x2d;
    double theta_err, wn, zeta;
    Target* target;

    Radar(Target* tgt, double theta_deg, double wn_, double zeta_) {
        target = tgt;
        x1 = theta_deg / R;
        x2 = 0.0;
        x1d = 0.0; x2d = 0.0;
        theta_err = 0.0;
        wn = wn_; zeta = zeta_;
        addIntegrator(x1, x1d);
        addIntegrator(x2, x2d);
    }
    void update() override {
        double theta_target = std::atan(target->y / target->x);
        theta_err = theta_target - x1;
        if (std::fabs(theta_err * R) < 1.0 && Sim::stop == 0) {
            Sim::stop = 1;   // launch the missile -> stage 1
        }
        x1d = x2;
        x2d = theta_err * wn * wn - 2.0 * zeta * wn * x2;
    }
    void rpt() override {
        if (State::sample(1.0) || State::tickfirst || State::ticklast) {
            std::printf("%12s %8.3f %8.3f\n",
                        "Radar", State::t, theta_err * R);
            std::printf("\n");
        }
    }
};

// ---------------------------------------------------------------------
int main() {
    double tmax = 10.0;
    double dt   = 0.01;

    Target*  target  = new Target(20.0, 5.0, -1.0, 0.0);
    Missile* missile = new Missile(target, 0.0, 0.0, 2.0);
    Radar*   radar   = new Radar(target, 60.0, 2.643, 0.7);

    std::vector<Block*> stage0 = { target, radar };
    std::vector<Block*> stage1 = { target, missile, radar };
    std::vector< std::vector<Block*> > stages = { stage0, stage1 };

    double dts[] = { dt, dt };
    Sim sim(dts, tmax, stages);
    sim.run();

    delete target;
    delete missile;
    delete radar;
    return 0;
}
