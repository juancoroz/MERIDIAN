// =====================================================================
//  shm.cpp  --  simple harmonic motion sanity check
//
//  Integrates  x'' + x = 0   with  x(0)=1, x'(0)=0.
//  Analytical solution: x(t) = cos(t), v(t) = -sin(t).
//
//  Verifies the RK4 integrator hits the analytical answer to within
//  the expected ~ dt^4 error.
// =====================================================================

#include "../osk/osk.h"
#include <cmath>
#include <cstdio>

using namespace osk;

class SHM : public Block {
public:
    double x, v, xd, vd;
    SHM() {
        x = 1.0; v = 0.0;
        xd = 0.0; vd = 0.0;
        addIntegrator(x, xd);
        addIntegrator(v, vd);
    }
    void update() override {
        xd =  v;
        vd = -x;
    }
    void rpt() override {
        if (State::sample(1.0) || State::ticklast) {
            double x_true = std::cos(State::t);
            std::printf("t=%6.2f   x=%9.6f   x_true=%9.6f   err=%10.2e\n",
                        State::t, x, x_true, x - x_true);
        }
    }
};

int main() {
    SHM* m = new SHM();
    std::vector<Block*> stage = { m };
    std::vector< std::vector<Block*> > stages = { stage };
    double dts[] = { 0.01 };
    Sim sim(dts, 10.0, stages);
    sim.run();
    delete m;
    return 0;
}
