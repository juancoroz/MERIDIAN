// =====================================================================
//  main.cpp  --  Section 4.6.3.4 driver (PDF page 73)
//
//  Same multi-stage train-of-objects as Section 4.5 but with the
//  nominal dt set to 0.1 (modification #1).  The autopilot uses the
//  async scheduler to force dt=0.01 inside [1.55, 1.63] and to capture
//  one-time events at 1.6001 / 1.6002.
// =====================================================================

#include "../../osk/osk.h"
#include "autopilot.h"
#include "missile.h"

using namespace osk;
using std::vector;

int main() {
    double tmax = 2.00;
    double dt   = 0.1;

    Autopilot *autopilot = new Autopilot(1000.0);
    Missile   *missile   = new Missile(0.0);

    autopilot->getsFrom(missile);
    missile->getsFrom(autopilot);

    vector<Block*> vObj0;
    vObj0.push_back(missile);

    vector<Block*> vObj1;
    vObj1.push_back(autopilot);
    vObj1.push_back(missile);

    vector< vector<Block*> > vStage;
    vStage.push_back(vObj0);
    vStage.push_back(vObj1);

    double dts[] = { dt, dt };

    Sim *sim = new Sim(dts, tmax, vStage);
    sim->run();

    delete sim;
    delete autopilot;
    delete missile;
    return 0;
}
