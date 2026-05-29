// =====================================================================
//  main.cpp  --  Section 4.5 driver (PDF page 57)
//
//  Two-stage configuration:
//    Stage 0: missile alone, flies open-loop (gammad = 0)
//    Stage 1: autopilot + missile, missile flies under closed-loop
//             control of the autopilot (gammad = a_cmd / v)
//
//  The missile object lives in BOTH stages -- it doesn't get
//  re-created when stage 1 begins; it just keeps propagating from
//  whatever state it had at the end of stage 0.  The staging event
//  is triggered by the missile itself at t = 1.0 s.
// =====================================================================

#include "../../osk/osk.h"
#include "autopilot.h"
#include "missile.h"

using namespace osk;
using std::vector;

int main() {
    double tmax = 3.00;
    double dt   = 0.01;

    Autopilot *autopilot = new Autopilot(1000.0);
    Missile   *missile   = new Missile(0.0);

    autopilot->getsFrom(missile);
    missile->getsFrom(autopilot);

    // Stage 0: missile only.
    vector<Block*> vObj0;
    vObj0.push_back(missile);

    // Stage 1: autopilot first, then missile.
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
