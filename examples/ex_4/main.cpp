// =====================================================================
//  main.cpp  --  Section 4.4 driver (PDF page 51)
//
//  Creates an autopilot + a missile, wires them with getsFrom() in
//  both directions, then runs them as a single stage with the
//  autopilot first in execution order.
// =====================================================================

#include "../../osk/osk.h"
#include "autopilot.h"
#include "missile.h"

using namespace osk;
using std::vector;

int main() {
    double tmax = 2.00;
    double dt   = 0.01;

    Autopilot *autopilot = new Autopilot(1000.0);
    Missile   *missile   = new Missile(0.0);

    // Wire the cross-references in both directions.
    autopilot->getsFrom(missile);
    missile->getsFrom(autopilot);

    // Single stage; autopilot executes before missile inside the stage.
    vector<Block*> vObj0;
    vObj0.push_back(autopilot);
    vObj0.push_back(missile);
    vector< vector<Block*> > vStage;
    vStage.push_back(vObj0);

    double dts[] = { dt };
    Sim *sim = new Sim(dts, tmax, vStage);
    sim->run();

    delete sim;
    delete autopilot;
    delete missile;
    return 0;
}
