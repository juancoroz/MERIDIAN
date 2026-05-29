// =====================================================================
//  main.cpp  --  Section 4.6.3.1 driver
//
//  tmax = 1.0, dt = 0.1.  The three EVENT samples in model.cpp's
//  update() fall at 0.15, 0.185, 0.1855 -- all at finer resolution
//  than the integrating step.  The OSK's asynchronous scheduler is
//  what makes them fire exactly.
// =====================================================================

#include "../../osk/osk.h"
#include "model.h"

using namespace osk;
using std::vector;

int main() {
    double tmax = 1.0;
    double dt   = 0.1;
    Model *model = new Model(0.0);
    vector<Block*> vObj0;
    vObj0.push_back(model);
    vector< vector<Block*> > vStage;
    vStage.push_back(vObj0);
    double dts[] = { dt };
    Sim *sim = new Sim(dts, tmax, vStage);
    sim->run();

    delete sim;
    delete model;
    return 0;
}
