// =====================================================================
//  main.cpp  --  Section 4.6.3.3 driver
//
//  Nominal dt=0.1 sec is used outside the [0.55, 0.60] window.  Inside
//  that window the model forces dt=0.01 sec by re-scheduling rolling
//  EVENT samples.
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
