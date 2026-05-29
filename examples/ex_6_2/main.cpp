// =====================================================================
//  main.cpp  --  Section 4.6.3.2 driver (PDF page 66)
//  Note: tmax = 1.0 but dt = 1.0 (coarse).  The model uses sample()
//  to actually force the integrator to step every 0.1 sec.
// =====================================================================

#include "../../osk/osk.h"
#include "model.h"

using namespace osk;
using std::vector;

int main() {
    double tmax = 1.0;
    double dt   = 1.0;
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
