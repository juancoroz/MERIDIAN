// =====================================================================
//  main.cpp  --  Section 4.1 driver (PDF page 41)
// =====================================================================

#include "../../osk/osk.h"
#include "model.h"

using namespace osk;
using std::vector;

int main() {
    double tmax = 2.00;
    double dt   = 0.01;
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
