// =====================================================================
//  main.cpp  --  Appendix 2 test driver (PDF page A-20)
//
//  tmax = 4.0, dt = 0.001, RK2 integrator.  Expected output:
//    starting Model...
//    0.000 1.000000
//    0.800 0.916113
//    1.600 0.405399
//    ...
//    4.000 0.219972
// =====================================================================

#include "../../osk/osk.h"
#include "model.h"

using namespace osk;
using std::vector;

int main() {
    double tmax = 4.00;
    double dt   = 0.001;

    Model *model = new Model(1.0, 0.0);

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
