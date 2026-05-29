// =====================================================================
//  main.cpp  --  driver for the parameter-driven Section 4.1 example
// =====================================================================

#include "../../osk/osk.h"
#include "model.h"

using namespace osk;
using std::vector;

int main() {
    Filer ff("input.txt");
    ff.setLine0("Model");

    double tmax = ff.getDouble("tmax");
    double dt   = ff.getDouble("dt");

    Model *model = new Model(ff);

    vector<Block*>  vObj0  = { model };
    vector< vector<Block*> > vStage = { vObj0 };
    double dts[] = { dt };

    Sim *sim = new Sim(dts, tmax, vStage);
    sim->run();

    delete sim;
    delete model;
    return 0;
}
