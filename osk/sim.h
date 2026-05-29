
#ifndef OSK_SIM_H
#define OSK_SIM_H

#include "block.h"
#include <vector>

namespace osk {

class Sim {
public:

    Sim(const double*                                dts,
        double                                       tmax,
        const std::vector< std::vector<Block*> >&    stages);

    void run();

    static int stop;

private:
    std::vector<double>                       dts_;
    double                                    tmax_;
    std::vector< std::vector<Block*> >        stages_;

    void step(const std::vector<Block*>& stage, double dt,
              bool allow_abort = true);
};

}

#endif
