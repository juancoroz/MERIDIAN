#ifndef EX1IO_MODEL_H
#define EX1IO_MODEL_H

#include "../../osk/osk.h"

class Model : public osk::Block {
public:
    double gamma, gammad;
    double gamma_cmd;
    double a_cmd, k, v;
    double gamma0;

    explicit Model(osk::Filer& ff);

    void init()   override;
    void update() override;
    void rpt()    override;
};

#endif
