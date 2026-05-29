
#ifndef OSK_BLOCK_H
#define OSK_BLOCK_H

#include "state.h"
#include <vector>
#include <memory>

namespace osk {

class Sim;

class Block {
public:
    Block() : initCount(0) {}
    virtual ~Block();

    int initCount;

    virtual void init()   {}
    virtual void update() {}
    virtual void rpt()    {}

    void addIntegrator(double& state, double& deriv);

    template <class StateT>
    void addIntegrator(double& state, double& deriv) {
        states_.push_back(std::unique_ptr<State>(new StateT(state, deriv)));
    }

    std::vector<std::unique_ptr<State>> states_;
};

}

#define ACCESS_FN(type, var)                              \
    public:                                               \
        type var##_ () const { return this->var; }        \
    public:

#endif
