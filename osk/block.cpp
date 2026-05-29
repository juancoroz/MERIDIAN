
#include "block.h"

namespace osk {

Block::~Block() = default;

void Block::addIntegrator(double& state, double& deriv) {
    states_.push_back(std::unique_ptr<State>(new State(state, deriv)));
}

}
