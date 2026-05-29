
#include "state.h"
#include <cmath>

namespace osk {

double State::t         = 0.0;
double State::dt        = 0.0;
bool   State::tickfirst = false;
bool   State::ticklast  = false;
bool   State::stepstart = false;
std::vector<double> State::pending_events_;

State::State(double& state, double& deriv)
    : state_(state),
      deriv_(deriv),
      save_(0.0),
      k_(4, 0.0)
{
}

int State::stages() const { return 4; }

double State::t_at_stage(int stage, double t0, double dt) const {
    switch (stage) {
        case 0:  return t0;
        case 1:  return t0 + 0.5 * dt;
        case 2:  return t0 + 0.5 * dt;
        case 3:  return t0 + dt;
        default: return t0 + dt;
    }
}

void State::propagate(int stage, double dt) {
    switch (stage) {
        case 0:

            save_   = state_;
            k_[0]   = deriv_;
            state_  = save_ + 0.5 * dt * k_[0];
            break;
        case 1:

            k_[1]   = deriv_;
            state_  = save_ + 0.5 * dt * k_[1];
            break;
        case 2:

            k_[2]   = deriv_;
            state_  = save_ + dt * k_[2];
            break;
        case 3:

            k_[3]   = deriv_;
            state_  = save_ + (dt / 6.0)
                        * (k_[0] + 2.0 * k_[1] + 2.0 * k_[2] + k_[3]);
            break;
        default:
            break;
    }
}

bool State::sample(double period) {
    if (!stepstart)      return false;
    if (period <= 0.0)   return true;

    double r = std::fmod(t, period);
    if (r < 0.0) r += period;
    double d = std::min(r, period - r);

    const double tol = std::max(1.0e-9, 1.0e-9 * std::fabs(period));

    return d <= tol;
}

bool State::sample(EventTag, double t_event) {
    if (!stepstart) return false;

    const double tol = std::max(1.0e-9, 1.0e-9 * std::fabs(dt));

    if (std::fabs(t - t_event) <= tol) {

        return true;
    }

    if (t_event > t + tol) {
        pending_events_.push_back(t_event);
        double dt_to_event = t_event - t;
        if (dt_to_event < dt) dt = dt_to_event;
    }
    return false;
}

bool State::sample() {
    return stepstart;
}

}
