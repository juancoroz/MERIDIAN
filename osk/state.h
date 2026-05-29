
#ifndef OSK_STATE_H
#define OSK_STATE_H

#include <vector>

namespace osk {

class State {
public:

    double& state_;
    double& deriv_;

    double             save_;
    std::vector<double> k_;

    State(double& state, double& deriv);
    virtual ~State() = default;

    virtual int    stages() const;
    virtual double t_at_stage(int stage, double t0, double dt) const;
    virtual void   propagate(int stage, double dt);

    static double t;
    static double dt;
    static bool   tickfirst;
    static bool   ticklast;

    static bool sample(double period);

    enum EventTag { EVENT };

    static bool sample(EventTag, double t_event);

    static bool sample();

    static bool stepstart;

    static std::vector<double> pending_events_;
};

}

#endif
