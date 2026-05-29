
#include "sim.h"
#include <algorithm>
#include <limits>
#include <cmath>

namespace osk {

int Sim::stop = 0;

Sim::Sim(const double*                                dts,
         double                                       tmax,
         const std::vector< std::vector<Block*> >&    stages)
    : dts_(dts, dts + stages.size()),
      tmax_(tmax),
      stages_(stages)
{
}

void Sim::run() {

    State::t         = 0.0;
    State::ticklast  = false;
    Sim::stop        = 0;

    std::size_t istage = 0;
    while (istage < stages_.size()) {
        const auto& stage = stages_[istage];
        const double dt   = dts_[istage];
        State::dt         = dt;

        for (Block* b : stage) {
            if (b) {
                b->init();
                b->initCount++;
            }
        }

        State::tickfirst = true;

        bool first_step_of_stage = true;
        while (true) {

            if (State::t >= tmax_ - 1.0e-12) {
                State::ticklast  = true;
                State::stepstart = true;
                for (Block* b : stage) if (b) b->update();
                for (Block* b : stage) if (b) b->rpt();
                for (Block* b : stage) if (b) b->rpt();
                return;
            }

            double step_dt = dt;
            if (State::t + dt >= tmax_ - 1.0e-12) {
                step_dt = tmax_ - State::t;
            }

            step(stage, step_dt,  !first_step_of_stage);
            first_step_of_stage = false;

            if (Sim::stop < 0) {

                State::ticklast  = true;
                State::stepstart = true;
                for (Block* b : stage) if (b) b->update();
                for (Block* b : stage) if (b) b->rpt();
                for (Block* b : stage) if (b) b->rpt();
                return;
            }
            if (Sim::stop > 0) {

                std::size_t next = static_cast<std::size_t>(Sim::stop);
                Sim::stop = 0;
                if (next >= stages_.size()) return;
                istage = next;
                break;
            }

        }
    }
}

void Sim::step(const std::vector<Block*>& stage, double dt,
               bool allow_abort) {

    int nstages = 4;
    for (Block* b : stage) {
        if (b && !b->states_.empty()) {
            nstages = b->states_[0]->stages();
            break;
        }
    }

    const double t0 = State::t;

    State::pending_events_.clear();

    State::t         = t0;
    State::stepstart = true;

    State::dt = dt;

    for (Block* b : stage) if (b) b->update();

    const double eps = 1.0e-9 * std::max(1.0, std::fabs(t0));
    double t_grid;
    {

        double k = std::floor(t0 / dt + 1.0 - 1.0e-12);
        t_grid = k * dt;

        while (t_grid <= t0 + eps) t_grid += dt;
    }

    double actual_dt = t_grid - t0;
    if (!State::pending_events_.empty()) {

        double t_next = std::numeric_limits<double>::infinity();
        for (double te : State::pending_events_) {
            if (te > t0 + eps && te < t_next) t_next = te;
        }
        if (t_next != std::numeric_limits<double>::infinity()) {
            double dt_to_event = t_next - t0;
            if (dt_to_event < actual_dt) actual_dt = dt_to_event;
        }
    }

    State::dt = actual_dt;

    for (Block* b : stage) if (b) b->rpt();

    if (allow_abort && Sim::stop != 0) {
        return;
    }

    Sim::stop = 0;

    for (Block* b : stage) {
        if (!b) continue;
        for (auto& s : b->states_) s->propagate(0, actual_dt);
    }
    State::tickfirst = false;

    State* time_ref = nullptr;
    for (Block* b : stage) {
        if (b && !b->states_.empty()) {
            time_ref = b->states_[0].get();
            break;
        }
    }
    for (int k = 1; k < nstages; ++k) {
        State::stepstart = false;
        if (time_ref) {
            State::t = time_ref->t_at_stage(k, t0, actual_dt);
        } else {
            State::t = t0 + actual_dt;
        }
        for (Block* b : stage) if (b) b->update();
        for (Block* b : stage) {
            if (!b) continue;
            for (auto& s : b->states_) s->propagate(k, actual_dt);
        }
    }

    State::t = t0 + actual_dt;
}

}
