//  peak_tracker.h  --  Per-trajectory peak monitoring
//
//  Tracks running maxima of selected quantities during a sim run so
//  they can be reported as outputs at the end.  The motivating use
//  case is Monte Carlo studies of launch-vehicle load envelope:
//  max(q_dyn) and max(|alpha|) per run characterize the structural
//  load the vehicle saw, neither of which is captured by reading the
//  end-of-run state alone.
//
//  Current tracked quantities:
//    * max_q_dyn   = max of |env->pdynmc|       [Pa]
//    * max_alpha   = max of |kin->alphax|       [deg]
//
//  Both initialize to 0 and update every step.  No reset semantics;
//  the block lives for one sim run.
//
//  Wiring:
//    PeakTracker* pk = new PeakTracker();
//    pk->getsFrom(env, kin);
//    // include pk in the stage block list, anywhere after env and kin
//
//  Outputs:
//    Add "max_q_dyn" or "max_alpha" to monte_carlo.outputs in the
//    config.  read_output() handles them via the Vehicle's pk pointer.

#ifndef ROCKET6DOF_PEAK_TRACKER_H
#define ROCKET6DOF_PEAK_TRACKER_H

#include "../osk/osk.h"
#include <cmath>

namespace rocket6dof {

class Environment;
class Kinematics;

class PeakTracker : public osk::Block {
public:
    Environment* env;
    Kinematics*  kin;
    void getsFrom(Environment* e, Kinematics* k) { env = e; kin = k; }

    // Cutoff time after which peaks are no longer updated.  Default
    // 200s covers stage 1 + a margin for the launcher class.  Adjust
    // for vehicles with longer or shorter ascent profiles.  This
    // excludes re-entry phenomena from sub-orbital failed-orbit
    // trajectories (where q_dyn can peak higher on descent than on
    // ascent, and alpha can be ~180 deg, neither of which is the
    // structural load envelope we usually want to characterize).
    double t_ascent_max;

    // Tracked peaks.  All initialized to 0 in init().
    double max_q_dyn;    // [Pa]  max |env->pdynmc| seen this run
    double max_alpha;    // [deg] max |kin->alphax| during ascent
    double max_q_alpha;  // [Pa*deg] max |q_dyn * alpha| -- bending-load proxy
    double max_mach;     // [-]   max env->vmach during ascent

    // Time at which each peak occurred (useful for understanding
    // when in the trajectory the load envelope was hit).
    double t_max_q_dyn;
    double t_max_alpha;
    double t_max_q_alpha;
    double t_max_mach;
    double alpha_at_max_q;   // [deg] alpha at the time q_dyn peaked

    PeakTracker() : env(nullptr), kin(nullptr),
                    t_ascent_max(200.0),
                    max_q_dyn(0.0), max_alpha(0.0), max_q_alpha(0.0),
                    max_mach(0.0),
                    t_max_q_dyn(0.0), t_max_alpha(0.0), t_max_q_alpha(0.0),
                    t_max_mach(0.0),
                    alpha_at_max_q(0.0) {}

    void init() override {
        max_q_dyn   = 0.0;
        max_alpha   = 0.0;
        max_q_alpha = 0.0;
        max_mach    = 0.0;
        t_max_q_dyn   = 0.0;
        t_max_alpha   = 0.0;
        t_max_q_alpha = 0.0;
        t_max_mach    = 0.0;
        alpha_at_max_q = 0.0;
    }

    void update() override;
    void rpt() override {}
};

} // namespace rocket6dof

#endif
