//  peak_tracker.cpp

#include "peak_tracker.h"
#include "environment.h"
#include "kinematics.h"

namespace rocket6dof {

void PeakTracker::update() {
    if (!env || !kin) return;

    double q = std::fabs(env->pdynmc);
    double a = std::fabs(kin->alphax);

    // Dynamic pressure: peak during ascent phase only (t < 200s).
    // Sub-orbital trajectories that re-enter can hit higher q_dyn
    // during descent than during their ascent max-Q -- e.g. a
    // vehicle that fails to reach orbit and comes back through
    // dense atmosphere at high Mach.  For load-envelope analysis
    // we want the structural q seen on the way up, not on the way
    // down (which is also informative but a different question).
    //
    // Record alpha at peak q simultaneously, so users can read off
    // the q*alpha at the structurally-critical moment.
    if (q > max_q_dyn && osk::State::t < t_ascent_max) {
        max_q_dyn      = q;
        t_max_q_dyn    = osk::State::t;
        alpha_at_max_q = a;
    }

    // q * alpha is the structural load proxy (bending moment scales
    // with q*alpha for a slender body).  Track this independently;
    // it may peak slightly off from max-Q if alpha is rising fast.
    // Gate to ascent phase (t < 200s) for the same reason as
    // max_alpha below.
    double qa = q * a;
    if (qa > max_q_alpha && osk::State::t < t_ascent_max) {
        max_q_alpha   = qa;
        t_max_q_alpha = osk::State::t;
    }

    // Max alpha: gate on (q_dyn > 1 kPa) AND (t < 200 s) to capture
    // only the ascent phase where structural alpha matters.  After
    // stage 1 burnout (~155s) the vehicle is in vacuum or near-
    // vacuum and bending loads are negligible.  This excludes the
    // re-entry pathology where a failed-orbit trajectory comes back
    // through 1 kPa with alpha ~ 180 deg.
    //
    // Picking t<200s rather than tying to stage transitions keeps
    // this simple.  A config where stage 1 burns much longer (e.g.
    // 500s) would need a different threshold; that's out of scope
    // for the current launcher class but could be exposed via a
    // config field later.
    if (q > 1000.0 && osk::State::t < t_ascent_max && a > max_alpha) {
        max_alpha   = a;
        t_max_alpha = osk::State::t;
    }

    // Max Mach: track during ascent (same gate as the other peaks).
    // Useful for thermal/aero envelope analysis.  Unlike alpha/q,
    // max Mach is well-defined at all velocities (no atan2 issues),
    // so no q_dyn floor needed.
    double m = env->vmach;
    if (m > max_mach && osk::State::t < t_ascent_max) {
        max_mach   = m;
        t_max_mach = osk::State::t;
    }
}

} // namespace rocket6dof
