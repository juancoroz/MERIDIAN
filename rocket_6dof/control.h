//  control.h  --  Two-channel acceleration autopilot
//
//  Public members match the "out"-tagged module variables of Zipfel,
//  "Modeling and Simulation of Aerospace Vehicle Dynamics" 3rd ed.,
//  Section 10.3.x (hyper[500-599]).  Implements the pole-placement
//  acceleration controller from Zipfel's control.cpp.
//
//  Two independent channels:
//
//    Pitch (mautp = 3):  takes commanded normal accel ancomx (in g's),
//                        outputs commanded pitch deflection delecx (deg).
//                        Feedback uses sensed body pitch rate (qqcx)
//                        and sensed body z-accel (FSPCB.z).
//
//    Yaw (mauty = 5):    takes commanded lateral accel alcomx (in g's),
//                        outputs commanded yaw deflection delrcx (deg).
//                        Feedback uses sensed body yaw rate (rrcx) and
//                        sensed body y-accel (FSPCB.y).
//
//  Each channel has ONE integrator state (zz for pitch, yy for yaw)
//  that accumulates the error signal.  Three feedback gains are
//  computed online based on:
//      * closed-loop pole placement parameters (waclp, zaclp, paclp)
//      * dimensional derivatives from Aero (dla, dma, dmq, dmde for pitch)
//
//  The dominant closed-loop poles are GAIN-SCHEDULED with dynamic
//  pressure (Zipfel lines 184-185, 265-266):
//      waclp = (0.1 + 0.5e-5*(q - 20e3)) * (1 + factwaclp)
//      paclp =  0.7 + 1e-5 *(q - 20e3)  * (1 + factwaclp)
//
//  Sensor source:
//      qqcx, rrcx -- body rates from Euler::WBEB
//      FSPCB      -- specific force from Newton::FSPB
//      pdynmc     -- dynamic pressure from Environment
//      thrust on  -- from Propulsion::mprop  (control only fires when
//                    rocket is thrusting; coast-phase has no TVC authority)

#ifndef ROCKET6DOF_CONTROL_H
#define ROCKET6DOF_CONTROL_H

#include "../osk/osk.h"

namespace rocket6dof {

class Environment;
class Newton;
class Euler;
class Kinematics;
class Aerodynamics;
class Propulsion;
class INS;

class Control : public osk::Block {
public:
    // ---- Inputs ----
    Environment*  env;
    Newton*       newton;
    Euler*        euler;
    Kinematics*   kin;        // for body attitude (used in cascaded modes)
    Aerodynamics* aero;
    Propulsion*   prop;
    INS*          ins;        // optional; when wired, sensor reads come from INS

    // Thrust attachment point from vehicle nose [m].  Used to compute
    // gimbal moment arm as (thrust_loc - xcg).  Default 5.0 reflects
    // the small-vehicle assumption inherited from Zipfel's reference
    // configurations; large vehicles (e.g. mission_launcher.json)
    // should override this via config.  Must match tvc.parm if TVC
    // is being used, since TVC also reads its own parm value.
    double thrust_loc = 5.0;

    // ---- Load relief (maut=60 only) ----
    // Reduce structural q*alpha load during atmospheric ascent by
    // blending theta_com_eff toward the velocity-vector pitch (FPA)
    // when dynamic pressure is high.  At full blend the body tracks
    // velocity (alpha=0) and produces no aerodynamic side-load; at
    // zero blend the autopilot follows the open-loop attitude
    // program unchanged.
    //
    // Blend factor uses a smooth sigmoid in q_dyn:
    //   q < q_threshold - q_width:  w = 0   (no relief)
    //   q = q_threshold:             w = 0.5
    //   q > q_threshold + q_width:  w = 1   (full relief)
    //   theta_com_eff = (1-w)*theta_com_inertial + w*fpa
    //
    // Default load_relief_q_threshold = 0 disables the feature
    // (existing configs unaffected).  Set to ~15000 Pa to engage
    // relief through max-Q (which peaks around 27 kPa for the
    // launcher).
    //
    // KNOWN LIMITATION: this velocity-frame-blending implementation
    // is NOT robust under dispersions (~5% of MC runs diverge).
    // See accel_relief_gain below for the production-grade
    // accelerometer-feedback version.
    double load_relief_q_threshold = 0.0;   // [Pa]  sigmoid center; 0 = disabled
    double load_relief_q_width     = 5000.0; // [Pa]  half-width of transition

    // ---- Accelerometer-feedback load relief (maut=60 only) ----
    // Production-style load relief: subtract a gimbal-deflection
    // term proportional to the HIGH-PASS-FILTERED body lateral
    // acceleration.  Real launchers use accelerometer feedback
    // with a washout filter for exactly this reason: the autopilot
    // needs nonzero steady-state lateral force to shape the
    // gravity-turn trajectory, so the feedback must reject only
    // TRANSIENT deviations (gusts, wind shear, q*alpha excursions),
    // not the DC component.
    //
    // Mechanism (inside the maut=50||60 branch before saturation):
    //   anx_lp <- anx_lp + alpha * (anx - anx_lp)   # low-pass state
    //   anx_hp = anx - anx_lp                        # high-pass output
    //   delecx <- delecx - w(q) * accel_relief_gain * anx_hp
    //   delrcx <- delrcx - w(q) * accel_relief_gain * ayx_hp
    // where alpha = dt / (tau + dt) and w(q) is a sigmoid in q_dyn.
    //
    // Time constant `accel_relief_tau` separates the load-relief
    // band from the trajectory-shaping band.  ~1 s isolates
    // sub-second transients while letting longer-timescale aero
    // loads pass to the autopilot.
    //
    // Default accel_relief_gain = 0 disables the feature.
    // Recommended starting values for a vehicle with meaningful
    // lateral loads:
    //   accel_relief_gain = 1.0 to 5.0    [deg/g]
    //   accel_relief_tau  = 1.0           [s]
    //   accel_relief_q_threshold = 10000  [Pa]
    //   accel_relief_q_width     = 3000   [Pa]
    // Test on heavy_launcher.json (peak anx ~0.10 g) before
    // applying to a new vehicle.
    double accel_relief_gain        = 0.0;     // [deg/g] gimbal per g of HP-filtered lateral accel; 0 = disabled
    double accel_relief_q_threshold = 10000.0; // [Pa]  q above which feedback engages
    double accel_relief_q_width     = 3000.0;  // [Pa]  half-width of engagement sigmoid
    double accel_relief_tau         = 1.0;     // [s]   washout time constant; transients faster than this are rejected

    // ---- Load-relief autopilot detuning (active when accel feedback engages) ----
    // When accel feedback is providing load relief, the outer attitude
    // loop should be SLOWED so the feedback can do its job without
    // fighting attitude tracking.  Multiplier on tau_att applied
    // through the same sigmoid in q_dyn used for the feedback term.
    //
    // Why this matters: the attitude autopilot in maut=60 uses
    // tau_att (outer loop) as its dominant time constant -- it
    // drives theta -> theta_com over ~tau_att seconds.  Load
    // relief, in contrast, accepts a brief attitude error in
    // exchange for reduced bending moment.  If tau_att is small
    // while load relief is also active, the autopilot will fight
    // the feedback term, pushing the closed loop toward instability
    // under dispersion (verified on heavy_launcher MC).
    //
    // The fix is a *blended* time constant: full tau_att at low q
    // (attitude tracking takes priority), increased tau_att at
    // high q where load relief takes priority.  Same q_dyn sigmoid
    // as the feedback term means the two are perfectly coordinated.
    //
    // Default 1.0 = no scaling = backward compatible.
    // Recommended 2.0-3.0 for vehicles where load relief is
    // operationally important (heavy_launcher class).
    //
    // Heavy-vehicle MC validation (n=100) -- see SESSION_NOTES.
    double load_relief_tau_factor   = 1.0;     // [-]   tau_att multiplier when load-relief sigmoid is active; 1=disabled

    // Internal state for the washout filter.  Not user-configurable;
    // initialized to 0 in init() and updated each step.
    double anx_lp_state             = 0.0;
    double ayx_lp_state             = 0.0;
    void getsFrom(Environment* e, Newton* n, Euler* eu,
                  Aerodynamics* a, Propulsion* p) {
        env = e; newton = n; euler = eu; aero = a; prop = p;
        kin = nullptr; ins = nullptr;
    }
    void getsFrom(Environment* e, Newton* n, Euler* eu,
                  Aerodynamics* a, Propulsion* p, INS* i) {
        env = e; newton = n; euler = eu; aero = a; prop = p;
        kin = nullptr; ins = i;
    }
    void getsFrom(Environment* e, Newton* n, Euler* eu, Kinematics* k,
                  Aerodynamics* a, Propulsion* p) {
        env = e; newton = n; euler = eu; kin = k; aero = a; prop = p;
        ins = nullptr;
    }
    void getsFrom(Environment* e, Newton* n, Euler* eu, Kinematics* k,
                  Aerodynamics* a, Propulsion* p, INS* i) {
        env = e; newton = n; euler = eu; kin = k; aero = a; prop = p;
        ins = i;
    }

    // ---- Mode selector (Zipfel maut = |mauty|mautp|) ----
    //   mautp = 0 -> no pitch control;  mautp = 3 -> pitch accel control
    //   mauty = 0 -> no yaw control;    mauty = 5 -> yaw accel control
    // Encoded as maut = mauty*10 + mautp (so 53 = both on).
    //
    //   maut = 99: VACUUM-MODE rate damper.  Pure rate-feedback; ignores
    //              ancomx/alcomx.  Robust in both atmosphere and vacuum
    //              but does no active attitude tracking.
    //
    //   maut = 60: CASCADED ATTITUDE-TRACKING autopilot.  Outer loop
    //              drives body attitude to (theta_com_inertial,
    //              psi_com_inertial).  Inner loop is the same rate
    //              damper as maut=99.  Works in atmosphere and vacuum.
    //
    //   maut = 50: CASCADED ACCELERATION-TRACKING autopilot.  Like
    //              maut=60 but takes ancomx/alcomx as input and
    //              internally converts to a target attitude using
    //              current aero (q*S*CN_a) to compute the alpha that
    //              would produce the commanded normal acceleration.
    //              Falls back to attitude-hold when q is too low for
    //              the conversion to be reliable.
    int maut;

    // ---- Vacuum-mode autopilot parameters (maut=99) ----
    double vac_rate_damp;   // [s]   inner-loop rate time constant (default 0.5)
    double vac_max_gain;    // [-]   max inner-loop gain magnitude (rad/rad)

    // ---- Cascaded autopilot parameters (maut=50, maut=60) ----
    double tau_att;         // [s]   outer attitude-loop time constant
    double q_max;           // [deg/s] saturation on commanded body rate
    double theta_com_inertial;  // [deg] commanded body pitch (geodetic)
    double psi_com_inertial;    // [deg] commanded body yaw   (geodetic)

    // ---- Closed-loop pole parameters (gain-scheduled with q) ----
    // Pitch channel:
    double waclp;     // [rad/s]  natural frequency
    double zaclp;     // [-]      damping ratio
    double paclp;     // [-]      real pole
    double factwaclp; // [-]      modifier on waclp
    // Yaw channel:
    double wacly;
    double zacly;
    double pacly;
    double factwacly;

    // ---- Output limiters ----
    double delimx;    // [deg]   max pitch deflection magnitude
    double drlimx;    // [deg]   max yaw deflection magnitude

    // ---- Acceleration capability limits (gnmax/gymax in Zipfel) ----
    double gnmax;     // [g's]   max commandable normal accel
    double gymax;     // [g's]   max commandable lateral accel

    // ---- Acceleration commands (from Guidance; settable for testing) ----
    double ancomx;    // [g's]   commanded normal (pitch-plane) acceleration
    double alcomx;    // [g's]   commanded lateral (yaw-plane) acceleration

    // ---- States (integrator accumulates accel error) ----
    double zz, zzd;   // pitch integrator state and its derivative
    double yy, yyd;   // yaw integrator state and its derivative

    // ---- Outputs ----
    double delecx;    // [deg]   commanded pitch deflection
    double delrcx;    // [deg]   commanded yaw deflection
    double ancomx_actual;  // diagnostic: limited normal accel command
    double alcomx_actual;  // diagnostic: limited lateral accel command

    // ---- Diagnostic gains ----
    double gainfp1, gainfp2, gainfp3;
    double gainfy1, gainfy2, gainfy3;

    Control();
    void init()   override;
    void update() override;
    void rpt()    override;

    ACCESS_FN(double, delecx)
    ACCESS_FN(double, delrcx)
};

} // namespace rocket6dof

#endif
