//  gps.h  --  GPS receiver block
//
//  Three operating modes:
//
//    mgps = 0  off
//    mgps = 1  perfect GPS (measurement = truth at gps_step interval)
//    mgps = 2  simple noisy GPS (truth + Gaussian noise at gps_step)
//    mgps = 3  full Zipfel Kalman filter GPS
//
//  Modes 1 and 2 are the simplified verification-friendly modes.
//  Mode 3 implements Zipfel's complete pseudorange-based Kalman filter:
//
//    * 24-satellite GPS constellation (6 orbital planes, 4 SVs each,
//      Yuma Almanac Week 787 data, 26560 km circular orbits, 55 deg
//      inclination)
//    * Satellite ephemeris propagation in time
//    * Visibility test (line-of-sight clearance over Earth)
//    * Quadriga selection (best 4 SVs by minimum GDOP)
//    * Pseudo-range and delta-range measurements with bias + noise
//      + user-clock-bias error
//    * 8-state Kalman filter (position 3, velocity 3, clock bias,
//      clock frequency) with covariance propagation and update
//    * Output: SXH (position correction), VXH (velocity correction)
//      for the INS to apply
//
//  Cross-reference to Zipfel gps.cpp (hyper[700-799]).

#ifndef ROCKET6DOF_GPS_H
#define ROCKET6DOF_GPS_H

#include "../osk/osk.h"
#include "kalman_matrix.h"

namespace rocket6dof {

class Newton;
class Euler;
class INS;

class GPS : public osk::Block {
public:
    // ---- Inputs ----
    Newton* newton;
    Euler*  euler;     // truth angular rate for satellite range-rate
    INS*    ins;       // INS-predicted state for filter residuals (mgps=3)
    void getsFrom(Newton* n) {
        newton = n; euler = nullptr; ins = nullptr;
    }
    void getsFrom(Newton* n, Euler* e, INS* i) {
        newton = n; euler = e; ins = i;
    }

    // ---- Mode selector ----
    int mgps;

    // ---- Common parameters ----
    double gps_step;       // [s]    interval between GPS updates
    double gps_epoch;      // [s]    time of last produced measurement
    double t_first;        // [s]    earliest time measurements can begin

    // ---- Modes 1/2 parameters ----
    double rpos;           // [m]    1-sigma position noise (mode 2)
    double rvel;           // [m/s]  1-sigma velocity noise (mode 2)
    unsigned long noise_seed;

    // ---- Mode 3 (Zipfel Kalman) parameters ----
    double almanac_time;   // [s]   time since GPS almanac epoch at sim start
    double del_rearth;     // [m]   added to Earth radius for LOS calculation
    double gps_acqtime;    // [s]   initial acquisition delay (first update)
    // User-clock error model
    double ucfreq_noise;   // [m/s] user-clock frequency Markov 1-sigma
    double ucbias_error;   // [m]   accumulated user-clock bias state
    double ucfreqm;        // [m/s] user-clock frequency state
    double uctime_cor;     // [s]   clock time-correlation constant
    // Per-satellite pseudorange and delta-range bias/noise (4 active SVs)
    double pr_bias[4];     // [m]   pseudorange bias 1-sigma
    double pr_noise[4];    // [m]   pseudorange noise 1-sigma
    double dr_noise[4];    // [m/s] delta-range noise 1-sigma
    // Filter initial covariance (P) diagonal 1-sigma values
    double ppos, pvel, pclockb, pclockf;
    // Process noise covariance (Q) diagonal 1-sigma values
    double qpos, qvel, qclockb, qclockf;
    // Measurement noise covariance (R) 1-sigma values (for KF mode)
    double rpos_kf, rvel_kf;
    // Tuning multipliers (initial: factp, process: factq, meas: factr)
    double factp, factq, factr;

    // ---- Kalman filter state ----
    MatN PP;          // 8x8 covariance matrix
    MatN FF;          // 8x8 dynamics matrix
    MatN PHI;         // 8x8 state-transition matrix (per gps_step)
    int  gps_acq;     // 1 = first acquisition delay still pending
    int  kf_phase;    // internal: 1 = init pending, 2 = extrapolate, 3 = update this cycle
    int  meas_count;
    double slotsum;   // diagnostic: sum of quadriga slot IDs (for change detection)

    // ---- Outputs ----
    osk::Vec SBII_meas;       // mode 1/2: GPS measurement of inertial position
    osk::Vec VBII_meas;       // mode 1/2: GPS measurement of inertial velocity
    osk::Vec SXH;             // mode 3: Kalman position correction state
    osk::Vec VXH;             // mode 3: Kalman velocity correction state
    osk::Vec CXH;             // mode 3: clock bias/freq state (2 elements, .x=bias, .y=freq)
    int      gps_update_avail;
    double   last_pos_err;
    double   last_vel_err;
    double   gdop;            // mode 3: GDOP of selected quadriga
    int      slot[4];         // mode 3: slot numbers (1..24) of active SVs

    GPS();
    void init()   override;
    void update() override;
    void rpt()    override;

    ACCESS_FN(osk::Vec, SBII_meas)
    ACCESS_FN(osk::Vec, VBII_meas)
    ACCESS_FN(osk::Vec, SXH)
    ACCESS_FN(osk::Vec, VXH)
    ACCESS_FN(int,      gps_update_avail)

private:
    // Internal helpers for mgps=3
    void kf_init();
    void kf_extrapolate();
    void kf_update();
    // Constellation: propagate to current time, find 4 best SVs.
    // Returns false if fewer than 4 satellites are visible.
    bool quadriga_select(double t,
                         osk::Vec ssii_quad[4], osk::Vec vsii_quad[4],
                         int slot_out[4], double& gdop_out);
};

} // namespace rocket6dof

#endif
