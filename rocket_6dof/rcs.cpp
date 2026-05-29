//  rcs.cpp -- Reaction Control System (vacuum-phase attitude actuator)
//
//  Ported from Zipfel rcs.cpp.  Implements proportional and Schmitt-
//  trigger thrusters in body axes for both moment (attitude) and force
//  (side-thrust) control.

#include "rcs.h"
#include "propulsion.h"
#include "newton.h"
#include "kinematics.h"
#include "ins.h"
#include "guidance.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double RAD = osk::PI / 180.0;
constexpr double DEG = 180.0 / osk::PI;
} // anon

RCS::RCS()
    : prop(nullptr), newton(nullptr), kin(nullptr),
      ins(nullptr), guidance(nullptr),
      mrcs_moment(0), mrcs_force(0),
      rcs_freq(1.0), rcs_zeta(0.707),
      dead_zone(0.1), hysteresis(0.05), rcs_tau(0.5),
      roll_mom_max(100.0), pitch_mom_max(100.0), yaw_mom_max(100.0),
      side_force_max(100.0), acc_gain(1.0),
      phibdcomx(0), thtbdcomx(0), psibdcomx(0),
      alphacomx(0), betacomx(0),
      aycomx(0), azcomx(0),
      roll_save(0), pitch_save(0), yaw_save(0),
      right_save(0), down_save(0),
      o_roll(0), o_pitch(0), o_yaw(0),
      o_right(0), o_down(0),
      roll_count(0), pitch_count(0), yaw_count(0),
      right_count(0), down_count(0),
      FMRCS(0,0,0), FARCS(0,0,0),
      e_roll(0), e_pitch(0), e_yaw(0),
      e_right(0), e_down(0)
{}

void RCS::init() {
    if (initCount == 0) {
        FMRCS = osk::Vec(0, 0, 0);
        FARCS = osk::Vec(0, 0, 0);
        roll_save = pitch_save = yaw_save = 0.0;
        right_save = down_save = 0.0;
        o_roll = o_pitch = o_yaw = 0;
        o_right = o_down = 0;
        roll_count = pitch_count = yaw_count = 0;
        right_count = down_count = 0;
        e_roll = e_pitch = e_yaw = 0.0;
        e_right = e_down = 0.0;
    }
}

// ---- Proportional thruster output with saturation ----
double RCS::rcs_prop(double input, double limit) {
    if (std::fabs(input) > limit) {
        return (input >= 0.0) ? limit : -limit;
    }
    return input;
}

// ---- Schmitt trigger: dead-zone + hysteresis ----
// Returns +1, 0, or -1 based on a state machine combining current input
// magnitude and the change-direction (trend).
int RCS::rcs_schmitt(double input_new, double input_prev,
                     double dz, double hyst) {
    int trend = 0;
    if      (input_new > input_prev) trend = 1;
    else if (input_new < input_prev) trend = -1;
    int side = 0;
    if      (input_prev > 0) side = 1;
    else if (input_prev < 0) side = -1;

    double trigger = (dz * side + hyst * trend) * 0.5;

    if (input_prev >= trigger && side == 1)  return 1;
    if (input_prev <= trigger && side == -1) return -1;
    return 0;
}

// UPDATE_MARKER

void RCS::update() {
    FMRCS = osk::Vec(0, 0, 0);
    FARCS = osk::Vec(0, 0, 0);

    // Decode mrcs_moment = 10*rcs_type + rcs_mode
    int rcs_type = mrcs_moment / 10;
    int rcs_mode = mrcs_moment % 10;

    // Pull state from INS (or Kinematics/Newton as fallback)
    double ppcx_loc = 0, qqcx_loc = 0, rrcx_loc = 0;
    osk::Vec FSPCB_loc(0, 0, 0);
    if (ins) {
        ppcx_loc = ins->ppcx;
        qqcx_loc = ins->qqcx;
        rrcx_loc = ins->rrcx;
        FSPCB_loc = ins->FSPCB;
    }
    // Euler angles and incidence from truth (Kinematics/Newton).  In a
    // production system INS would compute these from TBIC; for v1 we
    // use truth.
    double phibdcx_loc = 0, thtbdcx_loc = 0, psibdcx_loc = 0;
    double alphacx_loc = 0, betacx_loc = 0;
    if (kin) {
        phibdcx_loc = kin->phibdx;
        thtbdcx_loc = kin->thtbdx;
        psibdcx_loc = kin->psibdx;
        alphacx_loc = kin->alphax;   // already in deg
        betacx_loc  = kin->betax;    // already in deg
    }
    osk::Vec UTBC_loc(0, 0, 0);
    if (guidance) UTBC_loc = guidance->UTBC;

    // Moment of inertia from Propulsion (diagonal IBBB)
    double Ixx = 1.0, Iyy = 1.0, Izz = 1.0;
    if (prop) {
        osk::Mat I = prop->IBBB;
        Ixx = I[0][0];  Iyy = I[1][1];  Izz = I[2][2];
    }

    // ---- Proportional moment thrusters (rcs_type=1) ----
    if (rcs_type == 1) {
        double rgain_roll  = 2.0 * rcs_zeta * rcs_freq * Ixx;
        double rgain_pitch = 2.0 * rcs_zeta * rcs_freq * Iyy;
        double rgain_yaw   = 2.0 * rcs_zeta * rcs_freq * Izz;
        double pgain       = rcs_freq / (2.0 * rcs_zeta);

        // Roll is always controlled
        e_roll = rgain_roll * (pgain * (phibdcomx - phibdcx_loc) - ppcx_loc);
        FMRCS.x = rcs_prop(e_roll, roll_mom_max);

        if (rcs_mode == 1) {
            // Geodetic Euler-angle control
            e_pitch = rgain_pitch * (pgain * (thtbdcomx - thtbdcx_loc) - qqcx_loc);
            e_yaw   = rgain_yaw   * (pgain * (psibdcomx - psibdcx_loc) - rrcx_loc);
        }
        else if (rcs_mode == 2) {
            // Thrust-vector direction control (UTBC small-angle in body frame)
            // DEG appears because Zipfel keeps body rates in deg/s.
            e_pitch = rgain_pitch * (pgain * (-UTBC_loc.z) * DEG - qqcx_loc);
            e_yaw   = rgain_yaw   * (pgain *  (UTBC_loc.y) * DEG - rrcx_loc);
        }
        else if (rcs_mode == 3) {
            // Incidence (alpha, beta) control
            e_pitch = rgain_pitch * (pgain * (alphacomx - alphacx_loc) - qqcx_loc);
            e_yaw   = rgain_yaw   * (pgain * (-betacomx + betacx_loc)  - rrcx_loc);
        }
        FMRCS.y = rcs_prop(e_pitch, pitch_mom_max);
        FMRCS.z = rcs_prop(e_yaw,   yaw_mom_max);
    }

    // ---- On-off (Schmitt) moment thrusters (rcs_type=2) ----
    if (rcs_type == 2) {
        // Roll is always controlled
        e_roll = phibdcomx - (rcs_tau * ppcx_loc + phibdcx_loc);
        int o_roll_save = o_roll;
        o_roll = rcs_schmitt(e_roll, roll_save, dead_zone, hysteresis);
        roll_save = e_roll;
        if (o_roll != o_roll_save) roll_count++;

        if (rcs_mode == 1) {
            e_pitch = thtbdcomx - (rcs_tau * qqcx_loc + thtbdcx_loc);
            e_yaw   = psibdcomx - (rcs_tau * rrcx_loc + psibdcx_loc);
        }
        else if (rcs_mode == 2) {
            e_pitch = -rcs_tau * qqcx_loc - UTBC_loc.z * DEG;
            e_yaw   = -rcs_tau * rrcx_loc + UTBC_loc.y * DEG;
        }
        else if (rcs_mode == 3) {
            e_pitch =  alphacomx - (rcs_tau * qqcx_loc + alphacx_loc);
            e_yaw   = -betacomx  - (rcs_tau * rrcx_loc - betacx_loc);
        }
        int o_pitch_save = o_pitch;
        o_pitch = rcs_schmitt(e_pitch, pitch_save, dead_zone, hysteresis);
        pitch_save = e_pitch;
        if (o_pitch != o_pitch_save) pitch_count++;

        int o_yaw_save = o_yaw;
        o_yaw = rcs_schmitt(e_yaw, yaw_save, dead_zone, hysteresis);
        yaw_save = e_yaw;
        if (o_yaw != o_yaw_save) yaw_count++;

        FMRCS.x = o_roll  * roll_mom_max;
        FMRCS.y = o_pitch * pitch_mom_max;
        FMRCS.z = o_yaw   * yaw_mom_max;
    }

    // ---- Proportional side-force thrusters (mrcs_force=1) ----
    if (mrcs_force == 1) {
        const double G0 = 9.80675445;
        double ay = FSPCB_loc.y;
        double az = FSPCB_loc.z;
        e_right = acc_gain * (aycomx * G0 - ay);
        e_down  = acc_gain * (azcomx * G0 - az);
        FARCS.x = 0.0;
        FARCS.y = rcs_prop(e_right, side_force_max);
        FARCS.z = rcs_prop(e_down,  side_force_max);
    }

    // ---- On-off side-force thrusters (mrcs_force=2) ----
    if (mrcs_force == 2) {
        const double G0 = 9.80675445;
        double ay = FSPCB_loc.y;
        double az = FSPCB_loc.z;
        e_right = acc_gain * (aycomx * G0 - ay);
        e_down  = acc_gain * (azcomx * G0 - az);
        int o_right_save = o_right;
        o_right = rcs_schmitt(e_right, right_save, dead_zone, hysteresis);
        right_save = e_right;
        if (o_right != o_right_save) right_count++;
        int o_down_save = o_down;
        o_down = rcs_schmitt(e_down, down_save, dead_zone, hysteresis);
        down_save = e_down;
        if (o_down != o_down_save) down_count++;
        FARCS.x = 0.0;
        FARCS.y = o_right * side_force_max;
        FARCS.z = o_down  * side_force_max;
    }
}

void RCS::rpt() {
    if (osk::State::sample(1.0)) {
        if (mrcs_moment != 0 || mrcs_force != 0) {
            std::printf("RCS  t=%7.3f  mom=%d force=%d  "
                        "FMRCS=(%+.2f,%+.2f,%+.2f) Nm  "
                        "FARCS=(%+.2f,%+.2f,%+.2f) N\n",
                        osk::State::t, mrcs_moment, mrcs_force,
                        FMRCS.x, FMRCS.y, FMRCS.z,
                        FARCS.x, FARCS.y, FARCS.z);
        }
    }
}

} // namespace rocket6dof
