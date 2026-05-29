//  tvc.h  --  Thrust Vector Control block (instantaneous-actuator stub)
//
//  This block rotates the propulsion thrust vector by the deflection
//  angles commanded by Control, and produces the resulting body-frame
//  force (FPB) and moment (FMPB) about the CG.  The Forces block
//  aggregates FPB into FAPB and FMPB into FMB instead of using the
//  plain axial thrust from Propulsion.
//
//  V1 stub:
//    * Instantaneous actuator (no rate/position dynamics yet)
//    * delex = delecx_cmd, delrx = delrcx_cmd, clamped to mechanical
//      limits |del_max| once that's added
//
//  Cross-reference to Zipfel tvc.cpp (hyper[900-999]):
//      Zipfel name     OSK member       Meaning
//      -----------     ----------       -------
//      mtvc            mtvc             mode (0=off, 1=on)
//      gtvc            gtvc             gimbal effectiveness gain
//      parm            parm             gimbal location from nose [m]
//      eta, etax       eta, etax        pitch deflection [rad, deg]
//      zet, zetx       zet, zetx        yaw deflection [rad, deg]
//      FPB             FPB              rotated thrust force in body [N]
//      FMPB            FMPB             propulsion moment in body [N*m]

#ifndef ROCKET6DOF_TVC_H
#define ROCKET6DOF_TVC_H

#include "../osk/osk.h"

namespace rocket6dof {

class Propulsion;
class Control;

class TVC : public osk::Block {
public:
    // ---- Inputs ----
    Propulsion* prop;
    Control*    control;          // optional; when wired, auto-pulls
                                  // delecx_cmd and delrcx_cmd from
                                  // Control->delecx/delrcx each update
    void getsFrom(Propulsion* p) { prop = p; control = nullptr; }
    void getsFrom(Propulsion* p, Control* c) { prop = p; control = c; }

    // ---- Mode ----
    int mtvc;          // 0 = TVC off, 1 = on (stub instantaneous)

    // ---- Parameters ----
    double gtvc;       // [-]  gimbal effectiveness gain
    double parm;       // [m]  gimbal location from nose
    double del_max;    // [deg] mechanical limit on |deflection|

    // ---- Commands (set by Control) ----
    double delecx_cmd; // [deg]  commanded pitch deflection
    double delrcx_cmd; // [deg]  commanded yaw deflection

    // ---- Outputs ----
    double etax;       // [deg]  achieved pitch deflection
    double zetx;       // [deg]  achieved yaw deflection
    osk::Vec FPB;      // [N]    rotated thrust force, body frame
    osk::Vec FMPB;     // [N*m]  thrust-induced moment about CG, body frame

    TVC()
        : prop(nullptr), control(nullptr),
          mtvc(0),
          gtvc(1.0), parm(0.0), del_max(15.0),
          delecx_cmd(0.0), delrcx_cmd(0.0),
          etax(0.0), zetx(0.0),
          FPB(0, 0, 0), FMPB(0, 0, 0)
    {}

    void init() override {}
    void update() override;
    void rpt() override {}

    // Getters used by Aerodynamics (gtvc, parm) and Forces (FPB, FMPB)
    ACCESS_FN(double,   gtvc)
    ACCESS_FN(double,   parm)
    ACCESS_FN(double,   etax)
    ACCESS_FN(double,   zetx)
    ACCESS_FN(osk::Vec, FPB)
    ACCESS_FN(osk::Vec, FMPB)
};

} // namespace rocket6dof

#endif
