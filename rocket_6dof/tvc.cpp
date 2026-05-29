//  tvc.cpp  --  Rotate thrust vector, compute body-frame F/M

#include "tvc.h"
#include "propulsion.h"
#include "control.h"
#include <cmath>

namespace rocket6dof {

namespace {
constexpr double DEG = 180.0 / osk::PI;
constexpr double RAD = osk::PI / 180.0;
inline double sign(double x) { return (x > 0.0) - (x < 0.0); }
} // anon

void TVC::update() {
    if (mtvc == 0 || !prop) {
        etax = 0.0;
        zetx = 0.0;
        FPB  = osk::Vec(0, 0, 0);
        FMPB = osk::Vec(0, 0, 0);
        return;
    }

    // Auto-pull commands from Control when wired.  Falls back to the
    // direct delecx_cmd/delrcx_cmd fields when control is not connected
    // (allows external override for unit tests).
    if (control) {
        delecx_cmd = control->delecx;
        delrcx_cmd = control->delrcx;
    }

    // Apply mechanical limits to commanded deflection (Zipfel uses
    // separate delimx/drlimx in Control; we honor those here too as a
    // backup actuator limit)
    double eta_cmd_deg = delecx_cmd;
    double zet_cmd_deg = delrcx_cmd;
    if (std::fabs(eta_cmd_deg) > del_max) eta_cmd_deg = del_max * sign(eta_cmd_deg);
    if (std::fabs(zet_cmd_deg) > del_max) zet_cmd_deg = del_max * sign(zet_cmd_deg);

    // Stub: instantaneous actuator response (real TVC will lag this)
    etax = eta_cmd_deg;
    zetx = zet_cmd_deg;

    double eta = etax * RAD;
    double zet = zetx * RAD;

    // Rotated thrust force in body axes (Zipfel tvc.cpp lines 118-120)
    //   FPB = (T*cos(eta)*cos(zet), T*cos(eta)*sin(zet), -T*sin(eta))
    // Positive eta produces -z body force (nose-up turn).
    double thr = prop->thrust * gtvc;
    double seta = std::sin(eta);
    double ceta = std::cos(eta);
    double szet = std::sin(zet);
    double czet = std::cos(zet);

    FPB.x =  ceta * czet * thr;
    FPB.y =  ceta * szet * thr;
    FPB.z = -seta * thr;

    // Moment about CG from the gimbal point arm (Zipfel lines 123-126):
    //   arm = parm - xcg   (negative when gimbal is aft of CG)
    //   FMPB = (0, arm*FPB.z, -arm*FPB.y)
    double xcg = prop->xcg;
    double arm = parm - xcg;
    FMPB.x = 0.0;
    FMPB.y =  arm * FPB.z;
    FMPB.z = -arm * FPB.y;
}

} // namespace rocket6dof
