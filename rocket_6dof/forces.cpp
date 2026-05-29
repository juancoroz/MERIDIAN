//  forces.cpp  --  Aggregate body-frame forces and moments

#include "forces.h"
#include "propulsion.h"
#include "aerodynamics.h"
#include "environment.h"
#include "tvc.h"
#include "rcs.h"

namespace rocket6dof {

void Forces::update() {
    // Start from external override (zero in normal use)
    FAPB = FAPB_ext;
    FMB  = FMB_ext;

    // ---- Aerodynamic contribution ----
    if (aero && env) {
        double q = env->pdynmc;
        double S = aero->refa;
        double L = aero->refd;
        double qS = q * S;

        FAPB.x += qS * aero->cx;
        FAPB.y += qS * aero->cy;
        FAPB.z += qS * aero->cz;

        double qSL = qS * L;
        FMB.x += qSL * aero->cll;
        FMB.y += qSL * aero->clm;
        FMB.z += qSL * aero->cln;
    }

    // ---- Propulsion contribution ----
    bool tvc_active = (tvc && tvc->mtvc != 0);
    if (tvc_active) {
        FAPB.x += tvc->FPB.x;
        FAPB.y += tvc->FPB.y;
        FAPB.z += tvc->FPB.z;
        FMB.x  += tvc->FMPB.x;
        FMB.y  += tvc->FMPB.y;
        FMB.z  += tvc->FMPB.z;
    } else if (prop) {
        FAPB.x += prop->thrust;
    }

    // ---- RCS contribution ----
    // Body-frame side forces and torques produced by reaction thrusters.
    // Zero when RCS modes are off.
    if (rcs) {
        FAPB.x += rcs->FARCS.x;
        FAPB.y += rcs->FARCS.y;
        FAPB.z += rcs->FARCS.z;
        FMB.x  += rcs->FMRCS.x;
        FMB.y  += rcs->FMRCS.y;
        FMB.z  += rcs->FMRCS.z;
    }
}

} // namespace rocket6dof
