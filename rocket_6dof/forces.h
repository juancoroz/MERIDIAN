//  forces.h  --  Forces and moments aggregator block
//
//  Sums force and moment contributions from upstream blocks into FAPB
//  (applied body force) and FMB (applied body moment).  Inputs sources
//  in current build:
//
//      Aerodynamics -> coefficients (cx,cy,cz / cll,clm,cln), scaled
//                      by pdynmc * refa (and refd for moments)
//      Propulsion   -> thrust along body +x (used when TVC is OFF)
//      TVC          -> rotated thrust force + gimbal moment (used when ON)
//      RCS          (not yet implemented; passes zero)
//
//  When TVC is active, the propulsion thrust is rotated through the
//  TVC deflection angles; Forces reads FPB and FMPB from TVC INSTEAD
//  of adding the unrotated thrust from Propulsion.
//
//  Also accepts an external override (FAPB_ext, FMB_ext) ADDED on top of
//  the aggregated value, so tests can inject arbitrary forces or moments
//  without needing the upstream blocks.

#ifndef ROCKET6DOF_FORCES_H
#define ROCKET6DOF_FORCES_H

#include "../osk/osk.h"

namespace rocket6dof {

class Propulsion;
class Aerodynamics;
class Environment;
class TVC;
class RCS;

class Forces : public osk::Block {
public:
    // ---- Inputs ----
    Propulsion*   prop;
    Aerodynamics* aero;
    Environment*  env;
    TVC*          tvc;
    RCS*          rcs;
    void getsFrom(Propulsion* p) {
        prop = p; aero = nullptr; env = nullptr; tvc = nullptr; rcs = nullptr;
    }
    void getsFrom(Propulsion* p, Aerodynamics* a, Environment* e) {
        prop = p; aero = a; env = e; tvc = nullptr; rcs = nullptr;
    }
    void getsFrom(Propulsion* p, Aerodynamics* a, Environment* e, TVC* t) {
        prop = p; aero = a; env = e; tvc = t; rcs = nullptr;
    }
    void getsFrom(Propulsion* p, Aerodynamics* a, Environment* e, TVC* t, RCS* r) {
        prop = p; aero = a; env = e; tvc = t; rcs = r;
    }

    // ---- External overrides (additive on top of aggregated value) ----
    osk::Vec FAPB_ext;
    osk::Vec FMB_ext;

    // ---- Outputs ----
    osk::Vec FAPB;    // [N]    summed applied force in body frame
    osk::Vec FMB;     // [N*m]  summed applied moment in body frame

    Forces()
        : prop(nullptr), aero(nullptr), env(nullptr), tvc(nullptr), rcs(nullptr),
          FAPB_ext(0, 0, 0), FMB_ext(0, 0, 0),
          FAPB(0, 0, 0), FMB(0, 0, 0)
    {}

    void init() override {}
    void update() override;
    void rpt()    override {}

    ACCESS_FN(osk::Vec, FAPB)
    ACCESS_FN(osk::Vec, FMB)
};

} // namespace rocket6dof

#endif
