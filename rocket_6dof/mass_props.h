//  mass_props.h  --  Stub for vehicle mass properties, v2
//
//  Real implementation would track mass and inertia as they change
//  during fuel burn.  For testing we expose:
//    vmass  -- scalar mass [kg]
//    IBBB   -- inertia tensor in body frame [kg*m^2], a 3x3 matrix
//
//  Defaults: vmass = 1 kg, IBBB = identity.  Set externally for tests.
//  Zipfel stores this in hyper[15] (vmass) and hyper[18] (IBBB).

#ifndef ROCKET6DOF_MASS_PROPS_H
#define ROCKET6DOF_MASS_PROPS_H

#include "../osk/osk.h"

namespace rocket6dof {

class MassProps : public osk::Block {
public:
    double   vmass;   // [kg]
    osk::Mat IBBB;    // [kg*m^2] inertia tensor about CG, in body frame

    MassProps()
        : vmass(1.0),
          IBBB(1,0,0, 0,1,0, 0,0,1)    // default identity
    {}
    void init()   override {}
    void update() override {}
    void rpt()    override {}

    ACCESS_FN(double,   vmass)
    ACCESS_FN(osk::Mat, IBBB)
};

} // namespace rocket6dof

#endif
