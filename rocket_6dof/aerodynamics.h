//  aerodynamics.h  --  Aerodynamic force/moment coefficients for the
//                      three-stage rocket
//
//  Public members match the "out"-tagged module variables of Zipfel,
//  "Modeling and Simulation of Aerospace Vehicle Dynamics" 3rd ed.,
//  Section 10.3.x (hyper[100-199]).
//
//  This block produces NON-DIMENSIONAL coefficients (cx, cy, cz for
//  forces; cll, clm, cln for moments).  The Forces block multiplies
//  them by pdynmc * refa (and additionally by refd for moments) to
//  get N and N*m, then aggregates with propulsion and RCS contributions.
//
//  Selector 'maero':
//    0  = no aerodynamic forces (vacuum)
//    11 = configuration with 1 upper stage (last-stage solo)
//    12 = configuration with 2 upper stages (mid-stage)
//    13 = configuration with 3 upper stages (launch / first-stage)
//
//  Each stage uses its own set of aero tables in the data deck.
//  Coefficients are read via Table1 (Mach only) and Table2 (Mach +
//  total angle of attack).
//
//  Departures from Zipfel:
//    * No autopilot dimensional derivatives (dla, dma, dmq, ...) -- these
//      live in Control, not Aero.
//    * No analytical airframe stability roots (wnp, zetp, ...).
//    * No control-surface effectiveness terms (cmde, cnde, ...);
//      vehicles in this codebase steer via TVC, not control fins.
//    * No 'mfreeze' autopilot freeze logic.

#ifndef ROCKET6DOF_AERODYNAMICS_H
#define ROCKET6DOF_AERODYNAMICS_H

#include "../osk/osk.h"
#include "aero_deck.h"
#include <string>
#include <memory>

namespace rocket6dof {

class Environment;
class Kinematics;
class Propulsion;
class TVC;

class Aerodynamics : public osk::Block {
public:
    // ---- Inputs ----
    Environment* env;
    Kinematics*  kin;
    Propulsion*  prop;
    TVC*         tvc;        // optional; supplies gtvc, parm for TVC-mode derivatives
    void getsFrom(Environment* e, Kinematics* k, Propulsion* p) {
        env = e; kin = k; prop = p; tvc = nullptr;
    }
    void getsFrom(Environment* e, Kinematics* k, Propulsion* p, TVC* t) {
        env = e; kin = k; prop = p; tvc = t;
    }

    // ---- Stage / mode selector ----
    int    maero;   // 0 = off, 11/12/13 = stage configurations

    // ---- Geometric parameters (Zipfel "data" tag) ----
    double refa;      // [m^2]  reference area (e.g., body cross-section)
    double refd;      // [m]    reference length (e.g., body diameter)
    double xcg_ref;   // [m]    reference CG location from nose (where moments table is centered)

    // ---- Aero data file path and per-stage table tags ----
    // Filer parses the file; tags are looked up by name.  We allow per-
    // stage tag prefixes (e.g., "ca0slv3_vs_mach") rather than hardcoding,
    // so the same data file can serve multiple stages.
    std::string aero_file;       // path to data file (e.g., "aero.txt")
    std::string tag_ca0;         // axial force coefficient vs Mach (1-D)
    std::string tag_caa;         // dCa/dalpha vs Mach (1-D)
    std::string tag_ca0b;        // base bleed correction when thrust is on (1-D)
    std::string tag_cn0;         // normal force coefficient vs Mach, alpha (2-D)
    std::string tag_clm0;        // pitching moment coefficient vs Mach, alpha (2-D)
    std::string tag_clmq;        // pitch damping derivative vs Mach (1-D)

    // ---- Aero coefficient outputs (body axes; Zipfel "out") ----
    double cx, cy, cz;          // force coefficients (non-dimensional)
    double cll, clm, cln;       // moment coefficients (non-dimensional)

    // ---- Dimensional derivatives (Zipfel "out", used by Control) ----
    // Pitch plane:
    double dla;       // [1/s]    lift slope     = q*S/m * cla     [m/s^2 per rad alpha]
    double dma;       // [1/s^2]  pitch moment derivative wrt alpha
    double dmq;       // [1/s]    pitch damping derivative
    double dmde;      // [1/s^2]  pitch control derivative (TVC or fin)
    // Yaw plane:
    double dyb;       // [1/s]    sideforce slope wrt beta
    double dnb;       // [1/s^2]  yaw moment derivative wrt beta
    double dnr;       // [1/s]    yaw damping derivative
    double dndr;      // [1/s^2]  yaw control derivative (TVC or fin)

    // ---- Diagnostics (Zipfel "diag") ----
    double ca0, caa, ca0b;      // raw axial coefficients from tables
    double ca;                  // total axial coefficient
    double cn0, cna;            // raw and total normal coefficients
    double clm0_v, clmq_v;      // pitching moment table values
    double clma;                // total pitching moment coefficient (aerobalistic)

    Aerodynamics();
    void init()   override;
    void update() override;
    void rpt()    override;

    // ---- Getters for Forces ----
    ACCESS_FN(double, cx)
    ACCESS_FN(double, cy)
    ACCESS_FN(double, cz)
    ACCESS_FN(double, cll)
    ACCESS_FN(double, clm)
    ACCESS_FN(double, cln)

    // ---- Getters for Control (dimensional derivatives) ----
    ACCESS_FN(double, dla)
    ACCESS_FN(double, dma)
    ACCESS_FN(double, dmq)
    ACCESS_FN(double, dmde)
    ACCESS_FN(double, dyb)
    ACCESS_FN(double, dnb)
    ACCESS_FN(double, dnr)
    ACCESS_FN(double, dndr)

private:
    // Owned aero tables, loaded from disk in init().  The AeroTable
    // classes wrap both the OSK and Missile DATCOM file formats; the
    // file extension determines which parser is used.
    std::unique_ptr<AeroTable1> t_ca0_;
    std::unique_ptr<AeroTable1> t_caa_;
    std::unique_ptr<AeroTable1> t_ca0b_;
    std::unique_ptr<AeroTable2> t_cn0_;
    std::unique_ptr<AeroTable2> t_clm0_;
    std::unique_ptr<AeroTable1> t_clmq_;
    bool tables_loaded_;
};

} // namespace rocket6dof

#endif
