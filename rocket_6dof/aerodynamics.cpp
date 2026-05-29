//  aerodynamics.cpp  --  Aero coefficients for the three-stage rocket
//
//  Cross-reference to Zipfel "Modeling and Simulation of Aerospace
//  Vehicle Dynamics" 3rd ed., Section 10.3.x (hyper[100-199]):
//
//      Zipfel name      OSK member       Meaning
//      -----------      ----------       -------
//      maero            maero            stage selector
//      refa             refa             reference area      [m^2]
//      refd             refd             reference length    [m]
//      xcg_ref          xcg_ref          reference CG        [m from nose]
//      ca0              ca0              Ca vs Mach (table)
//      caa              caa              dCa/dalpha vs Mach
//      ca0b             ca0b             base bleed correction
//      cn0              cn0              Cn vs Mach,alpha
//      clm0_v           clm0             Cm vs Mach,alpha
//      clmq_v           clmq             pitch damping vs Mach
//      cx, cy, cz       cx, cy, cz       body-axis force coefficients
//      cll, clm, cln    cll, clm, cln    body-axis moment coefficients
//
//  Equations (Zipfel aerodynamics.cpp lines 252-352):
//
//    axial force coeff:
//        ca = ca0(M) + caa(M)*alppx + (thrust_on ? ca0b(M) : 0)
//
//    normal force coeff (aeroballistic):
//        cna = cn0(M, alppx)
//
//    pitching moment coeff (aeroballistic), with pitch damping and CG shift:
//        clmaref = clm0(M, alppx) + clmq(M) * qqax * refd / (2*V)
//        clma    = clmaref - cna * (xcg_ref - xcg) / refd
//
//    body axes transformation using aero roll angle phipx:
//        cx  = -ca
//        cy  = -cna * sin(phipx)
//        cz  = -cna * cos(phipx)
//        cll = 0
//        clm =  clma * cos(phipx)
//        cln = -clma * sin(phipx)
//
//  The body-rate transformation from body to aeroballistic frame is:
//        qqax =  qqx * cos(phipx) - rrx * sin(phipx)
//        rrax =  qqx * sin(phipx) + rrx * cos(phipx)
//  Only qqax is used (pitch rate about aerobalistic-y axis).

#include "aerodynamics.h"
#include "environment.h"
#include "kinematics.h"
#include "propulsion.h"
#include "tvc.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double RAD = osk::PI / 180.0;
} // anon

Aerodynamics::Aerodynamics()
    : env(nullptr), kin(nullptr), prop(nullptr), tvc(nullptr),
      maero(0),
      refa(1.0), refd(1.0), xcg_ref(0.0),
      cx(0), cy(0), cz(0), cll(0), clm(0), cln(0),
      dla(0), dma(0), dmq(0), dmde(0),
      dyb(0), dnb(0), dnr(0), dndr(0),
      ca0(0), caa(0), ca0b(0), ca(0),
      cn0(0), cna(0),
      clm0_v(0), clmq_v(0), clma(0),
      tables_loaded_(false)
{
    // No integrator states; Aero is pure feed-forward.
}

void Aerodynamics::init() {
    if (initCount == 0 && maero != 0 && !aero_file.empty() && !tables_loaded_) {
        // Load aero tables on first init.  The AeroTable classes are
        // format-agnostic: they detect the file extension (.asc =
        // Missile DATCOM, anything else = OSK Table) and use the
        // appropriate parser internally.
        t_ca0_  = std::make_unique<AeroTable1>(aero_file);
        t_caa_  = std::make_unique<AeroTable1>(aero_file);
        t_ca0b_ = std::make_unique<AeroTable1>(aero_file);
        t_cn0_  = std::make_unique<AeroTable2>(aero_file);
        t_clm0_ = std::make_unique<AeroTable2>(aero_file);
        t_clmq_ = std::make_unique<AeroTable1>(aero_file);

        bool ok = true;
        if (!tag_ca0.empty())  ok &= t_ca0_->read(tag_ca0);
        if (!tag_caa.empty())  ok &= t_caa_->read(tag_caa);
        if (!tag_ca0b.empty()) ok &= t_ca0b_->read(tag_ca0b);
        if (!tag_cn0.empty())  ok &= t_cn0_->read(tag_cn0);
        if (!tag_clm0.empty()) ok &= t_clm0_->read(tag_clm0);
        if (!tag_clmq.empty()) ok &= t_clmq_->read(tag_clmq);

        if (!ok) {
            std::fprintf(stderr, "*** Aerodynamics::init: failed to load tables from %s ***\n",
                         aero_file.c_str());
        }
        tables_loaded_ = ok;
    }
}

void Aerodynamics::update() {
    // ---- Early-out: no aero if maero=0 or no env/kin/prop ----
    if (maero == 0 || !tables_loaded_ || !env || !kin) {
        cx = cy = cz = 0;
        cll = clm = cln = 0;
        dla = dma = dmq = dmde = 0;
        dyb = dnb = dnr = dndr = 0;
        return;
    }

    // ---- Read inputs ----
    double vmach  = env->vmach;
    double dvba   = env->dvba;
    double alppx  = kin->alppx;       // total angle of attack [deg]
    double phipx  = kin->phipx;       // aerodynamic roll angle [deg]

    // Body pitch rate (deg/s).  The damping term clmq*qqx*refd/(2V)
    // is small unless the vehicle is maneuvering hard, so qqx is held
    // at zero unless a propulsion stage exposes a body pitch rate via
    // Euler.
    double qqx    = 0.0;   // [deg/s]
    (void)qqx;             // suppress unused warning

    bool thrust_on = (prop && prop->thrust > 0.0);

    // ---- Axial force coefficient ----
    ca0 = (*t_ca0_) (vmach);
    caa = (*t_caa_) (vmach);
    ca0b = tag_ca0b.empty() ? 0.0 : (*t_ca0b_)(vmach);
    ca = ca0 + caa * alppx + (thrust_on ? ca0b : 0.0);

    // ---- Normal force coefficient (aeroballistic) ----
    cn0 = (*t_cn0_)(vmach, alppx);
    cna = cn0;

    // ---- Pitching moment coefficient (aeroballistic) ----
    clm0_v = (*t_clm0_)(vmach, alppx);
    clmq_v = tag_clmq.empty() ? 0.0 : (*t_clmq_)(vmach);
    double xcg = (prop ? prop->xcg : xcg_ref);
    // Pitch damping term clmq*qqx contributes via qqx, which is held
    // at zero (see qqx above), so this reduces to clm0_v alone.
    double clmaref = clm0_v;   // + clmq_v * qqx * RAD * refd / (2.0 * dvba)
    // The (xcg_ref - xcg) shift moves the moment reference point from the
    // aero-table's xcg_ref to the actual current CG.
    if (refd > 0.0)
        clma = clmaref - cna * (xcg_ref - xcg) / refd;
    else
        clma = clmaref;
    (void)dvba;   // unused; needed if pitch damping (qqx) is enabled

    // ---- Body-axis transformation via aero roll phipx ----
    double phip = phipx * RAD;
    double sphip = std::sin(phip);
    double cphip = std::cos(phip);

    cx  = -ca;
    cy  = -cna * sphip;
    cz  = -cna * cphip;

    cll =  0.0;            // no roll moment for axisymmetric vehicle
    clm =  clma * cphip;
    cln = -clma * sphip;

    // ---- Dimensional derivatives for Control ----
    // Zipfel aerodynamics.cpp lines 503-534.  We compute these from
    // finite-difference table lookups around the current operating
    // point, scaled by q*S/m and q*S*d/I_ii to give units of m/s^2/rad
    // for force derivatives and 1/s^2/rad for moment derivatives.
    //
    // For the AXISYMMETRIC rocket (no cross-coupling roll moment,
    // alpha and beta affect the body identically):
    //     cla -> sideforce slope cyb has equal magnitude opposite sign
    //     cma -> yaw moment cnb has equal magnitude opposite sign
    //     cmq -> yaw damping cnr has same value
    //
    // TVC actuator: control derivatives dmde, dndr come from rotated
    // thrust force lever arm, NOT from aero (no fins on this booster).
    {
        double q     = env->pdynmc;
        double m     = (prop ? prop->vmass : 1.0);
        double V     = (env->dvba > 1.0 ? env->dvba : 1.0);  // guard div-by-zero
        double xcg   = (prop ? prop->xcg : xcg_ref);
        double thr   = (prop ? prop->thrust : 0.0);

        // IBBB from propulsion (time-varying); use diagonal entries
        osk::Mat I_(1,0,0, 0,1,0, 0,0,1);
        if (prop) I_ = prop->IBBB;
        double Iy = I_[1][1];   // pitch MOI
        double Iz = I_[2][2];   // yaw MOI

        // Finite-difference cla, cma via +/- 3 deg around current alpha.
        // Clamp lower bound to 0 to stay inside table.
        double alpp = kin->alppx;
        double alpl = alpp + 3.0;
        double alpm = alpp - 3.0;
        if (alpm < 0.0) alpm = 0.0;

        double cn_p   = (*t_cn0_)(vmach, alpl);
        double cn_m   = (*t_cn0_)(vmach, alpm);
        double clm0_p = (*t_clm0_)(vmach, alpl);
        double clm0_m = (*t_clm0_)(vmach, alpm);

        double dalp = alpl - alpm;       // typ. 6 deg
        double cla  = (dalp > 0.0) ? (cn_p   - cn_m  ) / dalp : 0.0;   // per deg
        double cma_table = (dalp > 0.0) ? (clm0_p - clm0_m) / dalp : 0.0;
        // Apply CG shift correction (Zipfel line 341)
        double cma  = cma_table - cla * (xcg_ref - xcg) / refd;       // per deg

        // Convert per-deg coefficients to dimensional derivatives (per rad)
        // via the /RAD factor (Zipfel lines 504, 508-510, etc.)
        double duml = (q * refa / m)         / RAD;   // [m/s^2 per rad]
        double dumm = (q * refa * refd / Iy) / RAD;   // [1/s^2 per rad]
        double dumn = (q * refa * refd / Iz) / RAD;   // [1/s^2 per rad]

        dla = duml * cla;
        dma = dumm * cma;
        // Pitch damping: cmq is per rad in our table, no /RAD scaling
        double cmq = clmq_v;                          // table value, per rad
        dmq = (q * refa * refd / Iy) * (refd / (2.0 * V)) * cmq;

        // Yaw plane (axisymmetric: yaw mirrors pitch)
        dyb = -duml * cla;        // sign convention: cyb = -cla per Zipfel
        dnb = -dumn * cma;
        dnr =  (q * refa * refd / Iz) * (refd / (2.0 * V)) * cmq;

        // ---- Control derivatives from TVC actuator path ----
        // dmde = -(parm - xcg) * gtvc * thrust / Iy   (Zipfel line 530)
        // dndr = same form, for yaw                   (Zipfel line 533)
        // dlde = gtvc * thrust / m                    (Zipfel line 529)
        // No aero control surfaces on this booster -> dmde, dndr come
        // exclusively from the TVC actuator.  If no TVC wired, zero.
        if (tvc && tvc->mtvc != 0 && Iy > 0 && Iz > 0) {
            double gtvc = tvc->gtvc;
            double parm = tvc->parm;
            // Note the sign:  positive TVC deflection should create
            // negative pitching moment (nose-up), so this is negative.
            dmde = -(parm - xcg) * gtvc * thr / Iy;
            dndr = -(parm - xcg) * gtvc * thr / Iz;
        } else {
            dmde = 0.0;
            dndr = 0.0;
        }
    }
}

void Aerodynamics::rpt() {
    if (osk::State::sample(1.0)) {
        std::printf("Aero t=%7.3f  M=%5.3f  alpha=%6.2f  phip=%6.2f  "
                    "ca=%6.3f  cn=%6.3f  cx=%6.3f  cz=%6.3f\n",
                    osk::State::t,
                    env ? env->vmach : 0.0,
                    kin ? kin->alppx : 0.0,
                    kin ? kin->phipx : 0.0,
                    ca, cna, cx, cz);
    }
}

} // namespace rocket6dof
