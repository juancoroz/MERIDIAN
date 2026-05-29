//  environment.cpp  --  Atmosphere + gravity for the 6DOF rocket
//
//  Cross-reference to Zipfel "Modeling and Simulation of Aerospace
//  Vehicle Dynamics" 3rd ed., Section 10.3.2:
//
//      Zipfel name      OSK member       Meaning
//      -----------      ----------       -------
//      press            press            atmospheric pressure   [Pa]
//      rho              rho              density                [kg/m^3]
//      tempk, tempc     tempk, tempc     temperature            [K], [C]
//      vsound           vsound           speed of sound         [m/s]
//      vmach            vmach            Mach number            [-]
//      pdynmc           pdynmc           dynamic pressure       [Pa]
//      grav             grav             gravity magnitude      [m/s^2]
//      GRAVG            GRAVG            gravity vector, ECI    [m/s^2]
//      VAED             VAED             wind in geodetic frame [m/s]
//      dvba             dvba             vehicle speed wrt air  [m/s]
//
//  Departures from Zipfel's CADAC source:
//    - No internal integrator states.  His VAEDS (smoothed wind) and
//      taux1, taux2 (Dryden turbulence filter states) are deferred.
//      In v1, VAED is identically zero, mturb is implicitly 0.
//    - Gravity defaults to inverse-square spherical, not WGS84-J2.
//      Hook is left for upgrade (mgrav switch).
//    - Atmosphere model selection via 'matmo' int member, not the
//      packed 'mair' integer Zipfel decodes at every step.

#include "environment.h"
#include "newton.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

//  Physical constants for the 1976 US Standard Atmosphere
namespace {

constexpr double G0      = 9.80665;       // [m/s^2]   standard gravity
constexpr double R_STAR  = 8.31432;       // [J/(mol K)] universal gas const
constexpr double M_AIR   = 0.0289644;     // [kg/mol]    dry-air molar mass
constexpr double R_SPEC  = R_STAR / M_AIR;// [J/(kg K)]  ~287.053
constexpr double GAMMA   = 1.4;           // [-]         air ratio of specific heats
constexpr double R_EARTH = 6356766.0;     // [m]         effective Earth radius for geopotential
constexpr double MU_E    = 3.986004418e14;// [m^3/s^2]   Earth GM (WGS84)
constexpr double R_E_EQ  = 6378137.0;     // [m]         Earth equatorial radius (for default SBII)

// Seven-layer parameters of the 1976 Standard Atmosphere.
// h_b is GEOPOTENTIAL height in meters.
struct Layer {
    double h_b;    // base geopotential altitude [m]
    double T_b;    // base temperature           [K]
    double L_b;    // lapse rate                 [K/m]   (note: per meter!)
    double p_b;    // base pressure              [Pa]
};

constexpr Layer LAYERS[7] = {
    {     0.0, 288.15, -6.5e-3, 101325.0       },
    { 11000.0, 216.65,  0.0,     22632.0596    },
    { 20000.0, 216.65, +1.0e-3,   5474.8890    },
    { 32000.0, 228.65, +2.8e-3,    868.01870   },
    { 47000.0, 270.65,  0.0,       110.90631   },
    { 51000.0, 270.65, -2.8e-3,     66.938874  },
    { 71000.0, 214.65, -2.0e-3,      3.9564206 },
};

constexpr double H_TOP_GP = 84852.0;   // top of layer 6 [m geopotential] (~86 km geometric)

// Geometric -> geopotential altitude conversion (US76 Eq. 19).
inline double geometric_to_geopotential(double z) {
    return R_EARTH * z / (R_EARTH + z);
}

} // anonymous namespace

//  atmosphere76: US 1976 Standard Atmosphere
//  Input:  alt_m  geometric altitude [m]
//  Output: rho [kg/m^3], press [Pa], tempk [K]
//  Returns 0 on success, 1 if outside [0, 86 km] (values are clamped
//  to the nearest boundary in that case).
int atmosphere76(double alt_m, double& rho, double& press, double& tempk) {
    int oor = 0;

    // Clamp below sea level to sea level
    double z = alt_m;
    if (z < 0.0) { z = 0.0; oor = 1; }

    // Geopotential altitude
    double h = geometric_to_geopotential(z);

    // Clamp above 86 km to layer 6 top
    if (h > H_TOP_GP) { h = H_TOP_GP; oor = 1; }

    // Find which layer h falls in
    int b = 0;
    for (int i = 0; i < 7; ++i) {
        if (h >= LAYERS[i].h_b) b = i;
        else break;
    }

    const Layer& L = LAYERS[b];
    double dh = h - L.h_b;

    if (L.L_b == 0.0) {
        // Isothermal layer
        tempk = L.T_b;
        press = L.p_b * std::exp(-G0 * dh / (R_SPEC * L.T_b));
    } else {
        // Gradient layer
        tempk = L.T_b + L.L_b * dh;
        double base_ratio = L.T_b / tempk;
        double exponent   = G0 / (R_SPEC * L.L_b);
        press = L.p_b * std::pow(base_ratio, exponent);
    }

    rho = press / (R_SPEC * tempk);
    return oor;
}

//  grav_inv_sq: simple inverse-square gravity at inertial position
osk::Vec grav_inv_sq(const osk::Vec& SBII) {
    double r = SBII.mag();
    if (r < 1.0) return osk::Vec(0, 0, 0);   // guard against r=0
    double g_mag = MU_E / (r * r);
    // Acceleration points from body toward Earth center, i.e. -SBII / r
    osk::Vec g = SBII * (-g_mag / r);
    return g;
}

//  Environment block
Environment::Environment() {
    matmo  = 0;          // US76 by default
    mgrav  = 0;          // inverse-square by default
    newton = nullptr;    // optional; if wired, pulls inputs from Newton

    // Default inputs: stationary on the equator at sea level.
    alt    = 0.0;
    SBII   = osk::Vec(R_E_EQ, 0.0, 0.0);
    VBED   = osk::Vec(0.0, 0.0, 0.0);

    // Initialise outputs to sea-level values so any block that reads
    // them before the first update() gets a sane number rather than 0.
    press  = 101325.0;
    rho    = 1.225;
    tempk  = 288.15;
    tempc  = tempk - 273.15;
    vsound = std::sqrt(GAMMA * R_SPEC * tempk);
    vmach  = 0.0;
    pdynmc = 0.0;
    grav   = MU_E / (R_E_EQ * R_E_EQ);
    GRAVG  = osk::Vec(-grav, 0.0, 0.0);
    VAED   = osk::Vec(0.0, 0.0, 0.0);
    dvba   = 0.0;

    // No addIntegrator() calls -- Environment is feed-forward in v1.
}

void Environment::init() {
    // Nothing stage-specific yet.  Reserved for parameter-file load
    // (Filer) once we wire one up.
    (void)initCount;
}

void Environment::update() {
    // ---- Pull inputs from Newton if wired ----
    if (newton) {
        alt  = newton->alt;
        SBII = newton->SBII;
        VBED = newton->VBED;
    }

    // ---- Atmosphere ----
    switch (matmo) {
    case 0:
    default:
        atmosphere76(alt, rho, press, tempk);
        break;
    }
    tempc  = tempk - 273.15;
    vsound = std::sqrt(GAMMA * R_SPEC * tempk);

    // ---- Wind: zero in v1 ----
    VAED = osk::Vec(0.0, 0.0, 0.0);

    // ---- Airspeed ----
    // With zero wind, vehicle airspeed equals geodetic velocity magnitude.
    osk::Vec VBAD = VBED - VAED;
    dvba = VBAD.mag();

    // ---- Mach & dynamic pressure ----
    vmach  = (vsound > 0.0) ? std::fabs(dvba / vsound) : 0.0;
    pdynmc = 0.5 * rho * dvba * dvba;

    // ---- Gravity ----
    // NOTE: cad_grav84 returns GRAVG in GEOCENTRIC (G) frame, not inertial.
    // Newton must rotate it via ~TGI before using as inertial accel.  In
    // mgrav=0 (spherical) mode we return GRAVG already in inertial frame.
    switch (mgrav) {
    case 1:
        GRAVG = cad_grav84(SBII, osk::State::t);
        break;
    case 0:
    default:
        GRAVG = grav_inv_sq(SBII);
        break;
    }
    grav = GRAVG.mag();
}

void Environment::rpt() {
    // Default reporter: one line per second, easy to silence by
    // overriding or by raising the period.
    if (osk::State::sample(1.0)) {
        std::printf("Env  t=%7.3f  alt=%8.1f m  rho=%9.4e  p=%9.4e  T=%6.2f K  M=%5.3f  g=%6.4f\n",
                    osk::State::t, alt, rho, press, tempk, vmach, grav);
    }
}

} // namespace rocket6dof
