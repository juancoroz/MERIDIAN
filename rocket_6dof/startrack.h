//  startrack.h  --  Star-tracker attitude reference sensor
//
//  Produces discrete-rate attitude corrections that INS applies to its
//  TBIC estimate.  The output URIC is a small-angle tilt-error vector:
//
//      URIC.x = roll-axis tilt error  (small angle, rad)
//      URIC.y = pitch-axis tilt error
//      URIC.z = yaw-axis tilt error
//
//  INS applies the correction by:
//      TBIC := (I - [URIC]_x) * TBIC
//  followed by Bar-Itzhack orthonormalization.
//
//  Operating modes:
//    mstar = 0  off, no measurements
//    mstar = 1  perfect attitude (URIC = exact INS->truth correction)
//    mstar = 2  noisy attitude (small Gaussian tilt error per axis)
//    mstar = 3  full Zipfel triad-based star tracker (25-star catalog,
//               selects 3 stars with maximum parallelepiped volume,
//               measures each in body frame with per-axis az/el noise,
//               computes RDIFF = TRIAD_MEAS * TRIAD_TRUE^-1, extracts
//               URIC from off-diagonal elements).
//
//  Operational constraint: in mode 3, measurements only occur when
//  vehicle altitude is above startrack_alt (above the atmosphere).

#ifndef ROCKET6DOF_STARTRACK_H
#define ROCKET6DOF_STARTRACK_H

#include "../osk/osk.h"

namespace rocket6dof {

class Newton;
class Kinematics;
class INS;

class Startrack : public osk::Block {
public:
    Newton*     newton;
    Kinematics* kin;
    INS*        ins;          // read TBIC for residual computation in mode 3
    void getsFrom(Newton* n, Kinematics* k) {
        newton = n; kin = k; ins = nullptr;
    }
    void getsFrom(Newton* n, Kinematics* k, INS* i) {
        newton = n; kin = k; ins = i;
    }

    int    mstar;
    double startrack_step;   // [s] update interval (1 Hz typical)
    double starfix_epoch;    // [s] time of last produced measurement
    double t_first;          // [s] earliest time measurements can begin
    double startrack_alt;    // [m] minimum altitude for star tracking (mode 3)
    double star_acqtime;     // [s] initial acquisition delay (mode 3)
    int    star_acq;         // 1 = first acquisition pending

    // Per-axis noise 1-sigma (modes 2 and 3)
    double tilt_noise;       // [rad] per-axis tilt error 1-sigma (mode 2)
    // Per-star azimuth/elevation noise 1-sigma (mode 3)
    double az_bias[3];       // [rad] per-star azimuth bias
    double az_noise[3];      // [rad] per-star azimuth noise 1-sigma
    double el_bias[3];       // [rad] per-star elevation bias
    double el_noise[3];      // [rad] per-star elevation noise 1-sigma

    double star_el_min;      // [deg] minimum star elevation from local horizon

    unsigned long noise_seed;

    // Outputs
    osk::Vec URIC;           // tilt error vector (small angle) for INS
    int      star_update_avail;
    int      meas_count;
    int      triad[3];       // mode 3: slot numbers of selected stars (1..25)
    double   star_volume;    // mode 3: parallelepiped volume of triad
    double   last_tilt_mag;  // [rad] |URIC| from last measurement
    double   slotsum;        // for change-detection diagnostic

    Startrack();
    void init()   override;
    void update() override;
    void rpt()    override;

    ACCESS_FN(osk::Vec, URIC)
    ACCESS_FN(int,      star_update_avail)

private:
    void measure_perfect();   // mstar=1
    void measure_noisy();     // mstar=2
    void measure_full();      // mstar=3
    bool select_triad(osk::Vec usii_triad[3], int slot_out[3],
                      double& volume_out);
};

} // namespace rocket6dof

#endif
