// startrack.cpp -- Star-tracker attitude reference sensor
//
// Implements four modes (0=off, 1=perfect, 2=noisy, 3=full Zipfel triad).
// Like GPS, only produces measurements at osk::State::stepstart so the
// INS attitude correction takes effect through RK4's save_ cache.

#include "startrack.h"
#include "newton.h"
#include "kinematics.h"
#include "ins.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>
#include <random>

namespace rocket6dof {

namespace {

// 25-bright-star catalog (J2000 inertial unit vectors), from Zipfel.
constexpr double STAR_CATALOG[25][3] = {
    {-.179457,  .947482, -.264715}, //  1 Sirius
    {-.062053,  .621699, -.780794}, //  2 Canopus
    {-.395709, -.321011, -.860446}, //  3 Rigil Kent
    {-.787739, -.518432,  .332710}, //  4 Arcturus
    { .119339, -.770884,  .625697}, //  5 Vega
    { .205167,  .969392, -.134851}, //  6 Rigel
    { .141589,  .680716,  .718733}, //  7 Capella
    {-.407741,  .908325,  .093239}, //  8 Procyon
    { .504096,  .224174, -.834046}, //  9 Achernar
    {-.434428, -.251576, -.864859}, // 10 Hadar
    { .449879, -.880088,  .151836}, // 11 Altair
    { .355451,  .890944,  .282620}, // 12 Aldebaran
    {-.479412, -.049965, -.876167}, // 13 Acrux
    { .032446,  .991140,  .128796}, // 14 Betelgeuse
    {-.357911, -.827083, -.433397}, // 15 Antares
    {-.923975, -.348220, -.158158}, // 16 Spica
    {-.380630,  .795326,  .471781}, // 17 Pollux
    { .846646, -.247176, -.471268}, // 18 Fomalhaut
    { .453017, -.541322,  .708340}, // 19 Deneb
    {-.511331, -.101246, -.853399}, // 20 Mimosa
    {-.858250,  .467448,  .211893}, // 21 Regulus
    {-.217999,  .863108, -.455545}, // 22 Adhara
    { .161912,  .980685,  .109734}, // 23 Bellatrix
    {-.103643, -.792588, -.600885}, // 24 Shaula
    { .140796,  .866902,  .478181}  // 25 El Nath
};
constexpr const char* STAR_NAMES[25] = {
    "Sirius", "Canopus", "RigilKent", "Arcturus", "Vega",
    "Rigel", "Capella", "Procyon", "Achernar", "Hadar",
    "Altair", "Aldebaran", "Acrux", "Betelgeuse", "Antares",
    "Spica", "Pollux", "Fomalhaut", "Deneb", "Mimosa",
    "Regulus", "Adhara", "Bellatrix", "Shaula", "ElNath"
};

inline osk::Mat skew_vec(const osk::Vec& w) {
    return osk::Mat(
         0.0, -w.z,  w.y,
         w.z,  0.0, -w.x,
        -w.y,  w.x,  0.0);
}

inline void cart_to_pol(const osk::Vec& u, double& az, double& el) {
    az = std::atan2(u.y, u.x);
    double r = u.mag();
    if (r < 1.0e-12) { el = 0.0; return; }
    el = std::asin(-u.z / r);
}
inline osk::Vec pol_to_cart(double az, double el) {
    double ce = std::cos(el);
    return osk::Vec(ce * std::cos(az), ce * std::sin(az), -std::sin(el));
}

// 3x3 inverse via cofactor expansion
osk::Mat mat3_inverse(const osk::Mat& M) {
    double a = M[0][0], b = M[0][1], c = M[0][2];
    double d = M[1][0], e = M[1][1], f = M[1][2];
    double g = M[2][0], h = M[2][1], i = M[2][2];
    double det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    if (std::fabs(det) < 1.0e-12) {
        return osk::Mat(1, 0, 0, 0, 1, 0, 0, 0, 1);
    }
    double inv = 1.0 / det;
    return osk::Mat(
         (e*i - f*h) * inv, -(b*i - c*h) * inv,  (b*f - c*e) * inv,
        -(d*i - f*g) * inv,  (a*i - c*g) * inv, -(a*f - c*d) * inv,
         (d*h - e*g) * inv, -(a*h - b*g) * inv,  (a*e - b*d) * inv);
}

} // anon

// STUB_MARKER

Startrack::Startrack()
    : newton(nullptr), kin(nullptr), ins(nullptr),
      mstar(0),
      startrack_step(1.0),
      starfix_epoch(-1.0e30),
      t_first(0.0),
      startrack_alt(0.0),
      star_acqtime(2.0),
      star_acq(1),
      tilt_noise(1.0e-4),
      star_el_min(20.0),
      noise_seed(42),
      URIC(0, 0, 0),
      star_update_avail(0),
      meas_count(0),
      star_volume(0.0),
      last_tilt_mag(0.0),
      slotsum(-1.0)
{
    for (int i = 0; i < 3; i++) {
        az_bias[i]  = 0.0;
        az_noise[i] = 0.0;
        el_bias[i]  = 0.0;
        el_noise[i] = 0.0;
        triad[i]    = 0;
    }
}

void Startrack::init() {
    if (initCount == 0) {
        starfix_epoch     = -1.0e30;
        star_update_avail = 0;
        meas_count        = 0;
        star_acq          = 1;
        URIC              = osk::Vec(0, 0, 0);
        last_tilt_mag     = 0.0;
        slotsum           = -1.0;
        for (int i = 0; i < 3; i++) triad[i] = 0;
    }
}

void Startrack::update() {
    star_update_avail = 0;
    if (mstar == 0) return;
    if (!newton)    return;
    if (!osk::State::stepstart) return;

    double t = osk::State::t;
    if (t < t_first) return;
    // Altitude gate: stars are only useful above the atmosphere.
    // Applies to all measurement modes (perfect, noisy, full).
    if (newton->alt < startrack_alt) return;

    double dtime_star = star_acq ? star_acqtime : startrack_step;
    if (t < starfix_epoch + dtime_star) return;

    if      (mstar == 1) measure_perfect();
    else if (mstar == 2) measure_noisy();
    else if (mstar == 3) measure_full();
    else                 return;

    starfix_epoch     = t;
    star_acq          = 0;
    star_update_avail = 1;
    meas_count++;
    last_tilt_mag = URIC.mag();
}

void Startrack::rpt() {
    if (osk::State::sample(1.0)) {
        if (mstar == 3) {
            std::printf("Star t=%7.3f  mstar=3  meas=%d  vol=%.3f  triad={%d,%d,%d}  |URIC|=%.2e rad  avail=%d\n",
                        osk::State::t, meas_count, star_volume,
                        triad[0], triad[1], triad[2],
                        last_tilt_mag, star_update_avail);
        } else if (mstar != 0) {
            std::printf("Star t=%7.3f  mstar=%d  meas=%d  |URIC|=%.3e rad  avail=%d\n",
                        osk::State::t, mstar, meas_count,
                        last_tilt_mag, star_update_avail);
        }
    }
}

// ---- Mode 1: perfect attitude --------------------------------------
// URIC = small-angle vector from TBIC to TBI_true.
// Extracted from skew-symmetric part of (TBI_true * TBIC^T - I).
void Startrack::measure_perfect() {
    if (!kin || !ins) {
        URIC = osk::Vec(0, 0, 0);
        return;
    }
    osk::Mat E = kin->TBI * ins->TBIC.transpose()
               - osk::Mat(1, 0, 0, 0, 1, 0, 0, 0, 1);
    URIC = osk::Vec(E[2][1], E[0][2], E[1][0]);
}

// ---- Mode 2: noisy attitude ----------------------------------------
void Startrack::measure_noisy() {
    std::mt19937 eng(noise_seed ^ static_cast<unsigned long>(meas_count));
    std::normal_distribution<double> N(0.0, tilt_noise);
    URIC = osk::Vec(N(eng), N(eng), N(eng));
}

// STUB_MODE3_MARKER

// ---- Mode 3: full Zipfel triad-based star tracking ------------------
// Algorithm:
//   1. Pick 3 visible stars maximizing |u1.(u2 x u3)| (parallelepiped vol)
//   2. For each star: compute true (az, el) in body via TBI_true
//   3. Add per-star bias + noise to (az, el)
//   4. Convert measured (az, el) back to inertial via TBIC
//   5. RDIFF = TRIAD_MEAS * TRIAD_TRUE^-1
//   6. URIC.x = RDIFF[2][1]; URIC.y = RDIFF[0][2]; URIC.z = RDIFF[1][0]
void Startrack::measure_full() {
    osk::Vec usii_triad[3];
    int      slot_loc[3];
    double   vol = 0.0;
    if (!select_triad(usii_triad, slot_loc, vol)) {
        URIC = osk::Vec(0, 0, 0);
        return;
    }
    star_volume = vol;
    for (int i = 0; i < 3; i++) triad[i] = slot_loc[i];

    double sl_sum = static_cast<double>(slot_loc[0] + slot_loc[1] + slot_loc[2]);
    if (slotsum != sl_sum) {
        slotsum = sl_sum;
        std::printf(" *** Star triad: %s, %s, %s  volume=%.4f ***\n",
                    STAR_NAMES[slot_loc[0]-1],
                    STAR_NAMES[slot_loc[1]-1],
                    STAR_NAMES[slot_loc[2]-1], vol);
    }

    if (!kin || !ins) {
        URIC = osk::Vec(0, 0, 0);
        return;
    }
    osk::Mat TBI_true = kin->TBI;
    osk::Mat TBIC_est = ins->TBIC;

    std::mt19937 eng(noise_seed ^ static_cast<unsigned long>(meas_count));

    osk::Mat TRIAD_TRUE(0,0,0, 0,0,0, 0,0,0);
    osk::Mat TRIAD_MEAS(0,0,0, 0,0,0, 0,0,0);

    for (int i = 0; i < 3; i++) {
        osk::Vec USII = usii_triad[i];
        osk::Vec USIB = TBI_true * USII;
        double az, el;
        cart_to_pol(USIB, az, el);
        std::normal_distribution<double> Naz(0.0, az_noise[i]);
        std::normal_distribution<double> Nel(0.0, el_noise[i]);
        double az_meas = az + az_bias[i] + Naz(eng);
        double el_meas = el + el_bias[i] + Nel(eng);
        osk::Vec USIBM = pol_to_cart(az_meas, el_meas);
        osk::Vec USIIM = TBIC_est.transpose() * USIBM;
        for (int j = 0; j < 3; j++) {
            TRIAD_TRUE[j][i] = USII[j];
            TRIAD_MEAS[j][i] = USIIM[j];
        }
    }

    osk::Mat TRIAD_TRUE_INV = mat3_inverse(TRIAD_TRUE);
    osk::Mat RDIFF = TRIAD_MEAS * TRIAD_TRUE_INV;
    URIC = osk::Vec(RDIFF[2][1], RDIFF[0][2], RDIFF[1][0]);
}

// STUB_TRIAD_MARKER

// ---- Triad selection: max parallelepiped volume ---------------------
bool Startrack::select_triad(osk::Vec usii_triad[3], int slot_out[3],
                             double& volume_out) {
    osk::Vec UBII_pos = newton->SBII;
    double r = UBII_pos.mag();
    if (r < 1.0) {
        volume_out = 0.0;
        return false;
    }
    osk::Vec UBII = UBII_pos * (1.0 / r);

    double el_min_rad = star_el_min * (osk::PI / 180.0);

    int visible_idx[25];
    int visible_count = 0;
    for (int i = 0; i < 25; i++) {
        osk::Vec USII(STAR_CATALOG[i][0], STAR_CATALOG[i][1], STAR_CATALOG[i][2]);
        double el = std::asin(USII.dot(UBII));
        if (el > el_min_rad) {
            visible_idx[visible_count++] = i;
        }
    }

    if (visible_count < 3) {
        volume_out = 0.0;
        return false;
    }

    double best_vol = -1.0;
    int    best[3] = {0, 0, 0};
    for (int a = 0; a < visible_count - 2; a++) {
        int ia = visible_idx[a];
        osk::Vec ua(STAR_CATALOG[ia][0], STAR_CATALOG[ia][1], STAR_CATALOG[ia][2]);
        for (int b = a + 1; b < visible_count - 1; b++) {
            int ib = visible_idx[b];
            osk::Vec ub(STAR_CATALOG[ib][0], STAR_CATALOG[ib][1], STAR_CATALOG[ib][2]);
            for (int c = b + 1; c < visible_count; c++) {
                int ic = visible_idx[c];
                osk::Vec uc(STAR_CATALOG[ic][0], STAR_CATALOG[ic][1], STAR_CATALOG[ic][2]);
                double v = std::fabs(ua.dot(ub.cross(uc)));
                if (v > best_vol) {
                    best_vol = v;
                    best[0] = a;
                    best[1] = b;
                    best[2] = c;
                }
            }
        }
    }

    for (int k = 0; k < 3; k++) {
        int idx = visible_idx[best[k]];
        usii_triad[k] = osk::Vec(STAR_CATALOG[idx][0],
                                  STAR_CATALOG[idx][1],
                                  STAR_CATALOG[idx][2]);
        slot_out[k] = idx + 1;
    }
    volume_out = best_vol;
    return true;
}

} // namespace rocket6dof
