//  gps.cpp  --  GPS receiver with simple modes and full Kalman filter
//
//  Modes 0-2: simple verification-friendly behavior (off / perfect / noisy)
//  Mode 3:    Zipfel's 24-satellite constellation + 8-state Kalman filter,
//             ported faithfully from his gps.cpp.
//
//  GPS measurements are only produced at osk::State::stepstart=true so
//  that any INS-state correction takes effect through the RK4 save_
//  cache.  See ins.cpp for the consumer side.

#include "gps.h"
#include "newton.h"
#include "euler.h"
#include "ins.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>
#include <random>

namespace rocket6dof {

namespace {

// ---- GPS satellite ephemeris -----------------------------------------
// 24 satellites in 6 orbital planes (A..F).  Each entry is
// {right_ascension_rad, argument_of_latitude_at_almanac_epoch_rad}.
// Source: Zipfel gps_sv_init() with Yuma Almanac Week 787 (21 Sep 2014).
constexpr double SV_DATA[24][2] = {
    {5.63,  -1.600}, {5.63,   2.115}, {5.63,  -2.309}, {5.63,   0.319},   // A-plane
    {0.40,   1.063}, {0.40,  -1.342}, {0.40,   0.543}, {0.40,   2.874},   // B-plane
    {1.45,   1.705}, {1.45,  -2.841}, {1.45,  -2.321}, {1.45,  -0.640},   // C-plane
    {2.45,   1.941}, {2.45,  -0.147}, {2.45,   1.690}, {2.45,   0.409},   // D-plane
    {3.48,  -0.571}, {3.48,  -2.988}, {3.48,   0.858}, {3.48,   2.705},   // E-plane
    {4.59,  -0.7180},{4.59,   2.666}, {4.59,  -2.977}, {4.59,  -0.2090}   // F-plane
};
constexpr double SV_RADIUS = 26560000.0;
constexpr double SV_INCL   = 0.95986;

inline double sv_angular_rate() {
    return std::sqrt(WGS84_GM / (SV_RADIUS * SV_RADIUS * SV_RADIUS));
}

inline osk::Vec sv_to_inertial(double ra, double au) {
    double si = std::sin(SV_INCL);
    double ci = std::cos(SV_INCL);
    double cra = std::cos(ra), sra = std::sin(ra);
    double cau = std::cos(au), sau = std::sin(au);
    return osk::Vec(
        SV_RADIUS * (cra * cau - sra * sau * ci),
        SV_RADIUS * (sra * cau + cra * sau * ci),
        SV_RADIUS *  sau * si
    );
}
inline osk::Vec sv_velocity_inertial(double ra, double au, double wsi) {
    double si = std::sin(SV_INCL);
    double ci = std::cos(SV_INCL);
    double cra = std::cos(ra), sra = std::sin(ra);
    double cau = std::cos(au), sau = std::sin(au);
    double v = SV_RADIUS * wsi;
    return osk::Vec(
        v * (-sau * cra - cau * sra * ci),
        v * (-sau * sra + cau * cra * ci),
        v *   cau * si
    );
}

inline osk::Mat skew_vec(const osk::Vec& w) {
    return osk::Mat(
         0.0, -w.z,  w.y,
         w.z,  0.0, -w.x,
        -w.y,  w.x,  0.0);
}

} // anon

GPS::GPS()
    : newton(nullptr), euler(nullptr), ins(nullptr),
      mgps(0),
      gps_step(1.0),
      gps_epoch(-1.0e30),
      t_first(0.0),
      rpos(5.0), rvel(0.1),
      noise_seed(42),
      almanac_time(0.0),
      del_rearth(0.0),
      gps_acqtime(2.0),
      ucfreq_noise(0.0),
      ucbias_error(0.0),
      ucfreqm(0.0),
      uctime_cor(3600.0),
      ppos(20.0), pvel(2.0), pclockb(10.0), pclockf(1.0),
      qpos(0.1),  qvel(0.01), qclockb(0.1), qclockf(0.01),
      rpos_kf(5.0), rvel_kf(0.1),
      factp(0.0), factq(0.0), factr(0.0),
      PP(8, 8), FF(8, 8), PHI(8, 8),
      gps_acq(1),
      kf_phase(1),
      meas_count(0),
      slotsum(-1.0),
      SBII_meas(0,0,0), VBII_meas(0,0,0),
      SXH(0,0,0), VXH(0,0,0), CXH(0,0,0),
      gps_update_avail(0),
      last_pos_err(0.0), last_vel_err(0.0),
      gdop(0.0)
{
    for (int i = 0; i < 4; i++) {
        pr_bias[i]  = 0.0;
        pr_noise[i] = 0.0;
        dr_noise[i] = 0.0;
        slot[i]     = 0;
    }
}

void GPS::init() {
    if (initCount == 0) {
        gps_epoch        = -1.0e30;
        gps_update_avail = 0;
        meas_count       = 0;
        last_pos_err     = 0.0;
        last_vel_err     = 0.0;
        SBII_meas        = osk::Vec(0, 0, 0);
        VBII_meas        = osk::Vec(0, 0, 0);
        SXH              = osk::Vec(0, 0, 0);
        VXH              = osk::Vec(0, 0, 0);
        CXH              = osk::Vec(0, 0, 0);
        ucbias_error     = 0.0;
        ucfreqm          = 0.0;
        gps_acq          = 1;
        kf_phase         = 1;
        slotsum          = -1.0;
        gdop             = 0.0;
    }
}

void GPS::update() {
    gps_update_avail = 0;

    if (mgps == 0) return;
    if (!newton)   return;
    if (!osk::State::stepstart) return;

    double t = osk::State::t;
    if (t < t_first) return;

    // ---- Modes 1 and 2: simple measurement ----
    if (mgps == 1 || mgps == 2) {
        if (t < gps_epoch + gps_step) return;
        osk::Vec SBII_truth = newton->SBII;
        osk::Vec VBII_truth = newton->VBII;
        if (mgps == 1) {
            SBII_meas = SBII_truth;
            VBII_meas = VBII_truth;
        } else {
            std::mt19937 eng(noise_seed ^ static_cast<unsigned long>(meas_count));
            std::normal_distribution<double> Npos(0.0, rpos);
            std::normal_distribution<double> Nvel(0.0, rvel);
            osk::Vec npos(Npos(eng), Npos(eng), Npos(eng));
            osk::Vec nvel(Nvel(eng), Nvel(eng), Nvel(eng));
            SBII_meas = SBII_truth + npos;
            VBII_meas = VBII_truth + nvel;
        }
        gps_epoch        = t;
        gps_update_avail = 1;
        meas_count++;
        last_pos_err = (SBII_meas - SBII_truth).mag();
        last_vel_err = (VBII_meas - VBII_truth).mag();
        return;
    }

    // ---- Mode 3: full Kalman filter ----
    if (mgps == 3) {
        if (kf_phase == 1) {
            kf_init();
            gps_epoch = t;
            kf_phase  = 2;
            return;
        }
        kf_extrapolate();
        double dtime_gps = gps_acq ? gps_acqtime : gps_step;
        if ((t - gps_epoch) >= dtime_gps) {
            kf_update();
            // Only mark update_avail if KF didn't fall back to re-acquire
            if (kf_phase == 2) {
                gps_epoch = t;
                gps_acq   = 0;
                gps_update_avail = 1;
                meas_count++;
                last_pos_err = SXH.mag();
                last_vel_err = VXH.mag();
            }
        }
        return;
    }
}

void GPS::rpt() {
    if (osk::State::sample(1.0)) {
        if (mgps == 3) {
            std::printf("GPS  t=%7.3f  mgps=3  meas=%d  GDOP=%.2f  "
                        "|SXH|=%.3f m  |VXH|=%.4f m/s  ucbias=%.4f m\n",
                        osk::State::t, meas_count, gdop,
                        SXH.mag(), VXH.mag(), ucbias_error);
        } else if (mgps != 0) {
            std::printf("GPS  t=%7.3f  mgps=%d  meas=%d  pos_err=%.3f m  vel_err=%.4f m/s  avail=%d\n",
                        osk::State::t, mgps, meas_count,
                        last_pos_err, last_vel_err, gps_update_avail);
        }
    }
}

// ===== Kalman filter sub-functions (mgps=3) defined below ============

void GPS::kf_init() {
    PP.zero();
    for (int i = 0; i < 3; i++) {
        PP.at(i,     i)     = std::pow(ppos * (1.0 + factp), 2);
        PP.at(i + 3, i + 3) = std::pow(pvel * (1.0 + factp), 2);
    }
    PP.at(6, 6) = std::pow(pclockb * (1.0 + factp), 2);
    PP.at(7, 7) = std::pow(pclockf * (1.0 + factp), 2);

    FF.zero();
    FF.at(0, 3) = 1.0;
    FF.at(1, 4) = 1.0;
    FF.at(2, 5) = 1.0;
    FF.at(6, 7) = 1.0;
    FF.at(7, 7) = -1.0 / uctime_cor;

    double dt = osk::State::dt;
    MatN I8 = MatN::identity(8);
    PHI = I8 + FF * dt + (FF * FF) * (0.5 * dt * dt);
}

void GPS::kf_extrapolate() {
    double dt = osk::State::dt;
    double ucfreq_now = ucfreq_noise;
    ucbias_error += (ucfreq_now + ucfreqm) * (dt / 2.0);
    ucfreqm = ucfreq_now;

    MatN QQ(8, 8);
    for (int i = 0; i < 3; i++) {
        QQ.at(i,     i)     = std::pow(qpos * (1.0 + factq), 2);
        QQ.at(i + 3, i + 3) = std::pow(qvel * (1.0 + factq), 2);
    }
    QQ.at(6, 6) = std::pow(qclockb * (1.0 + factq), 2);
    QQ.at(7, 7) = std::pow(qclockf * (1.0 + factq), 2);

    MatN PHIT = PHI.transpose();
    PP = PHI * (PP + QQ * (dt / 2.0)) * PHIT + QQ * (dt / 2.0);
}

void GPS::kf_update() {
    double t = osk::State::t;

    // ---- Select 4 best satellites ----
    osk::Vec ssii_quad[4], vsii_quad[4];
    int      slot_loc[4];
    double   gdop_loc;
    bool ok = quadriga_select(t, ssii_quad, vsii_quad, slot_loc, gdop_loc);
    if (!ok) {
        std::fprintf(stderr, " *** GPS: <4 SVs visible at t=%.3f, re-acquiring ***\n", t);
        kf_phase = 1;
        gps_acq  = 1;
        return;
    }
    gdop = gdop_loc;
    for (int i = 0; i < 4; i++) slot[i] = slot_loc[i];

    double sl_sum = static_cast<double>(slot[0] + slot[1] + slot[2] + slot[3]);
    if (slotsum != sl_sum) {
        slotsum = sl_sum;
        std::printf(" *** GPS Quadriga slots: %d %d %d %d  GDOP=%.2f m ***\n",
                    slot[0], slot[1], slot[2], slot[3], gdop);
    }

    // Get state references
    osk::Vec SBII  = newton->SBII;
    osk::Vec VBII  = newton->VBII;
    osk::Vec SBIIC = ins ? ins->SBIIC : SBII;
    osk::Vec VBIIC = ins ? ins->VBIIC : VBII;
    osk::Vec WBII  = euler ? euler->WBIB : osk::Vec(0, 0, 0);
    osk::Vec WBICI = ins   ? ins->WBICI  : osk::Vec(0, 0, 0);

    MatN ZZ(8, 1);
    MatN HH(8, 8);

    std::mt19937 eng(noise_seed ^ static_cast<unsigned long>(meas_count));

    for (int i = 0; i < 4; i++) {
        osk::Vec SSII = ssii_quad[i];
        osk::Vec VSII = vsii_quad[i];

        osk::Vec SSBI = SSII - SBII;
        double dsb = SSBI.mag();

        std::normal_distribution<double> Npr(0.0, pr_noise[i]);
        std::normal_distribution<double> Ndr(0.0, dr_noise[i]);
        double npr_sample = Npr(eng);
        double ndr_sample = Ndr(eng);

        double dsb_meas = dsb + pr_bias[i] + npr_sample + ucbias_error;

        osk::Vec VSBI = VSII - VBII - skew_vec(WBII) * SSBI;
        osk::Vec USSBI = SSBI * (1.0 / dsb);
        double dvsb = VSBI.dot(USSBI);
        double dvsb_meas = dvsb + ndr_sample + ucfreqm;

        osk::Vec SSBIC = SSII - SBIIC;
        double dsbc = SSBIC.mag();
        osk::Vec VSBIC = VSII - VBIIC - skew_vec(WBICI) * SSBIC;
        osk::Vec USSBIC = SSBIC * (1.0 / dsbc);
        double dvsbc = VSBIC.dot(USSBIC);

        ZZ.at(i,     0) = dsb_meas  - dsbc;
        ZZ.at(i + 4, 0) = dvsb_meas - dvsbc;

        for (int j = 0; j < 3; j++) {
            HH.at(i,     j)     = USSBI[j];
            HH.at(i + 4, j + 3) = USSBI[j] * gps_step;
        }
        HH.at(i,     6) = 1.0;
        HH.at(i + 4, 7) = gps_step;
    }

    MatN RR(8, 8);
    for (int i = 0; i < 4; i++) {
        RR.at(i,     i)     = std::pow(rpos_kf * (1.0 + factr), 2);
        RR.at(i + 4, i + 4) = std::pow(rvel_kf * (1.0 + factr), 2);
    }

    MatN HHT = HH.transpose();
    MatN S   = HH * PP * HHT + RR;
    MatN KK  = PP * HHT * S.inverse();
    MatN XH  = KK * ZZ;
    MatN I8  = MatN::identity(8);
    PP = (I8 - KK * HH) * PP;

    ucbias_error -= XH.at(6, 0);

    SXH = osk::Vec(XH.at(0, 0), XH.at(1, 0), XH.at(2, 0));
    VXH = osk::Vec(XH.at(3, 0), XH.at(4, 0), XH.at(5, 0));
    CXH = osk::Vec(XH.at(6, 0), XH.at(7, 0), 0.0);

    // For INS compatibility with the simple-mode interface: expose the
    // corrected position/velocity as SBII_meas/VBII_meas.  INS reads
    // these whenever gps_update_avail is set, regardless of mode.
    SBII_meas = SBIIC + SXH;
    VBII_meas = VBIIC + VXH;
}

bool GPS::quadriga_select(double t,
                          osk::Vec ssii_quad[4], osk::Vec vsii_quad[4],
                          int slot_out[4], double& gdop_out) {
    double wsi = sv_angular_rate();
    double t_eff = almanac_time + t;

    // Propagate all 24 SVs; determine which are visible.
    osk::Vec sv_r[24], sv_v[24];
    int visible_idx[24];
    int visible_count = 0;

    double r_earth_los = REARTH_SPHERICAL + del_rearth;
    double eps_grazing = std::acos(r_earth_los / SV_RADIUS);

    osk::Vec SBII = newton->SBII;
    double dbi = SBII.mag();

    for (int i = 0; i < 24; i++) {
        double ra  = SV_DATA[i][0];
        double au0 = SV_DATA[i][1];
        double au  = au0 + t_eff * wsi;
        sv_r[i] = sv_to_inertial(ra, au);
        sv_v[i] = sv_velocity_inertial(ra, au, wsi);

        double delta = angle(sv_r[i], SBII);
        bool vis = false;
        if (delta < eps_grazing) {
            vis = true;
        } else {
            double diff = delta - eps_grazing;
            double cdiff = std::cos(diff);
            if (cdiff > 1.0e-12) {
                double rmin = r_earth_los / cdiff;
                if (rmin > 0.0 && rmin < dbi) vis = true;
            }
        }
        if (vis) visible_idx[visible_count++] = i;
    }

    if (visible_count < 4) {
        gdop_out = 1.0e10;
        return false;
    }

    // Try all C(N, 4) combinations of visible SVs; pick the one with
    // minimum GDOP = sqrt(trace((H*H^T)^-1)), where H is the 4x4 matrix
    // whose first 3 columns are -unit-vectors-from-SV-to-user and last
    // column is all ones.
    double best_gdop = 1.0e10;
    int best_quad[4] = {0, 0, 0, 0};

    for (int a = 0; a < visible_count - 3; a++) {
        int ia = visible_idx[a];
        osk::Vec u_a = (SBII - sv_r[ia]) * (1.0 / (SBII - sv_r[ia]).mag());
        for (int b = a + 1; b < visible_count - 2; b++) {
            int ib = visible_idx[b];
            osk::Vec u_b = (SBII - sv_r[ib]) * (1.0 / (SBII - sv_r[ib]).mag());
            for (int c = b + 1; c < visible_count - 1; c++) {
                int ic = visible_idx[c];
                osk::Vec u_c = (SBII - sv_r[ic]) * (1.0 / (SBII - sv_r[ic]).mag());
                for (int d = c + 1; d < visible_count; d++) {
                    int id = visible_idx[d];
                    osk::Vec u_d = (SBII - sv_r[id]) * (1.0 / (SBII - sv_r[id]).mag());

                    // Build 4x4 H
                    MatN H(4, 4);
                    for (int j = 0; j < 3; j++) {
                        H.at(0, j) = u_a[j];
                        H.at(1, j) = u_b[j];
                        H.at(2, j) = u_c[j];
                        H.at(3, j) = u_d[j];
                    }
                    H.at(0, 3) = 1.0;
                    H.at(1, 3) = 1.0;
                    H.at(2, 3) = 1.0;
                    H.at(3, 3) = 1.0;

                    // GDOP = sqrt( trace((H * H^T)^-1) )
                    MatN HT = H.transpose();
                    MatN HHT = H * HT;
                    try {
                        MatN COV = HHT.inverse();
                        double tr = COV.at(0,0) + COV.at(1,1)
                                   + COV.at(2,2) + COV.at(3,3);
                        if (tr <= 0.0) continue;
                        double gd = std::sqrt(tr);
                        if (gd < best_gdop) {
                            best_gdop = gd;
                            best_quad[0] = a;
                            best_quad[1] = b;
                            best_quad[2] = c;
                            best_quad[3] = d;
                        }
                    } catch (...) {
                        continue;   // singular combination, skip
                    }
                }
            }
        }
    }

    if (best_gdop >= 1.0e10) {
        gdop_out = best_gdop;
        return false;
    }

    // Extract the chosen quadriga
    for (int k = 0; k < 4; k++) {
        int vi = visible_idx[best_quad[k]];
        ssii_quad[k] = sv_r[vi];
        vsii_quad[k] = sv_v[vi];
        slot_out[k]  = vi + 1;
    }
    gdop_out = best_gdop;
    return true;
}

} // namespace rocket6dof
