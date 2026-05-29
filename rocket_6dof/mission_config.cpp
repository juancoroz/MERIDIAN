//  mission_config.cpp  --  Apply JSON config to rocket_6dof blocks
//
//  Each block has its own helper that walks its config section and
//  applies the listed keys to the block.  Missing keys keep the
//  block's existing value (from its constructor default or previous
//  assignment).
//
//  Vector members (e.g. ins.bias_accel) expect a 3-element JSON array
//  [x, y, z].  If the value is not an array or has wrong arity, the
//  member is left unchanged.
//
//  Adding a new field to a block:
//    1. Add the public member in the block's .h file.
//    2. Add one line here in the relevant load_* helper:
//         if (s["new_field"].exists()) blk->new_field = s["new_field"].asNumber();

#include "mission_config.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "propulsion.h"
#include "aerodynamics.h"
#include "tvc.h"
#include "rcs.h"
#include "control.h"
#include "guidance.h"
#include "ins.h"
#include "gps.h"
#include "startrack.h"
#include "intercept.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"

namespace rocket6dof {

using json::Value;

// Small helpers to make the per-block functions short and readable.
// `sd(v, def)` returns the value as a double if it's a number, else `def`.
// Crucially, if the key didn't exist at all, the Value is MISSING and
// asNumber() returns the default, so a block keeps its current setting.
namespace {

void set_d(double& target, const Value& v) {
    if (v.isNumber()) target = v.asNumber();
}
void set_i(int& target, const Value& v) {
    if (v.isNumber()) target = v.asInt();
}
void set_b(bool& target, const Value& v) {
    if (v.isBool()) target = v.asBool();
}
void set_s(std::string& target, const Value& v) {
    if (v.isString()) target = v.asString();
}
void set_vec(osk::Vec& target, const Value& v) {
    if (v.isArray() && v.size() == 3) {
        target = osk::Vec(v[0].asNumber(target.x),
                          v[1].asNumber(target.y),
                          v[2].asNumber(target.z));
    }
}

//  Propulsion
void load_propulsion(const Value& s, Propulsion* p) {
    if (!s.isObject() || !p) return;
    set_i(p->mprop,          s["mprop"]);
    set_d(p->vmass0,         s["vmass0"]);
    set_d(p->fmass0,         s["fmass0"]);
    set_d(p->spi,             s["spi"]);
    set_d(p->fuel_flow_rate, s["fuel_flow_rate"]);
    set_d(p->moi_roll_0,     s["moi_roll_0"]);
    set_d(p->moi_roll_1,     s["moi_roll_1"]);
    set_d(p->moi_trans_0,    s["moi_trans_0"]);
    set_d(p->moi_trans_1,    s["moi_trans_1"]);
    set_d(p->xcg_0,          s["xcg_0"]);
    set_d(p->xcg_1,          s["xcg_1"]);

    // Multi-stage configuration via "stages" array of per-stage objects.
    // The kernel propulsion struct has fixed-size arrays sized to
    // Propulsion::MAX_STAGES (default 4, overridable at compile time via
    // -DROCKET6DOF_MAX_STAGES=N).  This loader truncates any extra
    // entries past that limit.
    // The guidance arrays (BOTN, BURNTN, TAUN, etc.) are still sized
    // for 3 stages, so configs with 4+ stages AND LTG guidance enabled
    // (mguide=1) get a config-time error from load_guidance below.
    // mguide=0 (disabled) and mguide=2 (attitude program) work fine
    // with any stage count.
    if (s["stages"].isArray()) {
        const auto& stages = s["stages"].arr();
        int requested_n = static_cast<int>(stages.size());
        int n = requested_n;
        if (n > Propulsion::MAX_STAGES) {
            n = Propulsion::MAX_STAGES;
            std::fprintf(stderr,
                "WARN: propulsion.stages has %d entries; truncating to %d "
                "(MAX_STAGES limit).  Rebuild with "
                "-DROCKET6DOF_MAX_STAGES=N to lift this cap.\n",
                requested_n, Propulsion::MAX_STAGES);
        }
        p->num_stages = n;
        for (int i = 0; i < n; i++) {
            const Value& st = stages[i];
            set_d(p->vmass0_stage[i],     st["vmass0"]);
            set_d(p->fmass0_stage[i],     st["fmass0"]);
            set_d(p->spi_stage[i],        st["spi"]);
            set_d(p->mdot_stage[i],       st["fuel_flow_rate"]);
            set_d(p->xcg_0_stage[i],      st["xcg_0"]);
            set_d(p->xcg_1_stage[i],      st["xcg_1"]);
            set_d(p->moi_roll_0_stage[i], st["moi_roll_0"]);
            set_d(p->moi_roll_1_stage[i], st["moi_roll_1"]);
            set_d(p->moi_trans_0_stage[i],st["moi_trans_0"]);
            set_d(p->moi_trans_1_stage[i],st["moi_trans_1"]);
            set_d(p->dry_mass_dropped[i], st["dry_mass_dropped"]);
            set_d(p->coast_to_next[i],    st["coast_to_next"]);
        }
    }
}

//  Aerodynamics
void load_aerodynamics(const Value& s, Aerodynamics* a) {
    if (!s.isObject() || !a) return;
    set_i(a->maero,       s["maero"]);
    set_d(a->refa,         s["refa"]);
    set_d(a->refd,         s["refd"]);
    set_d(a->xcg_ref,      s["xcg_ref"]);
    set_s(a->aero_file,    s["aero_file"]);
    set_s(a->tag_ca0,      s["tag_ca0"]);
    set_s(a->tag_caa,      s["tag_caa"]);
    set_s(a->tag_ca0b,     s["tag_ca0b"]);
    set_s(a->tag_cn0,      s["tag_cn0"]);
    set_s(a->tag_clm0,     s["tag_clm0"]);
    set_s(a->tag_clmq,     s["tag_clmq"]);
}

//  TVC
void load_tvc(const Value& s, TVC* t) {
    if (!s.isObject() || !t) return;
    set_i(t->mtvc,    s["mtvc"]);
    set_d(t->gtvc,    s["gtvc"]);
    set_d(t->parm,    s["parm"]);
    set_d(t->del_max, s["del_max"]);
}

//  RCS
void load_rcs(const Value& s, RCS* r) {
    if (!s.isObject() || !r) return;
    set_i(r->mrcs_moment,    s["mrcs_moment"]);
    set_i(r->mrcs_force,     s["mrcs_force"]);
    set_d(r->rcs_freq,       s["rcs_freq"]);
    set_d(r->rcs_zeta,       s["rcs_zeta"]);
    set_d(r->dead_zone,      s["dead_zone"]);
    set_d(r->hysteresis,     s["hysteresis"]);
    set_d(r->rcs_tau,        s["rcs_tau"]);
    set_d(r->roll_mom_max,   s["roll_mom_max"]);
    set_d(r->pitch_mom_max,  s["pitch_mom_max"]);
    set_d(r->yaw_mom_max,    s["yaw_mom_max"]);
    set_d(r->side_force_max, s["side_force_max"]);
    set_d(r->acc_gain,       s["acc_gain"]);
    set_d(r->phibdcomx,      s["phibdcomx"]);
    set_d(r->thtbdcomx,      s["thtbdcomx"]);
    set_d(r->psibdcomx,      s["psibdcomx"]);
    set_d(r->alphacomx,      s["alphacomx"]);
    set_d(r->betacomx,       s["betacomx"]);
    set_d(r->aycomx,         s["aycomx"]);
    set_d(r->azcomx,         s["azcomx"]);
}

//  Control
void load_control(const Value& s, Control* c) {
    if (!s.isObject() || !c) return;
    set_i(c->maut,                s["maut"]);
    set_d(c->vac_rate_damp,       s["vac_rate_damp"]);
    set_d(c->vac_max_gain,        s["vac_max_gain"]);
    set_d(c->tau_att,             s["tau_att"]);
    set_d(c->q_max,                s["q_max"]);
    set_d(c->theta_com_inertial,  s["theta_com_inertial"]);
    set_d(c->psi_com_inertial,    s["psi_com_inertial"]);
    set_d(c->delimx,              s["delimx"]);
    set_d(c->drlimx,              s["drlimx"]);
    set_d(c->gnmax,               s["gnmax"]);
    set_d(c->gymax,                s["gymax"]);
    set_d(c->ancomx,              s["ancomx"]);
    set_d(c->alcomx,              s["alcomx"]);
    set_d(c->thrust_loc,          s["thrust_loc"]);
    set_d(c->load_relief_q_threshold, s["load_relief_q_threshold"]);
    set_d(c->load_relief_q_width,     s["load_relief_q_width"]);
    set_d(c->accel_relief_gain,         s["accel_relief_gain"]);
    set_d(c->accel_relief_q_threshold,  s["accel_relief_q_threshold"]);
    set_d(c->accel_relief_q_width,      s["accel_relief_q_width"]);
    set_d(c->accel_relief_tau,          s["accel_relief_tau"]);
    set_d(c->load_relief_tau_factor,    s["load_relief_tau_factor"]);
    set_d(c->waclp,               s["waclp"]);
    set_d(c->zaclp,               s["zaclp"]);
    set_d(c->paclp,               s["paclp"]);
    set_d(c->factwaclp,           s["factwaclp"]);
    set_d(c->wacly,               s["wacly"]);
    set_d(c->zacly,               s["zacly"]);
    set_d(c->pacly,               s["pacly"]);
    set_d(c->factwacly,           s["factwacly"]);
}

//  Guidance
void load_guidance(const Value& s, Guidance* g) {
    if (!s.isObject() || !g) return;
    set_i(g->mguide,             s["mguide"]);
    set_i(g->ltg_drives_control, s["ltg_drives_control"]);
    set_d(g->phase2_yaw_sign,    s["phase2_yaw_sign"]);
    // mguide=1 (accel program)
    set_d(g->t_pitch_start,      s["t_pitch_start"]);
    set_d(g->t_pitch_end,        s["t_pitch_end"]);
    set_d(g->ancomx_program,     s["ancomx_program"]);
    set_d(g->alcomx_program,     s["alcomx_program"]);
    // mguide=2 (attitude program)
    set_d(g->t_att_start,        s["t_att_start"]);
    set_d(g->t_att_end,          s["t_att_end"]);
    set_d(g->theta_com_start,    s["theta_com_start"]);
    set_d(g->theta_com_end,      s["theta_com_end"]);
    set_d(g->psi_com_start,      s["psi_com_start"]);
    set_d(g->psi_com_end,        s["psi_com_end"]);
    // LTG
    set_d(g->ltg_step,           s["ltg_step"]);
    set_d(g->dbi_desired,        s["dbi_desired"]);
    set_d(g->dvbi_desired,       s["dvbi_desired"]);
    set_d(g->thtvdx_desired,     s["thtvdx_desired"]);
    set_i(g->num_stages,         s["num_stages"]);
    set_d(g->delay_ignition,     s["delay_ignition"]);
    set_d(g->amin,               s["amin"]);

    // Load per-stage LTG arrays.  Three input forms accepted:
    //   1) Array form (preferred):    "char_time": [120, 80, 60, 40]
    //   2) Stages-block form:         "stages": [ {"char_time":120,...}, ... ]
    //   3) Legacy scalar form:        "char_time1": 120, "char_time2": 80, ...
    // Last loaded wins on collision; the rule of thumb is "use the array
    // form for new configs; legacy scalars work but only address stages 1-3."
    const Value& ct_arr = s["char_time"];
    const Value& ve_arr = s["exhaust_vel"];
    const Value& be_arr = s["burnout_epoch"];
    if (ct_arr.isArray()) {
        int n = static_cast<int>(ct_arr.size());
        if (n > Guidance::MAX_LTG_STAGES) n = Guidance::MAX_LTG_STAGES;
        for (int k = 0; k < n; ++k) set_d(g->char_time[k], ct_arr[k]);
    }
    if (ve_arr.isArray()) {
        int n = static_cast<int>(ve_arr.size());
        if (n > Guidance::MAX_LTG_STAGES) n = Guidance::MAX_LTG_STAGES;
        for (int k = 0; k < n; ++k) set_d(g->exhaust_vel[k], ve_arr[k]);
    }
    if (be_arr.isArray()) {
        int n = static_cast<int>(be_arr.size());
        if (n > Guidance::MAX_LTG_STAGES) n = Guidance::MAX_LTG_STAGES;
        for (int k = 0; k < n; ++k) set_d(g->burnout_epoch[k], be_arr[k]);
    }

    // Legacy scalar accessors (back-compat).  Stage indices 1/2/3 map
    // to array indices 0/1/2.  These overwrite array values if both
    // are present in the same config -- last set wins -- so don't mix
    // in the same config.
    set_d(g->char_time[0],       s["char_time1"]);
    set_d(g->char_time[1],       s["char_time2"]);
    set_d(g->char_time[2],       s["char_time3"]);
    set_d(g->exhaust_vel[0],     s["exhaust_vel1"]);
    set_d(g->exhaust_vel[1],     s["exhaust_vel2"]);
    set_d(g->exhaust_vel[2],     s["exhaust_vel3"]);
    set_d(g->burnout_epoch[0],   s["burnout_epoch1"]);
    set_d(g->burnout_epoch[1],   s["burnout_epoch2"]);
    set_d(g->burnout_epoch[2],   s["burnout_epoch3"]);
    set_d(g->lamd_limit,         s["lamd_limit"]);

    // Cap guard: num_stages must not exceed MAX_LTG_STAGES IF LTG (mguide=5)
    // is enabled.  Other guidance modes (0=off, 1=accel program, 2=attitude
    // program) don't use the per-stage arrays and are safe with any
    // num_stages value.
    if (g->mguide == 5 && g->num_stages > Guidance::MAX_LTG_STAGES) {
        std::fprintf(stderr,
            "ERROR: guidance.mguide=5 (LTG) with num_stages=%d exceeds "
            "MAX_LTG_STAGES=%d.\n"
            "       Rebuild with -DROCKET6DOF_MAX_STAGES=N to lift this cap.\n",
            g->num_stages, Guidance::MAX_LTG_STAGES);
        std::exit(1);
    }
}

//  INS
void load_ins(const Value& s, INS* ins) {
    if (!s.isObject() || !ins) return;
    set_i(ins->mins,        s["mins"]);
    set_vec(ins->bias_accel, s["bias_accel"]);
    set_vec(ins->bias_gyro,  s["bias_gyro"]);
}

//  GPS
void load_gps(const Value& s, GPS* g) {
    if (!s.isObject() || !g) return;
    set_i(g->mgps,         s["mgps"]);
    set_d(g->gps_step,     s["gps_step"]);
    set_d(g->t_first,      s["t_first"]);
    set_d(g->rpos,         s["rpos"]);
    set_d(g->rvel,         s["rvel"]);
    set_d(g->almanac_time, s["almanac_time"]);
    set_d(g->del_rearth,   s["del_rearth"]);
    set_d(g->gps_acqtime,  s["gps_acqtime"]);
    if (s["noise_seed"].isNumber()) {
        g->noise_seed = static_cast<unsigned long>(s["noise_seed"].asNumber());
    }
}

//  Startrack
void load_startrack(const Value& s, Startrack* st) {
    if (!s.isObject() || !st) return;
    set_i(st->mstar,           s["mstar"]);
    set_d(st->startrack_step,  s["startrack_step"]);
    set_d(st->t_first,         s["t_first"]);
    set_d(st->startrack_alt,   s["startrack_alt"]);
    set_d(st->star_acqtime,    s["star_acqtime"]);
    set_d(st->tilt_noise,      s["tilt_noise"]);
    set_d(st->star_el_min,     s["star_el_min"]);
}

//  Intercept
void load_intercept(const Value& s, Intercept* ic) {
    if (!s.isObject() || !ic) return;
    set_b(ic->check_ground_impact, s["check_ground_impact"]);
    set_b(ic->check_apogee,        s["check_apogee"]);
    set_b(ic->check_alt_max,       s["check_alt_max"]);
    set_b(ic->check_aero_decay,    s["check_aero_decay"]);
    set_b(ic->check_time_max,      s["check_time_max"]);
    set_b(ic->check_tumble,        s["check_tumble"]);
    set_d(ic->alt_max,             s["alt_max"]);
    set_d(ic->trmach,              s["trmach"]);
    set_d(ic->trdynm,              s["trdynm"]);
    set_d(ic->t_max,                s["t_max"]);
    set_d(ic->rate_max,            s["rate_max"]);
    set_d(ic->t_no_check_until,    s["t_no_check_until"]);
}

//  Launch initial conditions (Newton + Kinematics + Euler)
void load_launch(const Value& s, Newton* n, Kinematics* k, Euler* e) {
    if (!s.isObject()) return;
    if (n) {
        set_d(n->lonx0,    s["lonx0"]);
        set_d(n->latx0,    s["latx0"]);
        set_d(n->alt0,     s["alt0"]);
        set_d(n->dvbe0,    s["dvbe0"]);
        set_d(n->psivdx0,  s["psivdx0"]);
        set_d(n->thtvdx0,  s["thtvdx0"]);
    }
    if (k) {
        set_d(k->psibdx0,  s["psibdx0"]);
        set_d(k->thtbdx0,  s["thtbdx0"]);
        set_d(k->phibdx0,  s["phibdx0"]);
    }
    if (e) {
        set_d(e->ppx0,     s["ppx0"]);
        set_d(e->qqx0,     s["qqx0"]);
        set_d(e->rrx0,     s["rrx0"]);
    }
}

//  External forces
void load_forces(const Value& s, Forces* f) {
    if (!s.isObject() || !f) return;
    set_vec(f->FAPB_ext, s["FAPB_ext"]);
    set_vec(f->FMB_ext,  s["FMB_ext"]);
}

} // anon

//  Top-level entry point
int apply_config(const Value& cfg,
                 Propulsion*   prop,
                 Aerodynamics* aero,
                 TVC*          tvc,
                 RCS*          rcs,
                 Control*      ctrl,
                 Guidance*     guid,
                 INS*          ins,
                 GPS*          gps,
                 Startrack*    star,
                 Intercept*    icpt,
                 Newton*       newt,
                 Kinematics*   kin,
                 Euler*        eul,
                 Forces*       forc)
{
    if (!cfg.isObject()) return 1;

    load_propulsion  (cfg["propulsion"],   prop);
    load_aerodynamics(cfg["aerodynamics"], aero);
    load_tvc         (cfg["tvc"],          tvc);
    load_rcs         (cfg["rcs"],          rcs);

    // Cross-check: long inter-stage coast without RCS is a common
    // foot-gun.  During coast, thrust is zero -> TVC has nothing to
    // gimbal -> attitude is controlled only by RCS.  If RCS is also
    // off (mrcs_moment=0) and the vehicle is at high q_dyn, any
    // small angular rate at start of coast grows under aero moments
    // until the vehicle tumbles, which then ruins the next stage's
    // burn.  This was observed in mission_3stage.json with the
    // default 1-second coast; the trajectory loses ~85% of its
    // Tsiolkovsky dv because the vehicle is sideways when stage 2
    // ignites.
    //
    // The check is heuristic: we warn (not abort) when any stage
    // has coast > 0.5s and RCS-moment is off.  Users running in
    // vacuum or with very brief coasts will be fine, but the warning
    // surfaces the risk for new users assembling configs.
    if (rcs && prop && rcs->mrcs_moment == 0 && prop->num_stages > 1) {
        for (int s = 0; s < prop->num_stages - 1; ++s) {
            if (prop->coast_to_next[s] > 0.5) {
                std::fprintf(stderr,
                    "WARN: propulsion.stages[%d].coast_to_next = %.2fs with "
                    "rcs.mrcs_moment = 0.\n"
                    "      Inter-stage coast >0.5s without RCS attitude "
                    "control can produce\n"
                    "      attitude divergence under aerodynamic moments at "
                    "high dynamic\n"
                    "      pressure.  Set mrcs_moment=1, or shorten the "
                    "coast, or accept\n"
                    "      that the vehicle may tumble during this coast "
                    "phase.\n",
                    s, prop->coast_to_next[s]);
            }
        }
    }

    load_control     (cfg["control"],      ctrl);

    // Cross-check: control.thrust_loc and tvc.parm must agree.
    // Both fields name the same physical point (the gimbal location
    // measured from vehicle nose, in meters).  TVC uses parm to
    // compute the moment it actually applies; Control uses
    // thrust_loc to compute the moment its autopilot expects from a
    // commanded deflection.  If they disagree, the autopilot's
    // model differs from the plant: closed-loop response will be
    // wrong-magnitude or wrong-sign and the vehicle may tumble.
    // (This was the root cause of the launcher's "untransferable
    // autopilot" symptom for several sessions before being traced.)
    //
    // The check fires only when TVC is active (mtvc != 0) and at
    // least one of the two values has been overridden from its
    // default 5.0 -- otherwise both are at default and trivially
    // agree.  We warn rather than abort because the user might
    // intentionally model a TVC/autopilot mismatch (e.g. for fault
    // injection studies).
    if (tvc && ctrl && tvc->mtvc != 0) {
        double diff = std::fabs(tvc->parm - ctrl->thrust_loc);
        // 1 cm tolerance: large enough to absorb floating-point
        // re-encoding through JSON; small enough that any real
        // configuration difference exceeds it.
        if (diff > 0.01) {
            std::fprintf(stderr,
                "WARN: tvc.parm = %.4f m differs from control.thrust_loc = "
                "%.4f m.\n"
                "      Both should name the same physical gimbal location "
                "(meters from nose).\n"
                "      TVC applies moment as if gimbal is at parm; "
                "Control's autopilot\n"
                "      expects gimbal at thrust_loc.  Disagreement produces "
                "wrong-magnitude or\n"
                "      wrong-sign closed-loop response.  Set both to the "
                "same value unless\n"
                "      you are intentionally modeling a plant/controller "
                "mismatch.\n",
                tvc->parm, ctrl->thrust_loc);
        }
    }

    // Semantic check: thrust_loc must be aft of the vehicle CG
    // throughout the burn, i.e. greater than both xcg_0 (full) and
    // xcg_1 (empty).  When the CG sits forward of the gimbal,
    // arm = thrust_loc - xcg > 0 and a positive gimbal deflection
    // produces a positive moment.  When the CG sits aft of the
    // gimbal, the arm goes negative -- the autopilot's sign
    // convention inverts and closed-loop response goes the wrong
    // direction.  This is a different failure mode from the parm-
    // vs-thrust_loc mismatch above: there, plant and controller
    // disagree; here, both agree but on a physically nonsensical
    // configuration.
    //
    // The arm in control.cpp clamps to 0.5m if the formula returns
    // a smaller value, which masks the sign problem but leaves the
    // magnitude wrong by orders of magnitude.  Either way the
    // autopilot misbehaves; this check surfaces the upstream config
    // error before run time.
    //
    // Checked per-stage: each stage has its own xcg_0/xcg_1.  Stages
    // dropped early may have CG well aft of the gimbal; stages that
    // ignite later have their own xcg values that need to be aft of
    // their own thrust attachment.  But the simplifying assumption
    // here is that the single control.thrust_loc applies across all
    // stages (the autopilot can't switch gimbal locations mid-burn).
    // So the check requires thrust_loc > max(xcg) across all stages.
    if (tvc && ctrl && prop && tvc->mtvc != 0) {
        double max_xcg = 0.0;
        int max_stage = 0;
        for (int s = 0; s < prop->num_stages; ++s) {
            double xcg_0 = prop->xcg_0_stage[s];
            double xcg_1 = prop->xcg_1_stage[s];
            double sm = (xcg_0 > xcg_1) ? xcg_0 : xcg_1;
            if (sm > max_xcg) { max_xcg = sm; max_stage = s; }
        }
        // 0.1m margin: prevents the arm clamping to 0.5m and ensures
        // the autopilot has some real moment authority.  Set tight
        // (10cm) because vehicles with thrust attachment essentially
        // at the CG are physically pathological -- the warning is
        // helpful in those cases, even if the user knows.
        if (max_xcg + 0.1 > ctrl->thrust_loc) {
            std::fprintf(stderr,
                "WARN: control.thrust_loc = %.4f m is not aft of vehicle CG "
                "(max xcg = %.4f m on stage %d).\n"
                "      The gimbal moment arm = thrust_loc - xcg would be "
                "%.4f m (negative or near-zero).\n"
                "      Autopilot will produce wrong-direction or "
                "weak-authority pitch/yaw commands.\n"
                "      Set control.thrust_loc (and tvc.parm) to the actual "
                "rear-of-vehicle thrust attachment\n"
                "      point in meters from nose.  For the launcher (xcg_0=8.5m, "
                "18m vehicle) use 18.0.\n",
                ctrl->thrust_loc, max_xcg, max_stage,
                ctrl->thrust_loc - max_xcg);
        }
    }

    load_guidance    (cfg["guidance"],     guid);
    load_ins         (cfg["ins"],          ins);
    load_gps         (cfg["gps"],          gps);
    load_startrack   (cfg["startrack"],    star);
    load_intercept   (cfg["intercept"],    icpt);
    load_launch      (cfg["launch"],       newt, kin, eul);
    load_forces      (cfg["forces"],       forc);
    return 0;
}

double config_sim_dt(const Value& cfg, double def) {
    return cfg["sim"]["dt"].asNumber(def);
}
double config_sim_tmax(const Value& cfg, double def) {
    return cfg["sim"]["t_max"].asNumber(def);
}

} // namespace rocket6dof
