//  sim_runner.cpp  --  Single-run lifecycle, shared by MC and Sobol
//
//  Extracted from mc_runner.cpp so the Sobol runner can reuse the
//  same build/run/teardown logic without duplicating it.

#include "sim_runner.h"
#include "mission_config.h"

#include "environment.h"
#include "newton.h"
#include "kinematics.h"
#include "euler.h"
#include "forces.h"
#include "propulsion.h"
#include "aerodynamics.h"
#include "tvc.h"
#include "control.h"
#include "guidance.h"
#include "ins.h"
#include "gps.h"
#include "startrack.h"
#include "rcs.h"
#include "intercept.h"
#include "peak_tracker.h"
#include "../osk/osk.h"

#include <cmath>
#include <cstdio>
#include <set>
#include <fcntl.h>
#include <unistd.h>
#include <vector>

namespace rocket6dof {

namespace {

struct Vehicle {
    Environment*  env;
    Propulsion*   prop;
    Aerodynamics* aero;
    TVC*          tvc;
    RCS*          rcs;
    Control*      ctrl;
    Guidance*     guid;
    Forces*       forc;
    Euler*        eul;
    Kinematics*   kin;
    Newton*       newt;
    INS*          ins;
    GPS*          gps;
    Startrack*    star;
    Intercept*    icpt;
    PeakTracker*  pk;
};

double read_output(const Vehicle& v, const std::string& name) {
    if (name == "alt")            return v.newt->alt;
    if (name == "dvbi")           return v.newt->dvbi;
    if (name == "dvbe")           return v.newt->dvbe;
    if (name == "fpa")            return v.newt->thtvdx;
    if (name == "heading")        return v.newt->psivdx;
    if (name == "lon")            return v.newt->lonx;
    if (name == "lat")            return v.newt->latx;
    if (name == "vmass")          return v.prop->vmass;
    if (name == "fmassr")         return v.prop->fmassr;
    if (name == "fmasse")         return v.prop->fmasse;
    if (name == "phi")            return v.kin->phibdx;
    if (name == "theta")          return v.kin->thtbdx;
    if (name == "psi")            return v.kin->psibdx;
    if (name == "alpha")          return v.kin->alphax;
    if (name == "beta")           return v.kin->betax;
    if (name == "p_deg")          return v.eul->ppx;
    if (name == "q_deg")          return v.eul->qqx;
    if (name == "r_deg")          return v.eul->rrx;
    if (name == "delecx")         return v.ctrl->delecx;
    if (name == "delrcx")         return v.ctrl->delrcx;
    if (name == "ins_pos_err")    return v.ins->ins_pos_err;
    if (name == "ins_vel_err")    return v.ins->ins_vel_err;
    if (name == "ins_att_err")    return v.ins->ins_att_err;
    if (name == "gps_updates")    return static_cast<double>(v.ins->gps_update_count);
    if (name == "startrack_updates") return static_cast<double>(v.ins->startrack_update_count);
    if (name == "t_end")          return osk::State::t;

    // ---- Intercept termination code ----
    // trcond: 0=still running, 1=ground, 2=apogee, 3=alt_max,
    // 4=aero_decay, 5=time_max, 6=tumble.  Useful for MC studies
    // to count which termination conditions fired.
    if (v.icpt) {
        if (name == "trcond")     return static_cast<double>(v.icpt->trcond);
        if (name == "t_terminate") return v.icpt->t_terminate;
    }

    // ---- Per-trajectory peaks (from PeakTracker) ----
    // These accumulate during sim.run() rather than read end-of-run
    // state.  Useful for load-envelope MC studies.
    if (v.pk) {
        if (name == "max_q_dyn")       return v.pk->max_q_dyn;
        if (name == "max_alpha")       return v.pk->max_alpha;
        if (name == "max_q_alpha")     return v.pk->max_q_alpha;
        if (name == "max_mach")        return v.pk->max_mach;
        if (name == "alpha_at_max_q")  return v.pk->alpha_at_max_q;
        if (name == "t_max_q_dyn")     return v.pk->t_max_q_dyn;
        if (name == "t_max_alpha")     return v.pk->t_max_alpha;
        if (name == "t_max_q_alpha")   return v.pk->t_max_q_alpha;
        if (name == "t_max_mach")      return v.pk->t_max_mach;
    }

    // ---- Geographic vs inertial heading ----
    // psi (kin->psibdx) is body psi in INERTIAL frame.  For a vehicle
    // that has traveled significant range, the geographic-frame heading
    // diverges from inertial due to great-circle convergence and Earth
    // rotation.  psi_geographic reports the velocity-frame heading in
    // local geographic coords (psivdx in Newton).
    if (name == "psi_inertial")    return v.kin->psibdx;
    if (name == "psi_geographic")  return v.newt->psivdx;
    if (name == "heading_drift")   return v.kin->psibdx - v.newt->psivdx;

    // Synthetic circular outputs: sin/cos of an angle (deg) for use
    // when running Sobol on a wrapped quantity (psi spans ±180).
    // Variance of the wrapped value is inflated by the wraparound;
    // sin/cos give well-defined Euclidean coordinates.
    const double D2R = 0.017453292519943295;  // pi/180
    if (name == "sin_psi")     return std::sin(v.kin->psibdx * D2R);
    if (name == "cos_psi")     return std::cos(v.kin->psibdx * D2R);
    if (name == "sin_theta")   return std::sin(v.kin->thtbdx * D2R);
    if (name == "cos_theta")   return std::cos(v.kin->thtbdx * D2R);
    if (name == "sin_phi")     return std::sin(v.kin->phibdx * D2R);
    if (name == "cos_phi")     return std::cos(v.kin->phibdx * D2R);
    if (name == "sin_heading") return std::sin(v.newt->psivdx * D2R);
    if (name == "cos_heading") return std::cos(v.newt->psivdx * D2R);

    // ---- Final-state environment quantities ----
    if (name == "mach")        return v.env->vmach;
    if (name == "q_dyn")       return v.env->pdynmc;
    if (name == "rho")         return v.env->rho;
    if (name == "press")       return v.env->press;
    if (name == "tempk")       return v.env->tempk;

    // ---- Propulsion final-state quantities ----
    // Useful for diagnosing burn duration, residual propellant,
    // or per-stage performance dispersion.
    if (v.prop) {
        if (name == "thrust")        return v.prop->thrust;
        if (name == "xcg")           return v.prop->xcg;
        if (name == "current_stage") return static_cast<double>(v.prop->current_stage);
        if (name == "mprop")         return static_cast<double>(v.prop->mprop);
    }

    // ---- Forces / aero loads (final-state) ----
    // Useful when probing peak vs final loads.  For load envelopes
    // use the PeakTracker's max_q_dyn / max_q_alpha; these are
    // instantaneous snapshots at termination.
    if (v.forc) {
        if (name == "FAPB_x")    return v.forc->FAPB[0];
        if (name == "FAPB_y")    return v.forc->FAPB[1];
        if (name == "FAPB_z")    return v.forc->FAPB[2];
    }
    if (v.aero && v.env) {
        // Body-frame axial / lateral / normal aero force in N.
        // refa = ref area, pdynmc = q.  cx/cy/cz are non-dim.
        // Sign convention: cx>0 = drag (rearward).
        double qS = v.env->pdynmc * v.aero->refa;
        if (name == "F_axial")   return qS * v.aero->cx;
        if (name == "F_lateral") return qS * v.aero->cy;
        if (name == "F_normal")  return qS * v.aero->cz;
        if (name == "cx")        return v.aero->cx;
        if (name == "cy")        return v.aero->cy;
        if (name == "cz")        return v.aero->cz;
    }

    // ---- Newton body acceleration (in g's) ----
    // These were *only* available via PeakTracker as anx_at_max_q
    // (a peak measure).  Final-state instantaneous values are
    // useful for diagnosing landing/re-entry envelopes.
    if (v.newt) {
        if (name == "anx")       return v.newt->anx;
        if (name == "ayx")       return v.newt->ayx;
        if (name == "axx")       return v.newt->axx;
        if (name == "vbi_x")     return v.newt->VBII[0];
        if (name == "vbi_y")     return v.newt->VBII[1];
        if (name == "vbi_z")     return v.newt->VBII[2];
    }

    // ---- Control commanded attitudes (diagnostic) ----
    // Useful for diagnosing autopilot tracking error: theta_act
    // minus theta_com_eff at termination.  Only meaningful for
    // maut=60 mode.
    if (v.ctrl) {
        if (name == "theta_com")      return v.ctrl->theta_com_inertial;
        if (name == "psi_com")        return v.ctrl->psi_com_inertial;
        if (name == "ancomx_actual")  return v.ctrl->ancomx_actual;
        if (name == "alcomx_actual")  return v.ctrl->alcomx_actual;
    }

    // ---- Unknown output: warn once per name and return NaN ----
    // Warn on stderr so misconfigured MC/Sobol configs (e.g., output
    // names not in read_output) are visible rather than producing
    // all-NaN results.  The warning fires once per unique unknown
    // name to avoid spamming stderr in a large MC.
    {
        static std::set<std::string> warned;
        if (warned.find(name) == warned.end()) {
            warned.insert(name);
            std::fprintf(stderr,
                "[read_output] WARNING: unknown output name '%s' "
                "-- returning NaN.  Check config 'outputs' list.\n",
                name.c_str());
        }
    }
    return std::nan("");
}

void destroy_vehicle(Vehicle& v) {
    delete v.pk;
    delete v.icpt; delete v.star; delete v.gps; delete v.ins;
    delete v.newt; delete v.kin; delete v.eul; delete v.forc;
    delete v.guid; delete v.ctrl; delete v.rcs; delete v.tvc;
    delete v.aero; delete v.prop; delete v.env;
}

// RAII helper: dup/dup2 stdout/stderr to /dev/null on construction,
// restore on destruction.  Necessary because each block's rpt()
// prints during the sim; otherwise the MC/Sobol output is unreadable.
class StdoutSilencer {
public:
    StdoutSilencer() : saved_stdout_(-1), saved_stderr_(-1) {
        std::fflush(stdout);
        std::fflush(stderr);
        saved_stdout_ = ::dup(1);
        saved_stderr_ = ::dup(2);
        int devnull = ::open("/dev/null", O_WRONLY);
        if (devnull >= 0) {
            ::dup2(devnull, 1);
            ::dup2(devnull, 2);
            ::close(devnull);
        }
    }
    ~StdoutSilencer() {
        std::fflush(stdout);
        std::fflush(stderr);
        if (saved_stdout_ >= 0) { ::dup2(saved_stdout_, 1); ::close(saved_stdout_); }
        if (saved_stderr_ >= 0) { ::dup2(saved_stderr_, 2); ::close(saved_stderr_); }
    }
private:
    int saved_stdout_;
    int saved_stderr_;
};

} // anon

SingleRunResult run_single_mission(const json::Value& cfg,
                                   const std::vector<std::string>& output_names)
{
    SingleRunResult r;
    r.outputs.assign(output_names.size(), std::nan(""));
    r.t_end    = 0;
    r.nan_seen = false;

    Vehicle v;
    v.env  = new Environment();
    v.prop = new Propulsion();
    v.aero = new Aerodynamics();
    v.tvc  = new TVC();
    v.rcs  = new RCS();
    v.ctrl = new Control();
    v.guid = new Guidance();
    v.forc = new Forces();
    v.eul  = new Euler();
    v.kin  = new Kinematics();
    v.newt = new Newton();
    v.ins  = new INS();
    v.gps  = new GPS();
    v.star = new Startrack();
    v.icpt = new Intercept();
    v.pk   = new PeakTracker();

    apply_config(cfg, v.prop, v.aero, v.tvc, v.rcs, v.ctrl, v.guid,
                 v.ins, v.gps, v.star, v.icpt, v.newt, v.kin, v.eul, v.forc);

    v.env ->getsFrom(v.newt);
    v.aero->getsFrom(v.env, v.kin, v.prop, v.tvc);
    v.tvc ->getsFrom(v.prop, v.ctrl);
    v.rcs ->getsFrom(v.prop, v.newt, v.kin, v.ins, v.guid);
    v.forc->getsFrom(v.prop, v.aero, v.env, v.tvc, v.rcs);
    v.newt->getsFrom(v.env, v.kin, v.forc, v.prop);
    v.kin ->getsFrom(v.env, v.newt, v.eul);
    v.eul ->getsFrom(v.forc, v.prop, v.kin);
    v.gps ->getsFrom(v.newt, v.eul, v.ins);
    v.star->getsFrom(v.newt, v.kin, v.ins);
    v.ins ->getsFrom(v.newt, v.eul, v.kin, v.gps, v.star);
    v.ctrl->getsFrom(v.env, v.newt, v.eul, v.kin, v.aero, v.prop, v.ins);
    v.guid->getsFrom(v.ctrl, v.newt, v.ins, v.prop);
    v.icpt->getsFrom(v.env, v.newt, v.eul);
    v.pk  ->getsFrom(v.env, v.kin);

    std::vector<osk::Block*> stage0 = {
        v.kin, v.env, v.prop, v.aero, v.tvc, v.rcs, v.forc,
        v.newt, v.eul,
        v.gps, v.star, v.ins,
        v.guid, v.ctrl,
        v.icpt,
        v.pk
    };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    double dt   = config_sim_dt  (cfg, 0.01);
    double tmax = config_sim_tmax(cfg, v.icpt->t_max + 1.0);
    double dts[] = { dt };

    {
        StdoutSilencer silence;
        osk::Sim sim(dts, tmax, stages);
        sim.run();
    }

    r.t_end = osk::State::t;
    for (size_t i = 0; i < output_names.size(); ++i) {
        double val = read_output(v, output_names[i]);
        r.outputs[i] = val;
        if (!std::isfinite(val)) r.nan_seen = true;
    }
    destroy_vehicle(v);
    return r;
}

} // namespace rocket6dof
