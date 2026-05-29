//  mission_runner.cpp -- End-to-end ascent mission, library form
//
//  This file contains run_mission(config_path), the function that
//  composes every block (Env, Prop, Aero, TVC, RCS, Forces, Newton,
//  Euler, Kin, INS, GPS, Startrack, Guidance, Control, Intercept) into
//  one simulation, runs it, and writes mission_log.csv.
//
//  It is the implementation that both the CLI driver (mission.cpp)
//  and the GUI driver (mission_gui.cpp) link against.  The function
//  has external linkage via `extern "C"` so the symbol name is stable
//  across translation units.

#include "../osk/osk.h"
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
#include "mission_config.h"
#include "mission_progress.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <memory>

using namespace rocket6dof;

namespace {

// LOGGER_MARKER

class Logger : public osk::Block {
public:
    Newton*       newt = nullptr;
    Euler*        eul  = nullptr;
    Kinematics*   kin  = nullptr;
    Propulsion*   prop = nullptr;
    Environment*  env  = nullptr;
    Aerodynamics* aero = nullptr;
    Control*      ctrl = nullptr;
    Guidance*     guid = nullptr;
    INS*          ins  = nullptr;
    GPS*          gps  = nullptr;
    Startrack*    star = nullptr;
    double sample_period = 0.1;
    std::FILE* fp = nullptr;
    const char* filename = "mission_log.csv";
    int rows_written = 0;

    // LOGGER_BODY_MARKER

    void getsFromAll(Newton* n, Euler* e, Kinematics* k, Propulsion* p,
                     Environment* en, Aerodynamics* a, Control* c,
                     Guidance* g, INS* i, GPS* gp, Startrack* s) {
        newt = n; eul = e; kin = k; prop = p; env = en; aero = a;
        ctrl = c; guid = g; ins = i; gps = gp; star = s;
    }

    void init() override {
        if (initCount != 0) return;
        fp = std::fopen(filename, "w");
        if (!fp) return;
        std::fprintf(fp, "t,alt,dvbe,dvbi,fpa,heading,lon,lat,");
        std::fprintf(fp, "mach,q_dyn,thrust,vmass,alpha,beta,");
        std::fprintf(fp, "phi,theta,psi,p_deg,q_deg,r_deg,");
        std::fprintf(fp, "delecx,delrcx,ancomx,alcomx,");
        std::fprintf(fp, "ins_pos_err,ins_vel_err,ins_att_err,");
        std::fprintf(fp, "gps_updates,star_updates\n");
    }

    void update() override {}
    // LOGGER_RPT_MARKER

    void rpt() override {
        if (!fp) return;
        if (!osk::State::sample(sample_period)) return;
        if (!newt || !env || !prop || !kin || !eul) return;
        // Skip the very first log call at t=0: kin's Euler-angle outputs
        // are stale from init() (they're not computed until update()
        // runs, but our init order makes the first sample fire before
        // the first update completes).
        if (osk::State::t < 1.0e-9) return;

        // Publish current sim time for the GUI's progress bar.  Cheap
        // store on every log call -- no synchronization needed because
        // the GUI side reads with relaxed semantics.
        rocket6dof::progress::current_t.store(osk::State::t);

        double mach = env->vmach;
        double q    = env->pdynmc;
        double thr  = prop->thrust;
        double vm   = prop->vmass;
        double a    = kin->alphax;
        double bb   = kin->betax;
        double pd   = eul->ppx;
        double qd   = eul->qqx;
        double rd   = eul->rrx;
        double dele = ctrl ? ctrl->delecx : 0.0;
        double delr = ctrl ? ctrl->delrcx : 0.0;
        double anc  = ctrl ? ctrl->ancomx : 0.0;
        double alc  = ctrl ? ctrl->alcomx : 0.0;
        double pe   = ins ? ins->ins_pos_err : 0.0;
        double ve   = ins ? ins->ins_vel_err : 0.0;
        double ae   = ins ? ins->ins_att_err : 0.0;
        int    gu   = ins ? ins->gps_update_count : 0;
        int    su   = ins ? ins->startrack_update_count : 0;

        std::fprintf(fp,
            "%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.6f,%.6f,"
            "%.4f,%.3f,%.3f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.5f,%.5f,"
            "%.4f,%.4f,%.6e,%d,%d\n",
            osk::State::t,
            newt->alt, newt->dvbe, newt->dvbi, newt->thtvdx, newt->psivdx,
            newt->lonx, newt->latx,
            mach, q, thr, vm, a, bb,
            kin->phibdx, kin->thtbdx, kin->psibdx, pd, qd, rd,
            dele, delr, anc, alc,
            pe, ve, ae, gu, su);
        rows_written++;
    }

    ~Logger() { if (fp) std::fclose(fp); }
};

// CONFIG_MARKER

// Forward declare run_mission so main can call it.
int run_mission(const char* config_path);

// RUN_MARKER

int run_mission(const char* config_path) {
    std::printf("================================================================\n");
    std::printf("    rocket_6dof  --  end-to-end mission demonstration\n");
    std::printf("    config: %s\n", config_path);
    std::printf("    build:  MAX_STAGES=%d\n", Propulsion::MAX_STAGES);
    std::printf("================================================================\n\n");

    Environment*  env  = new Environment();
    Propulsion*   prop = new Propulsion();
    Aerodynamics* aero = new Aerodynamics();
    TVC*          tvc  = new TVC();
    RCS*          rcs  = new RCS();
    Control*      ctrl = new Control();
    Guidance*     guid = new Guidance();
    Forces*       forc = new Forces();
    Euler*        eul  = new Euler();
    Kinematics*   kin  = new Kinematics();
    Newton*       newt = new Newton();
    INS*          ins  = new INS();
    GPS*          gps  = new GPS();
    Startrack*    star = new Startrack();
    Intercept*    icpt = new Intercept();
    Logger*       log_ = new Logger();

    // ---- Load config and apply to all blocks ----
    json::Value cfg;
    try {
        cfg = json::parse_file(config_path);
    } catch (const json::ParseError& e) {
        std::printf("ERROR parsing %s at line %d col %d: %s\n",
                    config_path, e.line, e.col, e.what());
        // Continue with defaults -- block constructors set reasonable
        // values, but the result will not match the expected mission.
        std::printf("Continuing with block defaults; results may differ.\n\n");
    }
    apply_config(cfg, prop, aero, tvc, rcs, ctrl, guid,
                 ins, gps, star, icpt, newt, kin, eul, forc);

    // Optional logger config
    {
        const json::Value& lc = cfg["logger"];
        if (lc.isObject()) {
            log_->sample_period = lc["sample_period"].asNumber(log_->sample_period);
            // Note: filename is a const char*, we hold a static buffer
            static std::string logger_filename;
            if (lc["filename"].isString()) {
                logger_filename = lc["filename"].asString();
                log_->filename = logger_filename.c_str();
            }
        } else {
            log_->sample_period = 0.1;
            log_->filename = "mission_log.csv";
        }
    }

    // SETUP_MARKER

    // ---- All block parameters loaded from the JSON config above ----
    // To customize the mission, edit mission.json (or any config path
    // passed on the command line).  See mission_config.h for the
    // schema.  Block constructors set sensible defaults for any key
    // absent from the config.

    // WIRING_MARKER

    env ->getsFrom(newt);
    aero->getsFrom(env, kin, prop, tvc);
    tvc ->getsFrom(prop, ctrl);
    rcs ->getsFrom(prop, newt, kin, ins, guid);
    forc->getsFrom(prop, aero, env, tvc, rcs);
    newt->getsFrom(env, kin, forc, prop);
    kin ->getsFrom(env, newt, eul);
    eul ->getsFrom(forc, prop, kin);
    gps ->getsFrom(newt, eul, ins);
    star->getsFrom(newt, kin, ins);
    ins ->getsFrom(newt, eul, kin, gps, star);
    ctrl->getsFrom(env, newt, eul, kin, aero, prop, ins);
    guid->getsFrom(ctrl, newt, ins, prop);
    icpt->getsFrom(env, newt, eul);
    log_->getsFromAll(newt, eul, kin, prop, env, aero, ctrl, guid, ins, gps, star);

    // STAGES_MARKER

    std::vector<osk::Block*> stage0 = {
        kin, env, prop, aero, tvc, rcs, forc,
        newt, eul,
        gps, star, ins,
        guid, ctrl,
        icpt, log_
    };
    std::vector<std::vector<osk::Block*>> stages = { stage0 };

    // EXECUTE_MARKER

    // Compute burn time for the summary banner
    double burn_time = (prop->fuel_flow_rate > 0.0)
                     ? prop->fmass0 / prop->fuel_flow_rate : 0.0;
    std::printf("Vehicle:    %.0f kg total, %.0f kg propellant, Isp %.0f s, %.1f s burn\n",
                prop->vmass0, prop->fmass0, prop->spi, burn_time);
    std::printf("Launch:     Lat %.1f%s, Lon %.1f%s, FPA %.1f deg\n",
                std::fabs(newt->latx0), newt->latx0 >= 0 ? " N" : " S",
                std::fabs(newt->lonx0), newt->lonx0 >= 0 ? " E" : " W",
                newt->thtvdx0);
    if (guid->mguide == 2) {
        std::printf("Guidance:   attitude program, theta %.1f deg at t=0 -> %.1f deg at t=%.1fs\n",
                    guid->theta_com_start, guid->theta_com_end, guid->t_att_end);
    } else if (guid->mguide == 1) {
        std::printf("Guidance:   accel program, pitch-over at t=%.1fs..%.1fs (%.2f g)\n",
                    guid->t_pitch_start, guid->t_pitch_end, guid->ancomx_program);
    } else {
        std::printf("Guidance:   mguide=%d\n", guid->mguide);
    }
    std::printf("Sensors:    INS mins=%d; GPS mgps=%d @%.1fHz; StarTrack mstar=%d (>%.0fkm)\n",
                ins->mins, gps->mgps,
                gps->gps_step > 0 ? 1.0/gps->gps_step : 0.0,
                star->mstar, star->startrack_alt * 1e-3);
    std::printf("Termination: alt>%.0fkm OR t>%.0fs\n\n",
                icpt->alt_max * 1e-3, icpt->t_max);

    double dt   = config_sim_dt  (cfg, 0.01);
    double tmax = config_sim_tmax(cfg, icpt->t_max + 1.0);
    double dts[] = { dt };
    osk::Sim sim(dts, tmax, stages);

    // Publish progress for the GUI (no-op for the CLI -- the CLI just
    // ignores these atomics).  current_t is reset to 0 here; the
    // Logger will start writing real values as the integration loop
    // ticks.  total_t is the expected end time, used by the GUI for
    // the progress bar denominator.
    rocket6dof::progress::current_t.store(0.0);
    rocket6dof::progress::total_t.store(tmax);
    rocket6dof::progress::is_running.store(true);

    sim.run();

    // Mission done; tell the GUI it can switch from progress bar to
    // post-run plots.  Set is_running false LAST so the GUI sees a
    // consistent (current_t, total_t, is_running=false) on the next
    // poll.
    rocket6dof::progress::current_t.store(tmax);
    rocket6dof::progress::is_running.store(false);

    // ---- Summary ----
    std::printf("\n================================================================\n");
    std::printf("MISSION COMPLETE\n");
    std::printf("================================================================\n");
    std::printf("Termination time: t = %.3f s\n", osk::State::t);
    std::printf("Final altitude:   %.1f m  (%.2f km)\n", newt->alt, newt->alt * 1e-3);
    std::printf("Final dvbi:       %.2f m/s\n", newt->dvbi);
    std::printf("Final dvbe:       %.2f m/s\n", newt->dvbe);
    std::printf("Final FPA:        %.2f deg\n", newt->thtvdx);
    std::printf("Final mass:       %.2f kg  (%.1f kg propellant left, %.1f kg burned)\n",
                prop->vmass, prop->fmassr, prop->fmasse);
    std::printf("Body attitude:    psi=%.2f  theta=%.2f  phi=%.2f deg\n",
                kin->psibdx, kin->thtbdx, kin->phibdx);
    std::printf("INS errors:       pos=%.2f m  vel=%.4f m/s  att=%.2e rad\n",
                ins->ins_pos_err, ins->ins_vel_err, ins->ins_att_err);
    std::printf("Sensor updates:   GPS=%d  StarTrack=%d\n",
                ins->gps_update_count, ins->startrack_update_count);
    std::printf("Log file:         %s (%d rows)\n", log_->filename, log_->rows_written);
    std::printf("================================================================\n");

    delete log_; delete icpt; delete star; delete gps; delete ins;
    delete newt; delete kin; delete eul; delete forc; delete guid;
    delete ctrl; delete rcs; delete tvc; delete aero; delete prop; delete env;
    return 0;
}

} // anon

// External-linkage entry point.  Both mission.cpp (CLI) and
// mission_gui.cpp (GUI) call this.  `extern "C"` is just for symbol-
// name stability; the function takes/returns plain C types so it's
// callable from anywhere without C++ name mangling concerns.
extern "C" int run_mission_ext(const char* config_path) {
    return run_mission(config_path);
}
