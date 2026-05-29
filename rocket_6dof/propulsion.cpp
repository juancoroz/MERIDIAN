//  propulsion.cpp  --  Multi-stage rocket propulsion + mass properties
//
//  State machine per stage:
//
//    BURNING (phase=0):
//      thrust = spi * mdot * g0
//      fmasse integrates at rate mdot
//      vmass = vmass0_stage[s] - fmasse
//      IBBB, xcg interpolate by fmasse/fmass0_stage[s]
//      transitions to COASTING when fmassr <= 0
//
//    COASTING (phase=1):
//      thrust = 0; mass/MoI/xcg held at burnout values
//      transitions to BURNING (next stage) when t >= t_coast_end
//      at transition: drop dry_mass_dropped[s], set s=s+1, fmasse <- 0
//
//    SPENT (phase=2):
//      thrust = 0; mass/MoI/xcg held at final values
//      no further transitions
//
//  Staging transitions only modify the integrator scalar fmasse at
//  osk::State::stepstart=true, matching the RK4 timing rule used by
//  GPS/INS/Startrack.

#include "propulsion.h"
#include "environment.h"
#include "cad_utils.h"
#include <cmath>
#include <cstdio>

namespace rocket6dof {

namespace {
constexpr double G0_PROP = 9.80675445;
} // anon

Propulsion::Propulsion()
    : env(nullptr)
{
    mprop = 0;
    vmass0 = 1.0;  fmass0 = 0.0;
    spi    = 0.0;  fuel_flow_rate = 0.0;
    xcg_0  = 0.0;  xcg_1  = 0.0;
    moi_roll_0  = 1.0;  moi_roll_1  = 1.0;
    moi_trans_0 = 1.0;  moi_trans_1 = 1.0;

    num_stages = 1;
    for (int s = 0; s < MAX_STAGES; s++) {
        vmass0_stage[s]      = 0.0;
        fmass0_stage[s]      = 0.0;
        spi_stage[s]         = 0.0;
        mdot_stage[s]        = 0.0;
        xcg_0_stage[s]       = 0.0;
        xcg_1_stage[s]       = 0.0;
        moi_roll_0_stage[s]  = 1.0;
        moi_roll_1_stage[s]  = 1.0;
        moi_trans_0_stage[s] = 1.0;
        moi_trans_1_stage[s] = 1.0;
        dry_mass_dropped[s]  = 0.0;
        coast_to_next[s]     = 0.0;
    }

    fmasse = 0.0;
    fmassd = 0.0;
    current_stage = 0;
    phase = 0;
    t_stage_start = 0.0;
    t_coast_end   = 0.0;

    thrust = 0.0;
    vmass  = vmass0;
    fmassr = fmass0;
    xcg    = xcg_0;
    IBBB   = osk::Mat(moi_roll_0, 0.0, 0.0,
                      0.0, moi_trans_0, 0.0,
                      0.0, 0.0, moi_trans_0);

    addIntegrator(fmasse, fmassd);
}

void Propulsion::init() {
    if (initCount == 0) {
        if (num_stages == 1) {
            vmass0_stage[0]      = vmass0;
            fmass0_stage[0]      = fmass0;
            spi_stage[0]         = spi;
            mdot_stage[0]        = fuel_flow_rate;
            xcg_0_stage[0]       = xcg_0;
            xcg_1_stage[0]       = xcg_1;
            moi_roll_0_stage[0]  = moi_roll_0;
            moi_roll_1_stage[0]  = moi_roll_1;
            moi_trans_0_stage[0] = moi_trans_0;
            moi_trans_1_stage[0] = moi_trans_1;
            dry_mass_dropped[0]  = 0.0;
            coast_to_next[0]     = 0.0;
        }
        fmasse = 0.0;
        fmassd = 0.0;
        current_stage = 0;
        phase = 0;
        t_stage_start = 0.0;
        t_coast_end   = 0.0;

        vmass  = vmass0_stage[0];
        fmassr = fmass0_stage[0];
        xcg    = xcg_0_stage[0];
        IBBB   = osk::Mat(moi_roll_0_stage[0],  0.0, 0.0,
                          0.0, moi_trans_0_stage[0], 0.0,
                          0.0, 0.0, moi_trans_0_stage[0]);
        thrust = 0.0;
    }
}

// UPDATE_MARKER

void Propulsion::update() {
    double t = osk::State::t;
    int s = current_stage;

    // ---- mprop == 0: explicitly off (no thrust, no integration) ----
    if (mprop == 0) {
        thrust = 0.0;
        fmassd = 0.0;
        return;
    }

    if (mprop != 3) {
        // Unknown mode; treat as off
        thrust = 0.0;
        fmassd = 0.0;
        return;
    }

    // ---- Phase 2 (spent): all stages exhausted, no more thrust ----
    if (phase == 2) {
        thrust = 0.0;
        fmassd = 0.0;
        mprop  = 0;     // persistent shutoff
        return;
    }

    // ---- Build current-stage MOI bookends ----
    osk::Mat IBBB0(moi_roll_0_stage[s],  0.0, 0.0,
                   0.0, moi_trans_0_stage[s], 0.0,
                   0.0, 0.0, moi_trans_0_stage[s]);
    osk::Mat IBBB1(moi_roll_1_stage[s],  0.0, 0.0,
                   0.0, moi_trans_1_stage[s], 0.0,
                   0.0, 0.0, moi_trans_1_stage[s]);

    // ---- Phase 0: BURNING ----
    if (phase == 0) {
        thrust = spi_stage[s] * mdot_stage[s] * G0_PROP;
        if (spi_stage[s] > 0.0) {
            fmassd = thrust / (spi_stage[s] * G0_PROP);
        } else {
            fmassd = 0.0;
        }
        vmass  = vmass0_stage[s] - fmasse;
        fmassr = fmass0_stage[s] - fmasse;

        double mass_ratio = (fmass0_stage[s] > 0.0)
                          ? (fmasse / fmass0_stage[s]) : 0.0;
        if (mass_ratio < 0.0) mass_ratio = 0.0;
        if (mass_ratio > 1.0) mass_ratio = 1.0;

        IBBB = IBBB0 + (IBBB1 - IBBB0) * mass_ratio;
        xcg  = xcg_0_stage[s] + (xcg_1_stage[s] - xcg_0_stage[s]) * mass_ratio;

        // Burnout transition: stage exhausted
        if (fmassr <= 0.0) {
            // Clamp to exact burnout values
            fmasse = fmass0_stage[s];
            vmass  = vmass0_stage[s] - fmass0_stage[s];
            fmassr = 0.0;
            IBBB   = IBBB1;
            xcg    = xcg_1_stage[s];
            thrust = 0.0;
            fmassd = 0.0;

            // Check if this is the last stage
            if (s + 1 >= num_stages) {
                phase = 2;        // SPENT
                mprop = 0;        // persistent shutoff
            } else {
                phase = 1;        // COAST
                t_coast_end = t + coast_to_next[s];
            }
        }
        return;
    }

    // ---- Phase 1: COASTING between stages ----
    if (phase == 1) {
        thrust = 0.0;
        fmassd = 0.0;
        // Mass/MoI/xcg are held at the burnout values from the previous
        // stage (already set when we entered phase 1).

        // End of coast: drop dry mass of prior stage, switch to next stage,
        // reset fmasse to 0.  This integrator-state modification must only
        // happen at stepstart=true (RK4 rule).
        if (t >= t_coast_end && osk::State::stepstart) {
            int next = s + 1;
            // Apply the dry-mass drop: the vehicle is lighter by
            // dry_mass_dropped[s].  The next stage's vmass0_stage[next]
            // is what the user said the vehicle weighs at start of stage
            // next.  We trust their parameterization here.  If it's
            // inconsistent with (vmass - dry_mass_dropped[s]), the
            // discrepancy is silent (a future enhancement could warn).
            current_stage = next;
            phase = 0;
            t_stage_start = t;
            fmasse = 0.0;     // reset (kernel will pick this up at stage 0)
            fmassd = 0.0;

            // Reset outputs to the new stage's initial values
            vmass  = vmass0_stage[next];
            fmassr = fmass0_stage[next];
            xcg    = xcg_0_stage[next];
            IBBB   = osk::Mat(moi_roll_0_stage[next],  0.0, 0.0,
                              0.0, moi_trans_0_stage[next], 0.0,
                              0.0, 0.0, moi_trans_0_stage[next]);
        }
        return;
    }
}

void Propulsion::rpt() {
    if (osk::State::sample(1.0)) {
        const char* phase_name = "burn";
        if (phase == 1) phase_name = "coast";
        if (phase == 2) phase_name = "spent";
        if (num_stages > 1) {
            std::printf("Prop t=%7.3f  mprop=%d  stage=%d/%d %-5s  "
                        "thrust=%9.1f N  vmass=%8.2f kg  fmassr=%8.2f kg\n",
                        osk::State::t, mprop, current_stage + 1, num_stages,
                        phase_name, thrust, vmass, fmassr);
        } else {
            std::printf("Prop t=%7.3f  mprop=%d  thrust=%9.1f N  vmass=%8.2f kg  "
                        "fmassr=%8.2f kg\n",
                        osk::State::t, mprop, thrust, vmass, fmassr);
        }
    }
}

} // namespace rocket6dof
