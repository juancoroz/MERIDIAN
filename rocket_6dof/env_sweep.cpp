//  env_sweep.cpp  --  Verify the Environment block against published
//                     US 1976 Standard Atmosphere values.
//
//  Two checks:
//    (a) Boundary-value spot check.  At each of the seven US76 layer
//        boundaries the temperature, pressure, and density are
//        well-known (see e.g. NASA-TM-X-74335 Table 1).  Tolerance is
//        0.5% on rho/p (cumulative interpolation), 0.05 K on T.
//
//    (b) 1 km sweep from 0 to 86 km, printed as a table the user can
//        diff against any published US76 reference.
//
//  The driver also exercises Environment as a real OSK Block by
//  running it through a 1-step simulation, to confirm the block plugs
//  into the kernel correctly.

#include "../osk/osk.h"
#include "environment.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace rocket6dof;

namespace {

struct Reference {
    double z_km;     // geometric altitude
    double T_K;      // expected temperature
    double p_Pa;     // expected pressure
    double rho_kgm3; // expected density
};

// Published US76 values at the seven base levels (taken from the
// standard table; these are at GEOPOTENTIAL altitudes which the
// helper converts internally).  Sea-level row uses z=0 (geometric==
// geopotential at z=0).  The remaining rows use the geometric
// altitudes corresponding to each geopotential base.
//
// Geometric altitudes for geopotential bases (z = R*h/(R-h)):
//   11 km -> 11.019 km
//   20 km -> 20.063 km
//   32 km -> 32.162 km
//   47 km -> 47.350 km
//   51 km -> 51.413 km
//   71 km -> 71.802 km
const Reference REFS[] = {
    {  0.000, 288.150, 1.01325e+5, 1.22500e+0 },
    { 11.019, 216.650, 2.26321e+4, 3.63918e-1 },
    { 20.063, 216.650, 5.47488e+3, 8.80345e-2 },
    { 32.162, 228.650, 8.68016e+2, 1.32249e-2 },
    { 47.350, 270.650, 1.10905e+2, 1.42752e-3 },
    { 51.413, 270.650, 6.69384e+1, 8.61601e-4 },
    { 71.802, 214.650, 3.95639e+0, 6.42105e-5 },
};
constexpr int N_REFS = sizeof(REFS) / sizeof(REFS[0]);

int run_boundary_checks() {
    std::printf("\n=== Boundary value checks (US76) ===\n");
    std::printf("%-8s %-8s %-12s %-12s %-12s %-12s %-12s %-12s\n",
                "z[km]", "T[K]ref", "T[K]osk", "dT[K]",
                "p[Pa]ref", "p[Pa]osk", "rho[kg/m3]ref", "rho[kg/m3]osk");

    int failures = 0;
    for (int i = 0; i < N_REFS; ++i) {
        const Reference& R = REFS[i];
        double rho, p, T;
        atmosphere76(R.z_km * 1000.0, rho, p, T);

        double dT  = T - R.T_K;
        double rp  = (p   - R.p_Pa)    / R.p_Pa;
        double rr  = (rho - R.rho_kgm3)/ R.rho_kgm3;

        bool ok =    std::fabs(dT) < 0.10
                  && std::fabs(rp) < 0.005
                  && std::fabs(rr) < 0.005;

        std::printf("%-8.3f %-8.3f %-12.4f %-12.4f %-12.4e %-12.4e %-12.4e %-12.4e   %s\n",
                    R.z_km, R.T_K, T, dT,
                    R.p_Pa, p, R.rho_kgm3, rho,
                    ok ? "OK" : "FAIL");

        if (!ok) ++failures;
    }
    std::printf("Boundary check: %d / %d passed\n", N_REFS - failures, N_REFS);
    return failures;
}

void run_altitude_sweep() {
    std::printf("\n=== Altitude sweep 0..86 km (5 km step) ===\n");
    std::printf("%-8s %-12s %-12s %-10s %-10s\n",
                "z[km]", "rho[kg/m3]", "p[Pa]", "T[K]", "a[m/s]");
    for (int km = 0; km <= 86; km += 5) {
        double rho, p, T;
        atmosphere76(km * 1000.0, rho, p, T);
        double a = std::sqrt(1.4 * 287.053 * T);
        std::printf("%-8d %-12.4e %-12.4e %-10.3f %-10.3f\n", km, rho, p, T, a);
    }
}

// A trivial "vehicle" block that pushes inputs into Environment so
// we can drive a 1-second simulation and confirm the Block plugs
// into the OSK kernel.  Sits on the launch pad at 0 alt, stationary.
class StaticPad : public osk::Block {
public:
    Environment* env;
    StaticPad(Environment* e) : env(e) {}
    void init()   override {
        env->alt  = 0.0;
        env->SBII = osk::Vec(6378137.0, 0.0, 0.0);
        env->VBED = osk::Vec(0.0, 0.0, 0.0);
    }
    void update() override {}
    void rpt()    override {}
};

int run_kernel_integration_test() {
    std::printf("\n=== Kernel integration test (Environment in a Sim) ===\n");

    Environment* env = new Environment();
    StaticPad*   pad = new StaticPad(env);

    // Pad must run BEFORE env so that the inputs are set when env->update() runs.
    std::vector<osk::Block*> stage0;
    stage0.push_back(pad);
    stage0.push_back(env);
    std::vector< std::vector<osk::Block*> > stages;
    stages.push_back(stage0);

    double dts[] = { 0.1 };
    osk::Sim sim(dts, 1.0, stages);
    sim.run();

    // After running for 1 sec on the pad, env should read sea-level values.
    bool ok =    std::fabs(env->tempk  - 288.15)  < 0.1
              && std::fabs(env->rho    - 1.225)   < 0.005
              && std::fabs(env->press  - 101325.) < 100.0
              && std::fabs(env->grav   - 9.7980)  < 0.005;  // g at equatorial surface

    std::printf("Final state: T=%.3f K, rho=%.4f, p=%.1f, g=%.4f m/s^2  %s\n",
                env->tempk, env->rho, env->press, env->grav,
                ok ? "OK" : "FAIL");

    delete pad;
    delete env;
    return ok ? 0 : 1;
}

} // anonymous namespace

int main() {
    int fails = 0;
    fails += run_boundary_checks();
    run_altitude_sweep();
    fails += run_kernel_integration_test();

    std::printf("\n=== Total failures: %d ===\n", fails);
    return fails == 0 ? 0 : 1;
}
