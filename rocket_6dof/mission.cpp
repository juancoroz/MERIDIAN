//  mission.cpp -- CLI driver for the rocket_6dof simulator
//
//  Thin wrapper around run_mission_ext() (in mission_runner.cpp).
//
//  Usage:
//    ./mission                 -- loads mission.json
//    ./mission config.json     -- loads named config
//
//  See mission_runner.cpp for the actual simulation pipeline.  Both
//  this binary and mission_gui link mission_runner.o.

extern "C" int run_mission_ext(const char* config_path);

int main(int argc, char** argv) {
    const char* config_path = (argc >= 2) ? argv[1] : "mission.json";
    return run_mission_ext(config_path);
}
