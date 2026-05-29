---
name: Bug report
about: Report a problem with the simulator
title: ''
labels: bug
assignees: ''
---

## Description
<!-- Brief description of the bug.  What did you do, what happened, what
     did you expect to happen? -->


## Steps to reproduce
<!-- A short, self-contained sequence anyone can follow.  If a particular
     mission config triggers the bug, paste its relevant fields here or
     attach the JSON file. -->

1.
2.
3.


## Expected behavior


## Actual behavior


## Environment

- **OS and version:** <!-- e.g. Ubuntu 24.04, native -->
- **WSL?** <!-- yes / no, and which Windows version if yes -->
- **Compiler and version:** <!-- output of `g++ --version` -->
- **MERIDIAN commit:** <!-- output of `git rev-parse --short HEAD` -->
- **Third-party versions used:** Dear ImGui <!-- e.g. v1.91.6 -->, ImPlot <!-- e.g. v0.16 -->


## Regression baseline status

The single most useful piece of information for triaging "broken on my
machine" reports is whether the regression baseline still passes on your
build.  Please run from `rocket_6dof/`:

```sh
make all
./mission mission_launcher.json
```

- [ ] All 15 regression tests pass (`env_sweep`, `dyn_test`, `cad_test`,
      `prop_test`, `aero_test`, `intercept_test`, `control_test`,
      `guidance_test`, `ins_test`, `gps_test`, `startrack_test`,
      `rcs_test`, `json_test`, `distributions_test`,
      `integration_test` all report `Total failures: 0`).
- [ ] `mission mission_launcher.json` produces `Final altitude:
      247089.7 m` and `Final dvbi: 9129.46 m/s`.

If either box is unchecked, the issue is likely environmental (toolchain,
third-party library version, file corruption) rather than a code bug;
please look at your local setup before filing.


## Additional context
<!-- Log output, screenshots, related issues, anything else. -->
