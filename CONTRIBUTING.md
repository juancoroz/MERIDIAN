# Contributing to MERIDIAN

Thanks for your interest. This document explains what kinds of contributions
fit the project, how to set up a development environment, and what the bar
is for getting a pull request merged.

If you're unsure whether a change is wanted, **open an issue first** to
discuss. That saves everyone time.

## Scope

MERIDIAN is a 6DOF rocket flight simulator built on an OSK-style kernel.
The project welcomes:

- **Bug fixes**, especially with a test case that reproduces the bug.
- **Documentation improvements** — README, BUILD_LINUX.md, in-code API
  documentation comments.
- **New regression tests** that cover code paths that aren't currently
  exercised.
- **Performance improvements** that don't change numerical behavior of the
  reference baseline (see "Regression baseline" below).
- **New mission configs** (`mission_*.json`) that exercise interesting
  trajectories or vehicle classes.
- **Additional sensitivity-analysis or sweep capability** (e.g. new sample
  designs in `sobol_runner.cpp` or `mc_runner.cpp`).

Things that are out of scope:

- **Windows-native ports.** The project is Linux-only; Windows users should
  use WSL. The parallel sweep runner uses POSIX `fork()`/`waitpid()`
  directly and we're not maintaining a `#ifdef _WIN32` path. WSL works
  cleanly and is documented in `BUILD_LINUX.md`.
- **Verbatim Zipfel/CADAC++ code.** The physics and subsystem models in
  this project follow Peter H. Zipfel's *Modeling and Simulation of
  Aerospace Vehicle Dynamics* (3rd ed., AIAA) modeling approach, but
  MERIDIAN's implementation is an independent reimplementation. Pull
  requests containing source code copied verbatim from CADAC++ or
  Zipfel's accompanying code will be declined for licensing reasons.
- **Proprietary aerodynamic data.** Only publicly-releasable aero decks
  should be added to the repo.

For anything substantial (new subsystem, new build target, public API
change, refactor across many files), please open an issue or a Discussion
to talk through the approach before writing the code.

## Quick start

1. Fork the repository on GitHub and clone your fork.
2. Set up dependencies per `rocket_6dof/BUILD_LINUX.md`.
3. Build and verify the regression baseline (see next section) **before**
   making changes. If the baseline doesn't pass on `main` for you, the
   problem is with your environment — fix that first.
4. Create a topic branch: `git checkout -b fix/short-description`.
5. Make your change. Run the regression again. Commit.
6. Push to your fork and open a pull request against `main`.

## Regression baseline

Every pull request must keep the existing regression intact:

- `make all` builds every CLI target without errors or new warnings.
- `make mission_gui` builds without errors (you don't need to *run* it
  unless your change touches GUI code, but it must compile).
- All 15 regression test binaries (`env_sweep`, `dyn_test`, `cad_test`,
  `prop_test`, `aero_test`, `intercept_test`, `control_test`,
  `guidance_test`, `ins_test`, `gps_test`, `startrack_test`, `rcs_test`,
  `json_test`, `distributions_test`, `integration_test`) report
  `Total failures: 0`.
- `./mission mission_launcher.json` produces:
  - `Final altitude:   247089.7 m  (247.09 km)`
  - `Final dvbi:       9129.46 m/s`

If your change deliberately alters numerical behavior, the baseline values
must change in lockstep — update them in `README.md`, `BUILD_LINUX.md`, and
the CI workflow in the same PR, and explain the physics reason in the PR
description.

CI runs this same regression on every PR. A green check is required for
merge.

## Code style

The project is plain C++14 — no exceptions, no RTTI, no Boost. Match the
existing style:

- **Indentation:** 4 spaces, no tabs.
- **Naming:** `snake_case` for functions, variables, and file names;
  `PascalCase` for classes and structs; `UPPER_CASE` for constants and
  preprocessor macros.
- **Headers:** include guards use `ROCKET6DOF_FILENAME_H` (or
  `OSK_FILENAME_H` for kernel files).
- **Comments:** describe **what the code is doing**, not how it was
  developed. No TODOs, no `(NEW)` / `(v1)` version markers, no "we'll
  refactor this later" notes. API documentation and algorithm-rationale
  comments (with citations to Zipfel section or equation numbers, etc.)
  are encouraged.
- **No `using namespace` in headers.** In `.cpp` files, prefer explicit
  `std::` qualification.
- **One class per `.h`/`.cpp` pair** where reasonable.

The OSK kernel under `osk/` is comment-free by design — keep it that way.
The OSK examples under `examples/` are pedagogical and byte-matched to
the OSK reference manual; do not modify them without coordinating in an
issue first.

## Tests

If you add a behavior or fix a bug, add a test that would have failed
before your change. The 15 existing `*_test.cpp` files show the pattern:
each test is a standalone binary that exercises a block or subsystem,
prints `Total failures: N` to stdout, and exits with `N`'s sign as the
process exit code.

For changes to subsystems that already have a test file, add a new test
function to that file. For genuinely new subsystems, add a new
`<name>_test.cpp` and wire it into the `Makefile`.

## Commit messages

- First line: imperative mood, 72 characters or fewer
  (e.g. `Fix sign error in yaw control derivation`).
- Blank line.
- Optional body explaining *why* the change is needed (the code shows
  *what* changed). Wrap at 72 characters.
- Reference issues by number where applicable: `Fixes #42`.

## Pull request review

- Keep PRs focused. One logical change per PR makes review tractable.
- Update relevant documentation in the same PR.
- Respond to review comments by pushing additional commits to the same
  branch — don't squash until the PR is approved.
- Be patient. Maintainer time is limited.

## Reporting bugs and asking questions

- **Bugs:** open an issue using the Bug Report template. Include the
  regression-baseline status — most "broken on my machine" reports turn
  out to be environment issues that the baseline check would have caught.
- **Feature ideas:** open an issue using the Feature Request template.
- **Design questions / general discussion:** use GitHub Discussions
  rather than issues.

## License

By contributing, you agree that your contributions will be licensed under
the same MIT License that covers the project (see `LICENSE`). You are
representing that you have the right to license the code you submit.
