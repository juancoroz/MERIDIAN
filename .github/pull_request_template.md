<!-- Thanks for contributing.  Please fill in the sections below. -->

## Summary
<!-- One or two sentences describing what this PR changes. -->


## Motivation
<!-- Why is this change needed?  Link any related issue with "Fixes #N"
     or "Refs #N". -->


## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds capability)
- [ ] Refactor (no functional change)
- [ ] Documentation only
- [ ] Build / CI / tooling
- [ ] Breaking change (numerical behavior changed, public API changed,
      or baseline values shift)


## Regression baseline

Required for any PR that touches code under `osk/` or `rocket_6dof/`:

- [ ] `make all` builds clean (no errors, no new warnings)
- [ ] `make mission_gui` builds clean
- [ ] All 15 regression tests pass (`Total failures: 0` for each)
- [ ] `./mission mission_launcher.json` produces `Final altitude:
      247089.7 m` and `Final dvbi: 9129.46 m/s`

If this PR **deliberately changes** the baseline (e.g. a physics fix
that shifts the trajectory), check this box instead and explain in the
"Numerical change" section below:

- [ ] Baseline values are updated in `README.md`, `BUILD_LINUX.md`, and
      `.github/workflows/build.yml`, and the change is justified below.


## Numerical change (only if baseline shifted)
<!-- What changed, by how much, and why is the new value correct? -->


## Tests
<!-- Did you add tests?  Which existing tests cover this code path?
     If the change is hard to test, explain why. -->


## Documentation
<!-- Did you update README.md, BUILD_LINUX.md, in-code API comments,
     or anything else user-facing?  If no docs needed, say so. -->


## Checklist

- [ ] My code matches the project style (see CONTRIBUTING.md)
- [ ] Comments describe what the code does, not how it was developed
      (no TODOs, version markers, refactoring history)
- [ ] If user-facing behavior changed, the relevant docs are updated
- [ ] If the third-party version requirements changed, the CI
      workflow's `IMGUI_VERSION` / `IMPLOT_VERSION` are updated
- [ ] I have read CONTRIBUTING.md and understand the licensing
      requirements (MIT; no Zipfel/CADAC++ verbatim code)
