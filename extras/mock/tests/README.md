# Mock tests

Test infrastructure for the mocked preCICE library (`extras/mock`).

A synthetic driver ([api/driver.cpp](api/driver.cpp)) exercises the Participant
API scenario by scenario and logs observable behavior as `T> ` lines.

1. **Golden tests (ctest, automated)** — the driver is built as
   `preciceMockedTestDriver` (linked directly against `libpreciceMocked`) and
   each scenario's output is compared with the files in `api/expected/`. The
   mock is deterministic, so any diff is a behavior change.

   ```bash
   ctest -L mock            # run from the build directory
   ```

   After an *intended* behavior change, regenerate the affected golden files
   and review the diff:

   ```bash
   UPDATE_EXPECTED=1 api/check.sh <build>/preciceMockedTestDriver "$PWD/api" \
       explicit-one SolverOne explicit lifecycle
   ```

   The scenario list lives in [api/tests.cmake](api/tests.cmake).

2. **Differential harness (manual)** — [api/run.sh](api/run.sh) compares the
   mock against the real `libprecice.so`: the same driver runs each scenario
   against the real library (as a coupled two-participant pair) and against the
   mock (`LD_PRELOAD`, standalone), and the `T> ` lines are diffed. Use this to
   (re)validate the golden files whenever the mock is meant to track a change
   in real preCICE behavior.

## Integration tests against real adapters

The mock has also been run against the preCICE tutorials and the course
material (OpenFOAM, the Python bindings, FMI, nutils, ...). Those cases live in
repositories of their own and evolve independently of preCICE, so the harness
for them is intentionally not part of this repository. It would be a good fit
for the [system tests](https://github.com/precice/tutorials) once they support
running a tutorial against the mock.
