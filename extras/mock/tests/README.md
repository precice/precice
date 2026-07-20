# Mock tests

Test infrastructure for the mocked preCICE library (`extras/mock`), in two
layers:

## `api/` — driver-based API tests

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

## `tutorials/` — tutorial and course integration tests

Runs real adapters (OpenFOAM, python bindings, FMI, nutils, ...) against the
mock.

- **Smoke test** — [tutorials/smoke.sh](tutorials/smoke.sh) runs a curated set
  of short tutorials plus the four course solution sets and compares the
  normalized output against the committed baselines in `tutorials/baselines/`
  (statuses *and* solver output; timestamps, paths, PIDs etc. are filtered).
  Course tasks that fail on purpose (course content) are simply part of the
  baseline. Runtime: ~4.5 minutes; courses run as parallel background jobs.

  ```bash
  tutorials/smoke.sh [--tutorials-only|--courses-only]
  ```

- **Baseline update** — [tutorials/smoke_update.sh](tutorials/smoke_update.sh)
  regenerates the baselines. It first validates every tutorial group as a real
  coupled pair (and the courses with concurrent participants) against the real
  `libprecice`, recording the outcomes in `tutorials/baselines/real-status.txt`,
  then rewrites the mock baselines. Run it whenever preCICE, the tutorials, or
  the course material changes intentionally; use `--skip-real` when only
  normalization changed.

- **Full sweeps (manual)** —
  [tutorials/run_all_tutorials.sh](tutorials/run_all_tutorials.sh) runs *all*
  tutorials against the mock (`--skip-missing`, `--skip-long`, `--timeout N`),
  and [tutorials/run_course_solutions.sh](tutorials/run_course_solutions.sh)
  runs the course solutions with either the mock (`--mock`, sequential) or real
  preCICE (concurrent). [tutorials/gen_mock_configs.py](tutorials/gen_mock_configs.py)
  generates the per-tutorial `precice-mock-config.xml` files (bounded-random
  read data with hand-tuned overrides for physics-sensitive cases).
