# Mock tests

Two layers, sharing one driver ([driver.cpp](driver.cpp)):

1. **Automated golden tests (ctest)** — the driver is built as
   `preciceMockedTestDriver` (linked directly against `libpreciceMocked`) and
   each scenario's `T> ` output is compared against the checked-in files in
   `expected/`. The mock is deterministic, so any diff is a behavior change.

   ```bash
   ctest -L mock            # run from the build directory
   ```

   After an *intended* behavior change, regenerate the affected golden files
   and review the diff:

   ```bash
   UPDATE_EXPECTED=1 ./check.sh <build>/preciceMockedTestDriver "$PWD" \
       explicit-one SolverOne explicit lifecycle
   ```

   The scenario list lives in [tests.cmake](tests.cmake).

2. **Differential harness (manual)** — `run.sh` compares the mock against the
   real `libprecice.so`: the same driver runs each scenario against the real
   library (as a coupled two-participant pair) and against the mock (via
   `LD_PRELOAD`, standalone), and the logged `T> ` lines are diffed. Use this
   to (re)validate the golden files whenever the mock is meant to track a
   change in real preCICE behavior.

## Usage

```bash
# build libprecice and libpreciceMocked first, then:
./run.sh            # run everything and print the diff summary
./run.sh lifecycle  # only the lifecycle scenarios
./run.sh errors     # only the error-probe scenarios
./run.sh diff       # re-compare the recorded runs in out/
```

Set `PRECICE_BUILD_DIR` if your build directory is not `<repo>/build`.

## Scenarios

- `explicit`, `explicit-sub`, `implicit`, `implicit-sub`, `maxtime`, `firstpart`:
  full coupling lifecycles (steps, windows, checkpoints, subcycling, max-time
  truncation, solver-prescribed time-window size).
- `err-*`, `emptymesh`: error probes comparing precondition checks and messages.

Lines starting with `T> V ` carry data values and are excluded from the diff:
the mock does not couple or interpolate data, so values differ by design.

## Expected differences

The diff summary is expected to be `OK` for all lifecycle scenarios except
`firstpart/two`: the second participant's `getMaxTimeStepSize()` is prescribed
by the partner in real preCICE, which the partner-less mock cannot know (window
completion and step counts still match).

Some error scenarios stay `DIFF` on wording only — both libraries throw for the
same misuse, but the mock's messages are deliberately more actionable (e.g.
naming the exact config tag to add). The real library segfaults (Release) or
asserts (Debug) on `advance()` before `initialize()` / after `finalize()`,
where the mock throws a proper error.
