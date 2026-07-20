Development Scripts
===================

This folder contains helper scripts for running preCICE tutorials and course solutions with both real preCICE coupling and the mock library as well as sample outputs.

Outputs
-------

### mock_all_tutorials.log

Complete run log from `run_all_tutorials.sh` with the mock library. Contains:
- Per-tutorial status (passed ✓, failed ✗, timeout ⏱)
- Error messages and exit codes
- Summary statistics: 22 passed, 67 failed, 3 timeout (out of 92 total)
- Lists of failed/timeout tutorials grouped by category

### mock_course.log

Course solutions run with `run_course_solutions.sh --mock`. Shows:
- Per-task status and outputs from generator/propagator processes
- Mock library version and configuration messages
- Summary: 12 passed, 5 failed, 0 timeout (out of 17 tasks)
- Failed solutions listed (e.g., B/solution/T1, MAP/solution/T3)

### nomock_course.log

Same course run using real preCICE (no mock library). Enables comparison:
- Real preCICE version and build info
- Socket communication and coupling iterations
- Identical pass/fail pattern: 12 passed, 5 failed (same failures as mock)
- Useful for verifying mock library equivalence

Scripts
-------

### run_all_tutorials.sh

Finds and executes all tutorial `run.sh` scripts in a tutorials directory,
optionally using the preCICE mock library.

**Usage:**
```bash
./run_all_tutorials.sh [OPTIONS]
```

**Options:**
- `--mock-lib PATH`: Path to mocked library (default: auto-detected from repo root)
- `--tutorials-dir PATH`: Path to tutorials directory (default: auto-detected)
- `--timeout SECONDS`: Timeout per tutorial (default: 300)
- `--mock-config PATH`: Mock config XML to copy into each tutorial
- `--help`: Show usage

**Features:**
- Automatic repo root detection
- Per-tutorial timeout enforcement
- Mock config distribution to tutorials
- Detailed pass/fail/timeout reporting
- CSV output of results

**Recommended usage:**
```bash
cd /path/to/precice/extras/mock/
./run_all_tutorials.sh &> mock_all_tutorials.log
```

---

### run_course_solutions.sh

Finds and executes course solution scripts, optionally with mock library.
Supports filtering by course number (B, T, CHT, MAP).

**Usage:**
```bash
./run_course_solutions.sh [OPTIONS]
```

**Options:**
- `--mock`: Use mocked preCICE (LD_PRELOAD)
- `--mock-lib PATH`: Path to mocked library
- `--course-number N`: Run specific course (1=B, 2=T, 3=CHT, 4=MAP)
- `--course-dir PATH`: Course directory path (default: auto-detected)
- `--timeout SECONDS`: Timeout per solution (default: 300)
- `--help`: Show usage

**Features:**
- Real vs. mock mode selection
- Course filtering
- Solution discovery patterns: `run_solution.sh`, `*/solutions/run.sh`, `*/solution*/run.sh`
- Timeout enforcement
- Summary statistics

**Recommended usage:**
```bash
./run_course_solutions.sh --mock &> mock_course.log
```

---

Configuration
--------------

Both scripts auto-detect the following from the repo structure:
- Mock library location: `<repo-root>/build/libpreciceMocked.so.3.3.0`
- Tutorials directory: `<repo-root>/examples/tutorials`
- Course directory: `<repo-root>/examples/course`

Override via command-line options if needed.

Notes
-----

- Scripts must be run from or reference a valid preCICE repository with the
  expected directory layout.
- Each script captures logs, pass/fail status, and timing for diagnostics.
- Mock config support in run_all_tutorials.sh allows consistent mock setup
  across all tutorials.
