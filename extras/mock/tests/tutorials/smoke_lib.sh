# Shared definitions for the tutorial/course smoke tests (sourced by smoke.sh
# and smoke_update.sh). Not executable on its own.

SMOKE_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SMOKE_LIB_DIR/../../../.." && pwd)"

MOCK_LIB="${MOCK_LIB:-$REPO_ROOT/build/libpreciceMocked.so.3.3.0}"
TUTORIALS_DIR="${TUTORIALS_DIR:-$REPO_ROOT/examples/tutorials}"
COURSE_DIR="${COURSE_DIR:-$REPO_ROOT/examples/course}"
BASELINE_DIR="$SMOKE_LIB_DIR/baselines"
VENV_DIR="$REPO_ROOT/.venv"

# Curated short, dependency-light tutorials. Each group is one coupled setup:
#   <slug>|<tutorial dir>|<participant case> <participant case> ...
# In mock mode the participants run sequentially; in real mode (update script)
# they run concurrently as an actual coupled pair.
SMOKE_GROUPS=(
    "elastic-tube-1d-cpp|elastic-tube-1d|fluid-cpp solid-cpp"
    "elastic-tube-1d-python|elastic-tube-1d|fluid-python solid-python"
    "oscillator-python|oscillator|mass-left-python mass-right-python"
    "oscillator-fmi|oscillator|mass-left-fmi mass-right-fmi"
    "oscillator-overlap|oscillator-overlap|mass-left-python mass-right-python"
    "quickstart|quickstart|fluid-openfoam solid-cpp"
    "perpendicular-flap|perpendicular-flap|fluid-fake solid-openfoam"
    "partitioned-heat-conduction|partitioned-heat-conduction|dirichlet-nutils neumann-nutils"
)

# Courses by the numbers used by run_course_solutions.sh
SMOKE_COURSES=("1:B" "2:T" "3:CHT" "4:MAP")

smoke_setup_env() {
    # Same environment the tutorial runners use: venv python + pip constraints
    if [ -x "$VENV_DIR/bin/python3" ]; then
        export PATH="$VENV_DIR/bin:$PATH"
        export VIRTUAL_ENV="$VENV_DIR"
    fi
    if [ -f "$SMOKE_LIB_DIR/pip-constraints.txt" ]; then
        export PIP_CONSTRAINT="$SMOKE_LIB_DIR/pip-constraints.txt"
    fi
    # Headless: several course/tutorial scripts end in plt.show(), which blocks
    # until the plot window is closed and would burn the whole timeout.
    export MPLBACKEND=Agg
}

# Strip everything that legitimately varies between identical runs: timestamps,
# durations, hosts/PIDs, absolute paths, pip/venv chatter, build-system noise.
smoke_normalize() {
    # Every line filter tolerates leading whitespace: the course runner indents
    # participant output by four spaces.
    sed -E \
        -e "s|$REPO_ROOT|<REPO>|g" \
        -e 's|/tmp/[A-Za-z0-9_./-]+|<TMP>|g' \
        -e 's|file:///[^ ]*|<LOGURL>|g' \
        -e '/^[[:space:]]*(Started on|Finished on|Duration):/d' \
        -e '/ExecutionTime = /d' \
        -e '/^[[:space:]]*(Date|Time|Host|PID|Case|nProcs|Build|Exec|Arch)[[:space:]]*:/d' \
        -e '/^[[:space:]]*(opened log at|log written to)/d' \
        -e '/^[[:space:]]*(start|finish) [A-Z][a-z]{2} [A-Z][a-z]{2} /d' \
        -e 's/\(git:[0-9a-f]+\+?\)/(git:<REV>)/g' \
        -e '/Initial memory [0-9]+ kB/d' \
        -e "/^[[:space:]]*removed (directory )?'/d" \
        -e 's/[ \t]+$//' \
        -e 's/[0-9]+\.[0-9]+ s, [0-9]+ kB/<TIME_MEM>/g' \
        -e '/^[[:space:]]*(Requirement already satisfied|Collecting |Installing collected|Successfully installed|Downloading |Using cached |Preparing metadata|Building wheel|Created wheel|Stored in directory|Attempting uninstall|Found existing installation|Uninstalling )/d' \
        -e '/^[[:space:]]*\[notice\]/d' \
        -e '/^Timeout per solution:/d' \
        -e '/^[[:space:]]*(-- |\[ *[0-9]+%\]|make\[|Scanning dependencies|Consolidate compiler)/d' \
        -e '/^[[:space:]]*(gmake|Linking CXX|Building CXX)/d'
}

smoke_clean_case() { # <tutorial dir> <case>
    local d="$TUTORIALS_DIR/$1/$2"
    if [ -x "$d/clean.sh" ]; then
        (cd "$d" && ./clean.sh) > /dev/null 2>&1 || true
    fi
}

# Run one participant against the mock; prints normalized output on stdout and
# returns the participant's exit code.
smoke_run_mock_case() { # <tutorial dir> <case> <timeout>
    local d="$TUTORIALS_DIR/$1/$2" tmo="$3" rc raw
    raw=$(mktemp) || return 1
    timeout "$tmo" bash -c "cd '$d' && LD_PRELOAD='$MOCK_LIB' ./run.sh" > "$raw" 2>&1
    rc=$?
    echo "# exit=$rc"
    smoke_normalize < "$raw"
    rm -f "$raw"
    return $rc
}

# Run one course with the mock via run_course_solutions.sh; normalized output on
# stdout. The course script's exit code is part of the printed "# exit=" header.
smoke_run_mock_course() { # <course number> <timeout>
    local num="$1" tmo="$2" rc raw
    raw=$(mktemp) || return 1
    bash "$SMOKE_LIB_DIR/run_course_solutions.sh" --mock --mock-lib "$MOCK_LIB" \
        --course-dir "$COURSE_DIR" --course-number "$num" --timeout "$tmo" > "$raw" 2>&1
    rc=$?
    echo "# exit=$rc"
    smoke_normalize < "$raw"
    rm -f "$raw"
    return $rc
}
