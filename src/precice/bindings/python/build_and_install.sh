#! /bin/env sh
set -e
echo "Building and installing precice bindings from a scons-build"
if [[ -z "$PRECICE_ROOT" ]]; then
    echo "ERROR: PRECICE_ROOT not set!"
    exit 1
fi
if [[ ! -d "$PRECICE_ROOT/src" ]]; then
    echo "ERROR: no source directory found! \"$PRECICE_ROOT/\""
    exit 1
fi
if [[ ! -d "$PRECICE_ROOT/build/last" ]]; then
    echo "ERROR: no scons build found! \"$PRECICE_ROOT/build/last\""
    exit 1
fi
python3 setup.py build_ext --library-dirs=$PRECICE_ROOT/build/last --include-dirs=$PRECICE_ROOT/src
python3 setup.py install --user
python3 setup.py clean --all
exit 0
