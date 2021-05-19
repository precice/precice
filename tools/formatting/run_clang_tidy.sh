#!/bin/bash
#
# This script runs clang-tidy on the precice code base
# It generates automatically a new buid directory called "clang-tidy" 
# where the actual script is executed
# ------------------------------------------------------------------
#
# Usage:
# ./formatting/tools/run_clang_tidy.sh $PRECICE_DIR OPTIONAL_CMAKE_ARGS
#   with:
#     PRECICE_DIR points to preCICE source directory
#     OPTIONAL_CMAKE_ARGS are optional arguments to pass to CMake
#
# -----------------------------------------------------------------

# Store the current PRECICE_DIR in SRC
SRC=$1
SRC=$(cd "$SRC";pwd)
shift

# Test the correct utilization of the script: SRC needs to point to the PRECICE_SRC
if test ! -d "$SRC/src" -o ! -d "$SRC/tools" -o ! -d "$SRC/examples" -o ! -f "$SRC/CMakeLists.txt" ; then
    echo "Usage:"
    echo "  run_clang_tidy.sh /path/to/precice"
    exit 1
fi
echo "SRC-DIR=$SRC"

# Specify the correct compile commands: 
# -EXPORT_COMPILE_COMMANDS: is required for clang itself
# -BUILD_TESTING: we don't want to analyze the test targets
# -PRECICE_ENABLE_C: naming conventions are different from the code base
# -PRECICE_ENABLE_FORTRAN: naming conventions are different from the code base
ARGS=("-D" "CMAKE_EXPORT_COMPILE_COMMANDS=ON" "-D" "BUILD_TESTING=OFF" "-D" "BUILD_SHARED_LIBS=ON" "-D" "CMAKE_BUILD_TYPE=Debug" "-D" "PRECICE_ENABLE_C=OFF" "-D" "PRECICE_ENABLE_FORTRAN=OFF" "-D" "PRECICE_Packages=OFF" "$@")


if ! [ -x "$(command -v run-clang-tidy)" ] || ! [ -x "$(command -v clang++)" ]; then
    echo "Unable to detect clang. Make sure run-clang-tidy (part of clang) and clang++ are in the path"
    exit 2
fi

mkdir -p ${SRC}/clang-tidy
cd ${SRC}/clang-tidy
echo `pwd`
CC=clang CXX=clang++ cmake "${ARGS[@]}" "$SRC" || (echo "cmake failed!"; false) || exit 2

# Generate versions.cpp file
cmake --build . --target GitRevision

# pipe output to output.txt apart from the warning summary
run-clang-tidy -p . -quiet -j 2 -header-filter "$SRC/src/precice/impl/*cpp"  2>error.txt >output.txt

# grep interesting errors and make sure we remove duplicates:
egrep '(warning|error): ' output.txt | sort | uniq >clang-tidy.log

# if we have errors, report them and set exit status to failure
if [ -s clang-tidy.log ]; then
    cat clang-tidy.log
    exit 3
fi

echo "All passed"
exit 0

