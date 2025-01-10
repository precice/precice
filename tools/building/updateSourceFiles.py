#!/usr/bin/env python3

import collections
import glob
import os
import pathlib
import subprocess
import sys

""" Files matching this pattern will be filtered out """
IGNORE_PATTERNS = ["drivers", "mapping/device"]

""" Configured files, which should be ignored by git, yet installed by CMake"""
CONFIGURED_PUBLIC = ["${PROJECT_BINARY_DIR}/src/precice/Version.h"]

""" Configured files, which should be ignored by git """
CONFIGURED_SOURCES = [
    "${PROJECT_BINARY_DIR}/src/precice/impl/versions.hpp",
    "${CMAKE_BINARY_DIR}/src/precice/impl/versions.cpp",
]


def get_gitfiles():
    ret = subprocess.run(
        ["git", "ls-files", "--full-name"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if ret.returncode != 0:
        return None
    else:
        return ret.stdout.decode().split()


def file_extension(name):
    _, ext = os.path.splitext(name)
    return ext


def get_cmake_file_paths(root):
    sources = os.path.join(root, "src", "sources.cmake")
    utests = os.path.join(root, "src", "tests.cmake")
    itests = os.path.join(root, "tests", "tests.cmake")
    benchmarks = os.path.join(root, "benchmarks", "sources.cmake")
    cmakepaths = collections.namedtuple(
        "CMakePaths", "sources utests itests benchmarks"
    )
    return cmakepaths(sources, utests, itests, benchmarks)


def is_precice_root(root):
    paths = get_cmake_file_paths(root)
    return all(map(os.path.exists, paths))


def get_file_lists(root):
    src_dir = os.path.join(root, "src")
    tests_dir = os.path.join(root, "tests")
    bench_dir = os.path.join(root, "benchmarks")

    # Find interface headers
    public = glob.glob(os.path.join(src_dir, "precice", "*.hpp"))
    public += CONFIGURED_PUBLIC
    public = [os.path.relpath(p, root) for p in public]

    # Find all test and source cpp files
    sources, utests = [], []
    exts = [".cpp", ".c", ".hpp", ".h"]
    for dir, _, filenames in os.walk(src_dir):
        if any([elem in dir for elem in IGNORE_PATTERNS]):
            continue
        files = [
            os.path.relpath(os.path.join(dir, name), root)
            for name in filenames
            if file_extension(name) in exts
        ]
        if "test" in dir:
            utests += files
        else:
            sources += files
    sources += CONFIGURED_SOURCES

    itests = []
    for dir, _, filenames in os.walk(tests_dir):
        if any([elem in dir for elem in IGNORE_PATTERNS]):
            continue
        files = [
            os.path.relpath(os.path.join(dir, name), root)
            for name in filenames
            if file_extension(name) in exts
        ]
        itests += files

    benchmarks = []
    for dir, _, filenames in os.walk(bench_dir):
        files = [
            os.path.relpath(os.path.join(dir, name), root)
            for name in filenames
            if file_extension(name) in exts
        ]
        benchmarks += files

    return (
        sorted(sources),
        sorted(public),
        sorted(utests),
        sorted(itests),
        sorted(benchmarks),
    )


def itest_path_to_suite(path):
    """
    Extracts the test suite from a path and translates it to the boost test name
    """
    dir = pathlib.PurePath(path).parts[1]
    parts = map(lambda s: s.capitalize(), dir.split("-"))
    return "".join(parts)


SOURCES_BASE = """#
# This file lists all sources that will be compiles into the precice library
#

target_sources(preciceCore
    PRIVATE
    {}
    )

#
# Select headers to install
#

set_property(TARGET precice PROPERTY PUBLIC_HEADER
    {}
    )
"""
TESTS_BASE = """#
# This file lists all tests sources that will be compiled into the test executable
#
target_sources(testprecice
    PRIVATE
    {}
    )
"""
ITESTS_BASE = """#
# This file lists all integration test sources and test suites
#
target_sources(testprecice
    PRIVATE
    {}
    )
"""
BENCHMARKS_BASE = """#
# This file lists all benchmarks that will be compiles into precice-bench
#

target_sources(precice-bench
    PRIVATE
    {}
    )
"""


def generate_lib_sources(sources, public):
    return SOURCES_BASE.format("\n    ".join(sources), "\n    ".join(public))


def generate_unit_tests(utests):
    return TESTS_BASE.format("\n    ".join(utests))


def generate_integration_tests(itests):
    return ITESTS_BASE.format("\n    ".join(itests))


def generate_benchmark_sources(sources):
    return BENCHMARKS_BASE.format("\n    ".join(sources))


def main():
    root = os.curdir
    if not is_precice_root(root):
        print("Current dir {} is not the root of the precice repository!".format(root))
        return 1
    sources, public, utests, itests, benchmarks = get_file_lists(root)
    print(
        "Detected files:\n  sources: {}\n  public headers: {}\n  unit tests: {}\n  integration tests: {}\n  benchmarks: {}".format(
            len(sources), len(public), len(utests), len(itests), len(benchmarks)
        )
    )

    gitfiles = get_gitfiles()
    if gitfiles:
        not_tracked = list(
            set(sources + public + utests + itests)
            - set(gitfiles + CONFIGURED_SOURCES + CONFIGURED_PUBLIC)
        )
        if not_tracked:
            print("The source tree contains files not tracked by git.")
            print("Please do one of the following with them:")
            print("  - track them using 'git add'")
            print("  - add them to IGNORE_PATTERNS in this script")
            print("  - add them to CONFIGURED_SOURCES in this script!")
            print("Files:")
            for file in not_tracked:
                print("  {}".format(file))
            print("Verification FAILED")
        else:
            print("Verification SUCCEEDED")
    else:
        print("Git did not run successfully.")
        print("Verification SKIPPED")

    print("Generating CMake files")
    # sources_file, tests_file = get_cmake_file_paths(root)
    files = get_cmake_file_paths(root)
    sources_content = generate_lib_sources(sources, public)
    utests_content = generate_unit_tests(utests)
    itests_content = generate_integration_tests(itests)
    benchmarks_content = generate_benchmark_sources(benchmarks)

    print("Writing Files")
    print(" {}".format(files.sources))
    with open(files.sources, "w") as f:
        f.write(sources_content)

    print(" {}".format(files.utests))
    with open(files.utests, "w") as f:
        f.write(utests_content)

    print(" {}".format(files.itests))
    with open(files.itests, "w") as f:
        f.write(itests_content)

    print(" {}".format(files.benchmarks))
    with open(files.benchmarks, "w") as f:
        f.write(benchmarks_content)

    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
