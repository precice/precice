#!/usr/bin/env python3

import os
import glob
import sys


def file_extension(name):
    _, ext = os.path.splitext(name)
    return ext


def get_cmake_file_paths(root):
    src = os.path.join(root, "src")
    sources = os.path.join(src, "sources.cmake")
    tests = os.path.join(src, "tests.cmake")
    return sources, tests


def is_precice_root(root):
    sources, tests = get_cmake_file_paths(root)
    return os.path.exists(sources) and os.path.exists(tests)


def get_file_lists(root):
    src_dir = os.path.join(root, "src")

    # Find interface headers
    public = glob.glob(os.path.join(src_dir, "precice", "*.hpp"))
    public = [os.path.relpath(p, root) for p in public]

    # Find interface c bindings
    bindings = glob.glob(os.path.join(src_dir, "precice", "bindings", "c", "*.h"))
    bindings = [os.path.relpath(b, root) for b in bindings]

    # Find all test and source cpp files
    sources, tests = [], []
    exts = [".cpp", ".c", ".hpp", ".h"]
    for dir, _, filenames in os.walk(src_dir):
        if "drivers" in dir:
            continue
        files = [
            os.path.relpath(os.path.join(dir, name), root)
            for name in filenames
            if file_extension(name) in exts
        ]
        if "test" in dir:
            tests += files
        else:
            sources += files

    # Remove bindings from the sources
    sources = [source for source in sources if source not in bindings]

    return sources, (public + bindings), tests


SOURCES_BASE = """#
# This file lists all sources that will be compiles into the precice library
#

target_sources(precice
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


def generate_cmake_files(sources, public, tests):
    sources = SOURCES_BASE.format(
        "\n    ".join(sources),
        "\n    ".join(public)
    )
    tests = TESTS_BASE.format(
        "\n    ".join(tests)
    )
    return sources, tests


def main():
    root = os.curdir
    if not is_precice_root(root):
        print("Current dir {} is not the root of the precice repository!".format(root))
        return 1
    sources, public, tests = get_file_lists(root)
    print("Detected:\n sources: {}\n public: {}\n tests: {}".format(len(sources), len(public), len(tests)))

    print("Generating CMake files")
    sources_file, tests_file = get_cmake_file_paths(root)
    sources_content, tests_content = generate_cmake_files(sources, public, tests)

    print("Writing Files")
    print(" {}".format(sources_file))
    with open(sources_file, "w") as f:
        f.write(sources_content)

    print(" {}".format(tests_file))
    with open(tests_file, "w") as f:
        f.write(tests_content)
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
