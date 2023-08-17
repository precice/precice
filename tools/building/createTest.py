#!/usr/bin/env python3

import argparse
import collections
import os
import pathlib
import re


def is_precice_root(dir):
    detect = ["CHANGELOG.md", "CMakeLists.txt", "LICENSE", "src", "tests", "cmake"]
    return all(map(lambda c: os.path.exists(os.path.join(dir, c)), detect))


def find_precice_root():
    search_depth = 10
    current = pathlib.Path(os.path.curdir).absolute()
    candidates = [current] + list(current.parents)[:search_depth]
    for dir in candidates:
        if is_precice_root(dir):
            return dir
    raise BaseException("Unable to find the root directory of precice")


ABBREVIATIONS = ["MPI", "QN", "RBF", "NN", "NP"]


def dirToSuite(dir):
    """
    Takes a kebab-case-directory and transforms it to a CamelCaseDirectory.
    Abbreviations defined above will be all upper case.
    """

    def toSuite(s):
        upper = s.upper()
        if upper in ABBREVIATIONS:
            return upper
        else:
            return s.capitalize()

    return "".join(map(toSuite, dir.split("-")))


def checkTestSuite(arg):
    assert arg
    if " " in arg:
        raise argparse.ArgumentTypeError(
            'The given suite name "{}" cannot contain spaces.'.format(arg)
        )
    if "." in arg:
        raise argparse.ArgumentTypeError(
            'The given suite name "{}" cannot contain the file extensions.'.format(arg)
        )
    if re.search(r"[^a-z-]", arg) is not None:
        raise argparse.ArgumentTypeError(
            'The given suite dir "{}" must be dashed-lower-case.'.format(arg)
        )
    if re.search(r"^[a-z]", arg) is None:
        raise argparse.ArgumentTypeError(
            'The given suite dir "{}" must start with a lowercase letter'.format(arg)
        )
    if re.search(r"[a-z]$", arg) is None:
        raise argparse.ArgumentTypeError(
            'The given suite dir "{}" must end with a lowercase letter'.format(arg)
        )
    return arg


def checkTestName(arg):
    assert arg
    if " " in arg:
        raise argparse.ArgumentTypeError(
            'The given test name "{}" cannot contain spaces.'.format(arg)
        )
    if "." in arg:
        raise argparse.ArgumentTypeError(
            'The given test name "{}" cannot contain the file extensions.'.format(arg)
        )
    if re.search(r"[^a-zA-Z0-9]", arg) is not None:
        raise argparse.ArgumentTypeError(
            'The given test name "{}" must use CamelCase.'.format(arg)
        )
    return arg


def testarg(arg):
    """
    Checks the given test argument and computes:
    - the location as pathlib.PurePath
    - the suites as a list of suite names
    - the name of the test
    """
    parts = pathlib.PurePath(arg).parts
    dirs, name = parts[:-1], parts[-1]
    checkTestName(name)
    [checkTestSuite(d) for d in dirs]

    # If the given path is inside the tests dir, then use the relative path
    full = pathlib.Path(arg).absolute()
    tests = find_precice_root().joinpath("tests")
    try:
        parts = full.relative_to(tests).parts
        dirs, name = parts[:-1], parts[-1]
    except ValueError:
        pass

    location = tests.joinpath(*dirs)
    if location.exists() and not location.is_dir():
        raise argparse.ArgumentTypeError(
            'The given test location "{}" exists, but is not a directory.'.format(
                location
            )
        )

    suites = [dirToSuite(dir) for dir in dirs]
    return collections.namedtuple("TestSetup", "location suites name")(
        location, suites, name
    )


PRECICE_TEST_BODY = """{
  PRECICE_TEST(TODO);

  // Implement your test here.
  BOOST_TEST(false);
  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;

  if (context.isNamed(TODO)) {
    auto meshID = interface.getMeshID(TODO);
    auto dataID = interface.getDataID(TODO, meshID);

    std::vector<double> coords;

    interface.setMeshVertices(meshID, TODO, coords.data(), vertexIDs.data());
  } else {
  }
}
"""


def generateTestSource(name, suite, filepath):
    if os.path.exists(filepath):
        raise BaseException('The test source at "{}" already exists.'.format(filepath))

    includes = ["<precice/precice.hpp>", "<vector>", '"testing/Testing.hpp"']
    suites = ["Integration"] + suite
    space = [""]
    lines = ["#ifndef PRECICE_NO_MPI"]
    lines += space
    lines += ["#include " + inc for inc in includes if inc[0] != "<"]
    lines += space
    lines += ["#include " + inc for inc in includes if inc[0] == "<"]
    lines += space
    lines += ["BOOST_AUTO_TEST_SUITE({})".format(s) for s in suites]
    lines += ["BOOST_AUTO_TEST_CASE({})".format(name), PRECICE_TEST_BODY]
    lines += ["BOOST_AUTO_TEST_SUITE_END() // " + s for s in suites]
    lines += space
    lines += ["#endif // PRECICE_NO_MPI"]

    with open(filepath, "w") as f:
        f.writelines([line + "\n" for line in lines])


def generateTestConfig(name, suite, filepath):
    if os.path.exists(filepath):
        print('The test config at "{}" already exists.'.format(filepath))
    else:
        open(filepath, "w").close()


def main():
    parser = argparse.ArgumentParser(
        description="preCICE integration test creation tool."
    )
    parser.add_argument(
        "test",
        metavar="[Suite/]TestName",
        type=testarg,
        help="The path to the test, the last component being the test name. "
        "If executed within tests/, then the test will be created relative to the local directory. "
        "Otherwise, the path will be assumed to be relative to the tests directory.",
    )
    parser.add_argument(
        "-n", "--dry-run", action="store_true", help="print actions only"
    )
    args = parser.parse_args()

    print("Create directory {}".format(args.test.location))
    if not args.dry_run:
        os.makedirs(args.test.location, exist_ok=True)

    source = args.test.name + ".cpp"
    config = args.test.name + ".xml"
    sourcePath = args.test.location.joinpath(source)
    configPath = args.test.location.joinpath(config)

    print("Create test source {}".format(source))
    if not args.dry_run:
        generateTestSource(args.test.name, args.test.suites, sourcePath)

    print("Create test config {}".format(config))
    if not args.dry_run:
        generateTestConfig(args.test.name, args.test.suites, configPath)
    print("Remember to run tools/building/updateSourceFiles.py or make sourcesIndex")


if __name__ == "__main__":
    main()
